#include <getopt.h>
#include <zlib.h>

#include <htswrapper/bc.h>
#include <htswrapper/bc_scanner.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <sys/stat.h>

namespace {

void print_help(FILE* stream) {
    std::fprintf(stream,
        "Usage: atac_fq_preprocess -1 R1 -2 R2 -3 R3 -o DIR -w FILE [-W FILE]\n"
        "\n"
        "Preprocess 10x Genomics scATAC-seq or Multiome ATAC FASTQs. The R2\n"
        "barcode is matched to the supplied whitelist with up to one substitution,\n"
        "and records without a valid barcode are omitted. The genomic R1 and R3\n"
        "reads are emitted as a paired FASTQ set with the corrected barcode added\n"
        "to each read header as a CB:Z tag.\n"
        "\n"
        "Required arguments:\n"
        "  -1, --r1 FILE          Forward genomic-read FASTQ (plain or gzip)\n"
        "  -2, --r2 FILE          Cell-barcode FASTQ (plain or gzip)\n"
        "  -3, --r3 FILE          Reverse genomic-read FASTQ (plain or gzip)\n"
        "  -o, --output_dir DIR  Directory for the two output FASTQs\n"
        "  -w, --whitelist FILE  scATAC whitelist; for Multiome, the RNA\n"
        "                         whitelist whose barcodes should be reported\n"
        "\n"
        "Multiome option:\n"
        "  -W, --whitelist2 FILE Multiome ATAC whitelist to match in R2. Entries\n"
        "                         must correspond by line to the RNA whitelist\n"
        "                         supplied with -w.\n"
        "\n"
        "Other options:\n"
        "  -h, --help            Display this help and exit\n"
        "\n"
        "Output naming:\n"
        "  The outputs retain the input R1 and R2 basenames under DIR. The second\n"
        "  output contains the reverse genomic sequence from R3 but intentionally\n"
        "  uses the R2 basename expected by the downstream paired-read workflow.\n");
}

int usage_error(const char* message) {
    std::fprintf(stderr, "ERROR: %s\n\n", message);
    print_help(stderr);
    return 1;
}

std::string filename_nopath(const std::string& filename) {
    const std::string::size_type pos = filename.find_last_of('/');
    return pos == std::string::npos ? filename : filename.substr(pos + 1);
}

std::string strip_trailing_slashes(std::string path) {
    while (path.size() > 1 && path[path.size() - 1] == '/') {
        path.erase(path.size() - 1);
    }
    return path;
}

bool ensure_directory(const std::string& path) {
    if (mkdir(path.c_str(), 0775) != 0 && errno != EEXIST) {
        std::fprintf(stderr, "ERROR: unable to create output directory %s: %s\n",
            path.c_str(), std::strerror(errno));
        return false;
    }

    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        std::fprintf(stderr, "ERROR: unable to inspect output path %s: %s\n",
            path.c_str(), std::strerror(errno));
        return false;
    }
    if (!S_ISDIR(info.st_mode)) {
        std::fprintf(stderr, "ERROR: output path is not a directory: %s\n", path.c_str());
        return false;
    }
    return true;
}

// Resolve an existing path, or resolve its existing parent and append the
// final component. The latter supports comparisons before an output exists.
std::string canonical_path_or_destination(const std::string& path) {
    char* resolved = realpath(path.c_str(), NULL);
    if (resolved != NULL) {
        const std::string result(resolved);
        std::free(resolved);
        return result;
    }

    const std::string::size_type pos = path.find_last_of('/');
    const std::string parent = pos == std::string::npos
        ? "."
        : (pos == 0 ? "/" : path.substr(0, pos));
    const std::string basename = pos == std::string::npos
        ? path
        : path.substr(pos + 1);

    if (basename.empty()) {
        return std::string();
    }

    resolved = realpath(parent.c_str(), NULL);
    if (resolved == NULL) {
        return std::string();
    }

    std::string result(resolved);
    std::free(resolved);
    if (result.empty() || result[result.size() - 1] != '/') {
        result += '/';
    }
    return result + basename;
}

bool same_path_or_destination(const std::string& left, const std::string& right) {
    struct stat left_info;
    struct stat right_info;
    if (stat(left.c_str(), &left_info) == 0 &&
        stat(right.c_str(), &right_info) == 0 &&
        left_info.st_dev == right_info.st_dev &&
        left_info.st_ino == right_info.st_ino) {
        return true;
    }

    const std::string left_canonical = canonical_path_or_destination(left);
    const std::string right_canonical = canonical_path_or_destination(right);
    return !left_canonical.empty() &&
        !right_canonical.empty() &&
        left_canonical == right_canonical;
}

void report_gzip_error(gzFile file, const std::string& path, const char* action) {
    int error_code = Z_OK;
    const char* detail = gzerror(file, &error_code);
    if (error_code == Z_ERRNO) {
        detail = std::strerror(errno);
    }
    std::fprintf(stderr, "ERROR: unable to %s %s: %s\n",
        action, path.c_str(), detail == NULL ? "gzip error" : detail);
}

bool write_gzip(gzFile file, const char* data, size_t length, const std::string& path) {
    while (length > 0) {
        const size_t maximum = std::numeric_limits<unsigned int>::max();
        const unsigned int chunk = static_cast<unsigned int>(length > maximum ? maximum : length);
        const int written = gzwrite(file, data, chunk);
        if (written <= 0) {
            report_gzip_error(file, path, "write");
            return false;
        }
        data += written;
        length -= static_cast<size_t>(written);
    }
    return true;
}

bool close_gzip(gzFile& file, const std::string& path) {
    if (file == NULL) {
        return true;
    }
    const int result = gzclose(file);
    file = NULL;
    if (result != Z_OK) {
        std::fprintf(stderr, "ERROR: unable to finish gzip output %s (zlib code %d)\n",
            path.c_str(), result);
        return false;
    }
    return true;
}

bool write_record(gzFile file,
                  const std::string& path,
                  const std::string& header,
                  const char* sequence,
                  size_t sequence_length,
                  const char* quality) {
    return write_gzip(file, header.data(), header.size(), path) &&
        write_gzip(file, sequence, sequence_length, path) &&
        write_gzip(file, "\n+\n", 3, path) &&
        write_gzip(file, quality, sequence_length, path) &&
        write_gzip(file, "\n", 1, path);
}

}  // namespace

int main(int argc, char* argv[]) {
    static const struct option long_options[] = {
        {"r1", required_argument, NULL, '1'},
        {"r2", required_argument, NULL, '2'},
        {"r3", required_argument, NULL, '3'},
        {"output_dir", required_argument, NULL, 'o'},
        {"whitelist", required_argument, NULL, 'w'},
        {"whitelist2", required_argument, NULL, 'W'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    std::string r1_path;
    std::string r2_path;
    std::string r3_path;
    std::string output_directory;
    std::string whitelist_path;
    std::string whitelist2_path;

    int option_index = 0;
    int option = 0;
    if (argc == 1) {
        print_help(stdout);
        return 0;
    }

    while ((option = getopt_long(argc, argv, "1:2:3:o:w:W:h", long_options,
                                  &option_index)) != -1) {
        switch (option) {
            case '1': r1_path = optarg; break;
            case '2': r2_path = optarg; break;
            case '3': r3_path = optarg; break;
            case 'o': output_directory = optarg; break;
            case 'w': whitelist_path = optarg; break;
            case 'W': whitelist2_path = optarg; break;
            case 'h': print_help(stdout); return 0;
            default: return 1;
        }
    }

    if (optind != argc) {
        return usage_error("unexpected positional argument");
    }
    if (r1_path.empty()) {
        return usage_error("--r1 / -1 is required");
    }
    if (r2_path.empty()) {
        return usage_error("--r2 / -2 is required");
    }
    if (r3_path.empty()) {
        return usage_error("--r3 / -3 is required");
    }
    if (output_directory.empty()) {
        return usage_error("--output_dir / -o is required");
    }
    if (whitelist_path.empty()) {
        return usage_error("--whitelist / -w is required");
    }

    output_directory = strip_trailing_slashes(output_directory);
    if (!ensure_directory(output_directory)) {
        return 1;
    }

    const std::string r1_basename = filename_nopath(r1_path);
    const std::string r2_basename = filename_nopath(r2_path);
    if (r1_basename.empty() || r2_basename.empty()) {
        std::fprintf(stderr, "ERROR: input FASTQ paths must end with a filename\n");
        return 1;
    }

    const std::string separator = output_directory == "/" ? "" : "/";
    const std::string output_paths[2] = {
        output_directory + separator + r1_basename,
        output_directory + separator + r2_basename
    };
    const std::string input_paths[3] = {r1_path, r2_path, r3_path};

    for (size_t output_index = 0; output_index < 2; ++output_index) {
        for (size_t input_index = 0; input_index < 3; ++input_index) {
            if (same_path_or_destination(output_paths[output_index], input_paths[input_index])) {
                std::fprintf(stderr,
                    "ERROR: output %s resolves to input %s; choose a different --output_dir\n",
                    output_paths[output_index].c_str(), input_paths[input_index].c_str());
                return 1;
            }
        }
    }
    if (same_path_or_destination(output_paths[0], output_paths[1])) {
        std::fprintf(stderr,
            "ERROR: R1 and R2 would use the same output path %s; input basenames must differ\n",
            output_paths[0].c_str());
        return 1;
    }

    // Initialize inputs and whitelist before opening outputs, so invalid input
    // cannot leave behind newly truncated output files.
    bc_scanner scanner(r1_path, r2_path, r3_path);
    if (whitelist2_path.empty()) {
        scanner.init_10x_ATAC(whitelist_path);
    } else {
        scanner.init_10x_multiome_ATAC(whitelist_path, whitelist2_path);
    }
    scanner.trim_barcodes(true);

    gzFile outputs[2] = {NULL, NULL};
    for (size_t i = 0; i < 2; ++i) {
        outputs[i] = gzopen(output_paths[i].c_str(), "wb");
        if (outputs[i] == NULL) {
            std::fprintf(stderr, "ERROR: unable to open %s for writing: %s\n",
                output_paths[i].c_str(), std::strerror(errno));
            close_gzip(outputs[0], output_paths[0]);
            close_gzip(outputs[1], output_paths[1]);
            return 1;
        }
    }

    bool write_ok = true;
    while (write_ok && scanner.next()) {
        const std::string barcode = bc2str(scanner.barcode);
        std::string header("@");
        header.append(scanner.seq_id, static_cast<size_t>(scanner.seq_id_len));
        header += " CB:Z:";
        header += barcode;
        header += '\n';

        write_ok = write_record(outputs[0], output_paths[0], header,
            scanner.read_f, static_cast<size_t>(scanner.read_f_len), scanner.read_f_qual) &&
            write_record(outputs[1], output_paths[1], header,
            scanner.read_r, static_cast<size_t>(scanner.read_r_len), scanner.read_r_qual);
    }

    const bool close_first_ok = close_gzip(outputs[0], output_paths[0]);
    const bool close_second_ok = close_gzip(outputs[1], output_paths[1]);
    return write_ok && close_first_ok && close_second_ok ? 0 : 1;
}
