#include <getopt.h>
#include <zlib.h>

#include <htslib/kseq.h>

#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>
#include <vector>

KSEQ_INIT(gzFile, gzread)

namespace {

void print_help(FILE* stream) {
    std::fprintf(stream,
        "Usage: split_read_files (-1 R1 -2 R2 [-3 R3] | -s READS) -o DIR -n N\n"
        "\n"
        "Split a single FASTQ or synchronized paired/triplet FASTQs into exactly N\n"
        "gzip-compressed chunks. Records are assigned round-robin, so corresponding\n"
        "records from paired or triplet inputs always go to the same chunk. A chunk\n"
        "can be empty when N is greater than the number of records.\n"
        "\n"
        "Input mode (choose exactly one):\n"
        "  -1, --r1 FILE                 Forward paired-read FASTQ\n"
        "  -2, --r2 FILE                 Reverse paired-read FASTQ\n"
        "  -3, --r3 FILE                 Optional third FASTQ synchronized to R1/R2\n"
        "  -s, --single FILE             Unpaired FASTQ\n"
        "\n"
        "Required arguments:\n"
        "  -o, --output_directory DIR   Directory for split FASTQs\n"
        "  -n, --num_chunks N           Positive number of output chunks per input\n"
        "\n"
        "Other options:\n"
        "  -h, --help                   Display this help and exit\n"
        "\n"
        "Inputs may be plain or gzip-compressed. Output files are named\n"
        "<input-stem>.<1..N>.fastq.gz in DIR.\n");
}

int usage_error(const char* message) {
    std::fprintf(stderr, "ERROR: %s\n\n", message);
    print_help(stderr);
    return 1;
}

std::string strip_trailing_slashes(std::string path) {
    while (path.size() > 1 && path[path.size() - 1] == '/') {
        path.erase(path.size() - 1);
    }
    return path;
}

std::string filename_nopath(const std::string& filename) {
    const std::string::size_type pos = filename.find_last_of('/');
    return pos == std::string::npos ? filename : filename.substr(pos + 1);
}

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
        value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string fastq_stem(const std::string& path) {
    std::string name = filename_nopath(path);
    static const char* const suffixes[] = {
        ".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"
    };
    for (size_t i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
        const std::string suffix(suffixes[i]);
        if (ends_with(name, suffix)) {
            name.erase(name.size() - suffix.size());
            break;
        }
    }
    return name;
}

bool parse_positive_integer(const char* text, int& value) {
    errno = 0;
    char* end = NULL;
    const long parsed = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed <= 0 || parsed > INT_MAX) {
        return false;
    }
    value = static_cast<int>(parsed);
    return true;
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

bool existing_file_identity(const std::string& path, std::pair<dev_t, ino_t>& identity) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        return false;
    }
    identity = std::make_pair(info.st_dev, info.st_ino);
    return true;
}

class GzipFile {
public:
    GzipFile() : file_(NULL) {}

    ~GzipFile() {
        if (file_ != NULL) {
            gzclose(file_);
        }
    }

    bool open(const std::string& path, const char* mode) {
        path_ = path;
        file_ = gzopen(path.c_str(), mode);
        return file_ != NULL;
    }

    gzFile get() const {
        return file_;
    }

    const std::string& path() const {
        return path_;
    }

    bool close() {
        if (file_ == NULL) {
            return true;
        }
        const int result = gzclose(file_);
        file_ = NULL;
        if (result != Z_OK) {
            std::fprintf(stderr,
                "ERROR: unable to finish gzip output %s (zlib code %d)\n",
                path_.c_str(), result);
            return false;
        }
        return true;
    }

private:
    GzipFile(const GzipFile&);
    GzipFile& operator=(const GzipFile&);

    gzFile file_;
    std::string path_;
};

class FastqInput {
public:
    FastqInput() : sequence_(NULL) {}

    ~FastqInput() {
        if (sequence_ != NULL) {
            kseq_destroy(sequence_);
        }
    }

    bool open(const std::string& path) {
        if (!file_.open(path, "rb")) {
            return false;
        }
        sequence_ = kseq_init(file_.get());
        return sequence_ != NULL;
    }

    int read() {
        return kseq_read(sequence_);
    }

    bool clean_eof() const {
        return gzeof(file_.get()) != 0;
    }

    kseq_t* sequence() const {
        return sequence_;
    }

    gzFile file() const {
        return file_.get();
    }

    const std::string& path() const {
        return file_.path();
    }

private:
    FastqInput(const FastqInput&);
    FastqInput& operator=(const FastqInput&);

    GzipFile file_;
    kseq_t* sequence_;
};

void report_gzip_error(gzFile file, const std::string& path, const char* action) {
    int error_code = Z_OK;
    const char* detail = gzerror(file, &error_code);
    if (error_code == Z_ERRNO) {
        detail = std::strerror(errno);
    }
    std::fprintf(stderr, "ERROR: unable to %s %s: %s\n",
        action, path.c_str(), detail == NULL ? "gzip error" : detail);
}

bool write_gzip(GzipFile& output, const char* data, size_t length) {
    while (length > 0) {
        const size_t maximum = std::numeric_limits<unsigned int>::max();
        const unsigned int chunk = static_cast<unsigned int>(length > maximum ? maximum : length);
        const int written = gzwrite(output.get(), data, chunk);
        if (written <= 0) {
            report_gzip_error(output.get(), output.path(), "write");
            return false;
        }
        data += written;
        length -= static_cast<size_t>(written);
    }
    return true;
}

bool write_fastq(const kseq_t* sequence, GzipFile& output) {
    if (sequence->qual.s == NULL || sequence->qual.l != sequence->seq.l) {
        std::fprintf(stderr,
            "ERROR: record %s in an input is not a complete FASTQ record\n",
            sequence->name.s == NULL ? "<unnamed>" : sequence->name.s);
        return false;
    }

    return write_gzip(output, "@", 1) &&
        write_gzip(output, sequence->name.s, sequence->name.l) &&
        (sequence->comment.l == 0 ||
            (write_gzip(output, " ", 1) &&
             write_gzip(output, sequence->comment.s, sequence->comment.l))) &&
        write_gzip(output, "\n", 1) &&
        write_gzip(output, sequence->seq.s, sequence->seq.l) &&
        write_gzip(output, "\n+\n", 3) &&
        write_gzip(output, sequence->qual.s, sequence->qual.l) &&
        write_gzip(output, "\n", 1);
}

void report_read_error(const FastqInput& input, int result) {
    if (result == -2) {
        std::fprintf(stderr, "ERROR: truncated FASTQ quality string in %s\n",
            input.path().c_str());
        return;
    }

    int error_code = Z_OK;
    const char* detail = gzerror(input.file(), &error_code);
    if (error_code == Z_ERRNO) {
        detail = std::strerror(errno);
    }
    if (error_code != Z_OK && error_code != Z_STREAM_END) {
        std::fprintf(stderr, "ERROR: unable to read %s: %s\n",
            input.path().c_str(), detail == NULL ? "gzip error" : detail);
    } else {
        std::fprintf(stderr, "ERROR: malformed FASTQ record in %s (reader code %d)\n",
            input.path().c_str(), result);
    }
}

std::string output_path(const std::string& directory,
                        const std::string& stem,
                        int chunk_number) {
    std::ostringstream path;
    path << directory;
    if (directory != "/") {
        path << '/';
    }
    path << stem << '.' << chunk_number << ".fastq.gz";
    return path.str();
}

}  // namespace

int main(int argc, char* argv[]) {
    static const struct option long_options[] = {
        {"r1", required_argument, NULL, '1'},
        {"r2", required_argument, NULL, '2'},
        {"r3", required_argument, NULL, '3'},
        {"single", required_argument, NULL, 's'},
        {"output_directory", required_argument, NULL, 'o'},
        {"num_chunks", required_argument, NULL, 'n'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    std::string r1_path;
    std::string r2_path;
    std::string r3_path;
    std::string single_path;
    std::string output_directory;
    int number_of_chunks = 0;
    bool chunks_set = false;

    if (argc == 1) {
        print_help(stdout);
        return 0;
    }

    int option_index = 0;
    int option = 0;
    while ((option = getopt_long(argc, argv, "1:2:3:s:o:n:h", long_options,
                                  &option_index)) != -1) {
        switch (option) {
            case '1': r1_path = optarg; break;
            case '2': r2_path = optarg; break;
            case '3': r3_path = optarg; break;
            case 's': single_path = optarg; break;
            case 'o': output_directory = optarg; break;
            case 'n':
                chunks_set = parse_positive_integer(optarg, number_of_chunks);
                if (!chunks_set) {
                    return usage_error("--num_chunks / -n must be a positive integer");
                }
                break;
            case 'h': print_help(stdout); return 0;
            default: return 1;
        }
    }

    if (optind != argc) {
        return usage_error("unexpected positional argument");
    }
    if (!chunks_set) {
        return usage_error("--num_chunks / -n is required");
    }
    if (output_directory.empty()) {
        return usage_error("--output_directory / -o is required");
    }

    const bool any_paired_input = !r1_path.empty() || !r2_path.empty() || !r3_path.empty();
    if (!single_path.empty() && any_paired_input) {
        return usage_error("--single cannot be combined with --r1, --r2, or --r3");
    }
    if (single_path.empty() && (r1_path.empty() || r2_path.empty())) {
        return usage_error("provide either --single, or both --r1 and --r2");
    }
    if (!r3_path.empty() && (r1_path.empty() || r2_path.empty())) {
        return usage_error("--r3 requires both --r1 and --r2");
    }

    std::vector<std::string> input_paths;
    if (!single_path.empty()) {
        input_paths.push_back(single_path);
    } else {
        input_paths.push_back(r1_path);
        input_paths.push_back(r2_path);
        if (!r3_path.empty()) {
            input_paths.push_back(r3_path);
        }
    }

    output_directory = strip_trailing_slashes(output_directory);
    if (!ensure_directory(output_directory)) {
        return 1;
    }

    std::vector<std::string> stems;
    for (size_t i = 0; i < input_paths.size(); ++i) {
        const std::string stem = fastq_stem(input_paths[i]);
        if (stem.empty()) {
            std::fprintf(stderr, "ERROR: unable to derive an output stem from %s\n",
                input_paths[i].c_str());
            return 1;
        }
        stems.push_back(stem);
    }

    std::vector<std::unique_ptr<FastqInput> > inputs;
    for (size_t i = 0; i < input_paths.size(); ++i) {
        std::unique_ptr<FastqInput> input(new FastqInput());
        if (!input->open(input_paths[i])) {
            std::fprintf(stderr, "ERROR: unable to open %s for reading: %s\n",
                input_paths[i].c_str(), std::strerror(errno));
            return 1;
        }
        inputs.push_back(std::move(input));
    }

    std::set<std::string> canonical_inputs;
    std::set<std::pair<dev_t, ino_t> > input_identities;
    for (size_t i = 0; i < input_paths.size(); ++i) {
        const std::string canonical = canonical_path_or_destination(input_paths[i]);
        std::pair<dev_t, ino_t> identity;
        const bool duplicate_path = !canonical.empty() && !canonical_inputs.insert(canonical).second;
        const bool duplicate_file = existing_file_identity(input_paths[i], identity) &&
            !input_identities.insert(identity).second;
        if (duplicate_path || duplicate_file) {
            std::fprintf(stderr, "ERROR: input paths resolve to the same file: %s\n",
                input_paths[i].c_str());
            return 1;
        }
    }

    std::vector<std::string> output_paths;
    std::set<std::string> canonical_outputs;
    std::set<std::pair<dev_t, ino_t> > output_identities;
    for (int chunk = 1; chunk <= number_of_chunks; ++chunk) {
        for (size_t stream = 0; stream < stems.size(); ++stream) {
            const std::string path = output_path(output_directory, stems[stream], chunk);
            const std::string canonical = canonical_path_or_destination(path);
            if (canonical.empty()) {
                std::fprintf(stderr, "ERROR: unable to resolve output destination %s\n",
                    path.c_str());
                return 1;
            }
            if (canonical_inputs.count(canonical) != 0) {
                std::fprintf(stderr, "ERROR: output %s resolves to an input file\n", path.c_str());
                return 1;
            }
            if (!canonical_outputs.insert(canonical).second) {
                std::fprintf(stderr, "ERROR: multiple chunks would use output path %s\n",
                    path.c_str());
                return 1;
            }

            std::pair<dev_t, ino_t> identity;
            if (existing_file_identity(path, identity)) {
                if (input_identities.count(identity) != 0) {
                    std::fprintf(stderr, "ERROR: output %s is a hard link to an input file\n",
                        path.c_str());
                    return 1;
                }
                if (!output_identities.insert(identity).second) {
                    std::fprintf(stderr,
                        "ERROR: multiple output paths refer to the same existing file: %s\n",
                        path.c_str());
                    return 1;
                }
            }
            output_paths.push_back(path);
        }
    }

    std::vector<std::unique_ptr<GzipFile> > outputs;
    for (size_t i = 0; i < output_paths.size(); ++i) {
        std::unique_ptr<GzipFile> output(new GzipFile());
        if (!output->open(output_paths[i], "wb")) {
            std::fprintf(stderr, "ERROR: unable to open %s for writing: %s\n",
                output_paths[i].c_str(), std::strerror(errno));
            return 1;
        }
        outputs.push_back(std::move(output));
    }

    bool processing_ok = true;
    size_t record_number = 0;
    int chunk_index = 0;

    while (processing_ok) {
        const int first_result = inputs[0]->read();
        if (first_result < 0) {
            if (first_result != -1 || !inputs[0]->clean_eof()) {
                report_read_error(*inputs[0], first_result);
                processing_ok = false;
            } else {
                for (size_t stream = 1; stream < inputs.size(); ++stream) {
                    const int extra_result = inputs[stream]->read();
                    if (extra_result >= 0) {
                        std::fprintf(stderr,
                            "ERROR: %s contains more records than %s\n",
                            inputs[stream]->path().c_str(), inputs[0]->path().c_str());
                        processing_ok = false;
                    } else if (extra_result != -1 || !inputs[stream]->clean_eof()) {
                        report_read_error(*inputs[stream], extra_result);
                        processing_ok = false;
                    }
                }
            }
            break;
        }

        for (size_t stream = 1; stream < inputs.size(); ++stream) {
            const int result = inputs[stream]->read();
            if (result < 0) {
                if (result != -1 || !inputs[stream]->clean_eof()) {
                    report_read_error(*inputs[stream], result);
                } else {
                    std::fprintf(stderr,
                        "ERROR: %s ended before %s at record %zu\n",
                        inputs[stream]->path().c_str(), inputs[0]->path().c_str(),
                        record_number + 1);
                }
                processing_ok = false;
                break;
            }
        }
        if (!processing_ok) {
            break;
        }

        for (size_t stream = 0; stream < inputs.size(); ++stream) {
            const size_t output_index = static_cast<size_t>(chunk_index) * inputs.size() + stream;
            if (!write_fastq(inputs[stream]->sequence(), *outputs[output_index])) {
                processing_ok = false;
                break;
            }
        }
        if (!processing_ok) {
            break;
        }

        ++record_number;
        chunk_index = (chunk_index + 1) % number_of_chunks;
    }

    bool close_ok = true;
    for (size_t i = 0; i < outputs.size(); ++i) {
        if (!outputs[i]->close()) {
            close_ok = false;
        }
    }

    return processing_ok && close_ok ? 0 : 1;
}
