// One-pass RNA BAM evidence profiler.
//
// The only BAM traversal is the htslib loop below. STARsolo's ordinary matrix
// contains molecules assigned to one gene, including genomically multimapping
// reads when their union of gene annotations is a singleton. NH therefore is
// not an ordinary-versus-EM discriminator. Ordinary molecules are recovered
// from all mapped alignment records and deduplicated by the
// STARsolo (CB,GX,post-filter corrected UB) identity. Logical countedU reads
// are counted once per (RG,QNAME), before requiring an accepted UB.

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <zlib.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <sys/resource.h>

namespace {

const char* const VERSION = "2.3.0";
const char* const HASH_ALGORITHM = "fnv1a64_seeded_v1";
const uint8_t CURRENT_RAW = 1U << 0;
const uint8_t CURRENT_FILTERED = 1U << 1;
const uint8_t OLD_RAW = 1U << 2;
const uint8_t OLD_FILTERED = 1U << 3;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

struct Options {
    std::string bam;
    std::string output_dir;
    std::string library;
    std::string raw_barcodes;
    std::string raw_matrix;
    std::string filtered_barcodes;
    std::string filtered_matrix;
    std::string features;
    std::string old_raw_barcodes;
    std::string old_filtered_barcodes;
    std::string rg_metadata;
    std::string source_order;
    std::string class_manifest;
    std::string starsolo_feature;
    std::string starsolo_umi_filtering;
    std::string starsolo_umi_dedup;
    std::string starsolo_multimappers;
    uint64_t hash_seed = 1469598103934665603ULL;
    uint32_t hash_bins = 100;
    uint32_t threads = 1;
    uint64_t expected_molecules = 0;
    uint64_t max_memory_bytes = 8ULL * 1024ULL * 1024ULL * 1024ULL;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> result;
    std::string field;
    std::istringstream input(line);
    while (std::getline(input, field, '\t')) result.push_back(field);
    if (!line.empty() && line.back() == '\t') result.push_back("");
    return result;
}

bool file_exists_nonempty(const std::string& path) {
    if (path.empty()) return false;
    std::ifstream handle(path.c_str(), std::ios::binary | std::ios::ate);
    return handle.good() && handle.tellg() > 0;
}

void require_file(const std::string& path, const std::string& label) {
    if (!file_exists_nonempty(path)) fail("missing or empty " + label + ": " + path);
}

uint64_t parse_u64(const std::string& value, const std::string& option) {
    if (value.empty() || value[0] == '-') fail("invalid value for " + option + ": " + value);
    char* end = NULL;
    errno = 0;
    const unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
    if (errno || !end || *end != '\0') fail("invalid value for " + option + ": " + value);
    return static_cast<uint64_t>(parsed);
}

void usage(std::ostream& out) {
    out
        << "rna_bam_evidence " << VERSION << "\n"
        << "Usage: rna_bam_evidence --bam FILE --output-dir DIR --library NAME\n"
        << "  --raw-barcodes FILE --raw-matrix FILE --filtered-barcodes FILE\n"
        << "  --filtered-matrix FILE --features FILE\n"
        << "  --rg-metadata FILE --source-order FILE [--old-raw-barcodes FILE]\n"
        << "  --starsolo-feature GeneFull_Ex50pAS\n"
        << "  --starsolo-umi-filtering MultiGeneUMI_CR\n"
        << "  --starsolo-umi-dedup 1MM_CR --starsolo-multimappers EM\n"
        << "  [--old-filtered-barcodes FILE] [--class-manifest FILE]\n"
        << "  [--threads N] [--hash-bins N] [--hash-seed N]\n"
        << "  [--expected-molecules N] [--max-memory-bytes N]\n";
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--help" || arg == "-h") {
            usage(std::cout);
            std::exit(0);
        }
        if (arg == "--version") {
            std::cout << "rna_bam_evidence " << VERSION << "\n";
            std::exit(0);
        }
        if (i + 1 >= argc) fail("missing value after " + arg);
        const std::string value(argv[++i]);
        if (arg == "--bam") options.bam = value;
        else if (arg == "--output-dir") options.output_dir = value;
        else if (arg == "--library") options.library = value;
        else if (arg == "--raw-barcodes") options.raw_barcodes = value;
        else if (arg == "--raw-matrix") options.raw_matrix = value;
        else if (arg == "--filtered-barcodes") options.filtered_barcodes = value;
        else if (arg == "--filtered-matrix") options.filtered_matrix = value;
        else if (arg == "--features") options.features = value;
        else if (arg == "--old-raw-barcodes") options.old_raw_barcodes = value;
        else if (arg == "--old-filtered-barcodes") options.old_filtered_barcodes = value;
        else if (arg == "--rg-metadata") options.rg_metadata = value;
        else if (arg == "--source-order") options.source_order = value;
        else if (arg == "--class-manifest") options.class_manifest = value;
        else if (arg == "--starsolo-feature") options.starsolo_feature = value;
        else if (arg == "--starsolo-umi-filtering") options.starsolo_umi_filtering = value;
        else if (arg == "--starsolo-umi-dedup") options.starsolo_umi_dedup = value;
        else if (arg == "--starsolo-multimappers") options.starsolo_multimappers = value;
        else if (arg == "--hash-seed") options.hash_seed = parse_u64(value, arg);
        else if (arg == "--hash-bins") options.hash_bins = static_cast<uint32_t>(parse_u64(value, arg));
        else if (arg == "--threads") options.threads = static_cast<uint32_t>(parse_u64(value, arg));
        else if (arg == "--expected-molecules") options.expected_molecules = parse_u64(value, arg);
        else if (arg == "--max-memory-bytes") options.max_memory_bytes = parse_u64(value, arg);
        else fail("unknown option: " + arg);
    }
    if (options.bam.empty() || options.output_dir.empty() || options.library.empty() ||
        options.raw_barcodes.empty() || options.raw_matrix.empty() ||
        options.filtered_barcodes.empty() || options.filtered_matrix.empty() ||
        options.features.empty() || options.rg_metadata.empty() ||
        options.source_order.empty() || options.starsolo_feature.empty() ||
        options.starsolo_umi_filtering.empty() ||
        options.starsolo_umi_dedup.empty() ||
        options.starsolo_multimappers.empty()) {
        usage(std::cerr);
        fail("all required arguments must be supplied");
    }
    if (options.starsolo_feature != "GeneFull_Ex50pAS" ||
        options.starsolo_umi_filtering != "MultiGeneUMI_CR" ||
        options.starsolo_umi_dedup != "1MM_CR" ||
        options.starsolo_multimappers != "EM") {
        fail("this release is semantically bound to STARsolo "
             "GeneFull_Ex50pAS, MultiGeneUMI_CR, 1MM_CR, and EM outputs");
    }
    if (!options.hash_bins || !options.threads) fail("threads and hash bins must be positive");
    if (options.max_memory_bytes < 32ULL * 1024ULL * 1024ULL) {
        fail("--max-memory-bytes must be at least 33554432");
    }
    require_file(options.bam, "BAM");
    require_file(options.raw_barcodes, "raw barcode roster");
    require_file(options.raw_matrix, "ordinary raw matrix");
    require_file(options.filtered_barcodes, "filtered barcode roster");
    require_file(options.filtered_matrix, "ordinary filtered matrix");
    require_file(options.features, "ordinary feature metadata");
    require_file(options.rg_metadata, "RG metadata");
    require_file(options.source_order, "source-order file");
    if (!options.old_raw_barcodes.empty()) require_file(options.old_raw_barcodes, "historical raw barcode roster");
    if (!options.old_filtered_barcodes.empty()) require_file(options.old_filtered_barcodes, "historical filtered barcode roster");
    if (!options.class_manifest.empty()) require_file(options.class_manifest, "classification manifest");
    return options;
}

uint64_t configure_process_memory_limit(uint64_t requested) {
    struct rlimit current;
    if (getrlimit(RLIMIT_AS, &current) != 0) fail("getrlimit(RLIMIT_AS) failed");
    rlim_t effective = static_cast<rlim_t>(requested);
    if (current.rlim_max != RLIM_INFINITY && effective > current.rlim_max) {
        effective = current.rlim_max;
    }
    struct rlimit replacement = current;
    replacement.rlim_cur = effective;
    if (setrlimit(RLIMIT_AS, &replacement) != 0) {
        fail("could not enforce the total-process RLIMIT_AS memory ceiling");
    }
    return static_cast<uint64_t>(effective);
}

class MemoryBudget {
public:
    explicit MemoryBudget(uint64_t process_limit)
        : process_limit_(process_limit),
          admission_limit_(process_limit * 7ULL / 10ULL) {}

    void claim(uint64_t bytes, const std::string& component) {
        if (bytes > admission_limit_ || charged_ > admission_limit_ - bytes) {
            std::ostringstream message;
            message << "total memory admission rejected " << component << " allocation of "
                    << bytes << " bytes: charged " << charged_
                    << ", conservative admission limit " << admission_limit_
                    << ", total process limit " << process_limit_;
            fail(message.str());
        }
        charged_ += bytes;
        peak_ = std::max(peak_, charged_);
        current_[component] += bytes;
        component_peak_[component] =
            std::max(component_peak_[component], current_[component]);
    }

    void release(uint64_t bytes, const std::string& component) {
        if (bytes > charged_ || bytes > current_[component]) {
            fail("internal memory-accounting underflow for " + component);
        }
        charged_ -= bytes;
        current_[component] -= bytes;
    }

    uint64_t process_limit() const { return process_limit_; }
    uint64_t admission_limit() const { return admission_limit_; }
    uint64_t charged() const { return charged_; }
    uint64_t peak() const { return peak_; }
    const std::map<std::string, uint64_t>& component_peak() const {
        return component_peak_;
    }

private:
    uint64_t process_limit_;
    uint64_t admission_limit_;
    uint64_t charged_ = 0;
    uint64_t peak_ = 0;
    std::map<std::string, uint64_t> current_;
    std::map<std::string, uint64_t> component_peak_;
};

uint64_t fnv1a64(const char* data, std::size_t length, uint64_t seed) {
    uint64_t hash = 14695981039346656037ULL ^ seed;
    for (std::size_t i = 0; i < length; ++i) {
        hash ^= static_cast<unsigned char>(data[i]);
        hash *= 1099511628211ULL;
    }
    return hash;
}

uint64_t mix64(uint64_t value) {
    value ^= value >> 30;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27;
    value *= 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

uint64_t stable_qname_hash(const char* qname, uint64_t seed) {
    return fnv1a64(qname, std::strlen(qname), seed);
}

uint64_t stable_qname_hash(const char* qname, std::size_t length, uint64_t seed) {
    return fnv1a64(qname, length, seed);
}

class Interner {
public:
    Interner(MemoryBudget& budget, const std::string& component)
        : budget_(budget), component_(component) {
        ids_.max_load_factor(0.70f);
    }

    void reserve(uint64_t expected) {
        if (expected <= reserved_) return;
        const uint64_t extra = expected - reserved_;
        // Covers vector storage, hash nodes/buckets, both string objects, and
        // a rehash safety margin. RLIMIT_AS remains the definitive process cap.
        budget_.claim(extra * 160ULL, component_);
        try {
            ids_.reserve(static_cast<std::size_t>(expected));
            values_.reserve(static_cast<std::size_t>(expected));
        } catch (...) {
            budget_.release(extra * 160ULL, component_);
            throw;
        }
        reserved_ = expected;
    }

    std::pair<uint32_t, bool> get(const std::string& value) {
        const auto found = ids_.find(value);
        if (found != ids_.end()) return std::make_pair(found->second, false);
        if (values_.size() >= std::numeric_limits<uint32_t>::max()) {
            fail(component_ + " exceeded 32-bit ID space");
        }
        if (values_.size() >= reserved_) reserve(std::max<uint64_t>(1024, reserved_ * 2ULL));
        const uint32_t id = static_cast<uint32_t>(values_.size());
        values_.push_back(value);
        ids_.insert(std::make_pair(values_.back(), id));
        return std::make_pair(id, true);
    }

    bool find(const std::string& value, uint32_t& id) const {
        const auto found = ids_.find(value);
        if (found == ids_.end()) return false;
        id = found->second;
        return true;
    }

    const std::string& value(uint32_t id) const { return values_.at(id); }
    std::size_t size() const { return values_.size(); }

private:
    MemoryBudget& budget_;
    std::string component_;
    uint64_t reserved_ = 0;
    std::vector<std::string> values_;
    std::unordered_map<std::string, uint32_t> ids_;
};

uint64_t roster_record_count(const std::string& path) {
    if (path.empty()) return 0;
    gzFile handle = gzopen(path.c_str(), "rb");
    if (!handle) fail("could not open roster: " + path);
    std::vector<char> buffer(1 << 20);
    uint64_t count = 0;
    while (gzgets(handle, buffer.data(), static_cast<int>(buffer.size())) != NULL) {
        std::string line(buffer.data());
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        if (!line.empty()) ++count;
    }
    if (!gzeof(handle)) {
        gzclose(handle);
        fail("error reading roster: " + path);
    }
    gzclose(handle);
    return count;
}

struct BarcodeRegistry {
    explicit BarcodeRegistry(MemoryBudget& budget)
        : strings(budget, "barcode_registry") {}
    Interner strings;
    std::vector<uint8_t> flags;
};

void load_barcode_roster(const std::string& path, uint8_t flag,
                         BarcodeRegistry& registry,
                         std::vector<uint32_t>* ordered_ids = NULL) {
    if (path.empty()) return;
    gzFile handle = gzopen(path.c_str(), "rb");
    if (!handle) fail("could not open barcode roster: " + path);
    std::vector<char> buffer(1 << 20);
    while (gzgets(handle, buffer.data(), static_cast<int>(buffer.size())) != NULL) {
        std::string line(buffer.data());
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        const std::size_t tab = line.find('\t');
        if (tab != std::string::npos) line.resize(tab);
        if (line.empty()) continue;
        const std::pair<uint32_t, bool> item = registry.strings.get(line);
        if (item.second) registry.flags.push_back(0);
        if (registry.flags.at(item.first) & flag) {
            gzclose(handle);
            fail("duplicate barcode in roster " + path + ": " + line);
        }
        registry.flags[item.first] |= flag;
        if (ordered_ids) ordered_ids->push_back(item.first);
    }
    if (!gzeof(handle)) {
        gzclose(handle);
        fail("error reading barcode roster: " + path);
    }
    gzclose(handle);
}

void load_feature_roster(const std::string& path, Interner& features) {
    gzFile handle = gzopen(path.c_str(), "rb");
    if (!handle) fail("could not open ordinary feature roster: " + path);
    std::vector<char> buffer(1 << 20);
    while (gzgets(handle, buffer.data(), static_cast<int>(buffer.size())) != NULL) {
        std::string line(buffer.data());
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        const std::size_t tab = line.find('\t');
        if (tab != std::string::npos) line.resize(tab);
        if (line.empty()) continue;
        if (!features.get(line).second) {
            gzclose(handle);
            fail("duplicate feature identifier in " + path + ": " + line);
        }
    }
    if (!gzeof(handle)) {
        gzclose(handle);
        fail("error reading ordinary feature roster: " + path);
    }
    gzclose(handle);
    if (!features.size()) fail("ordinary feature roster is empty: " + path);
}

std::vector<std::string> load_source_order(const std::string& path) {
    std::ifstream handle(path.c_str());
    if (!handle) fail("could not open source-order file: " + path);
    std::vector<std::string> order;
    std::unordered_set<std::string> seen;
    std::string line;
    while (std::getline(handle, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        const std::string value = split_tab(line)[0];
        if (value == "source_id" || value == "bp_id") continue;
        if (!seen.insert(value).second) fail("duplicate source in source order: " + value);
        order.push_back(value);
    }
    if (order.empty()) fail("source-order file contains no sources: " + path);
    if (order.size() > 64) fail("source-order file has more than 64 sources");
    return order;
}

struct RGInfo {
    std::string rg_id;
    std::string source_id;
    uint32_t source_index = 0;
    bool in_bam_header = false;
};

struct RGManifest {
    std::vector<RGInfo> rows;
    std::unordered_map<std::string, uint32_t> by_id;
};

RGManifest load_rg_manifest(const std::string& path, const std::string& library,
                            const std::vector<std::string>& source_order) {
    std::ifstream handle(path.c_str());
    if (!handle) fail("could not open RG metadata: " + path);
    std::string line;
    if (!std::getline(handle, line)) fail("empty RG metadata: " + path);
    const std::vector<std::string> header = split_tab(line);
    std::unordered_map<std::string, std::size_t> column;
    for (std::size_t i = 0; i < header.size(); ++i) column[header[i]] = i;
    for (const char* required : {"library", "bp_id", "rg_id"}) {
        if (!column.count(required)) fail("RG metadata lacks column " + std::string(required));
    }
    std::unordered_map<std::string, uint32_t> source_rank;
    for (uint32_t i = 0; i < source_order.size(); ++i) source_rank[source_order[i]] = i;
    std::vector<std::tuple<uint32_t, std::size_t, RGInfo> > pending;
    std::size_t row_number = 0;
    while (std::getline(handle, line)) {
        ++row_number;
        if (line.empty()) continue;
        std::vector<std::string> fields = split_tab(line);
        fields.resize(header.size());
        if (fields[column["library"]] != library) continue;
        RGInfo info;
        info.source_id = fields[column["bp_id"]];
        info.rg_id = fields[column["rg_id"]];
        const auto source = source_rank.find(info.source_id);
        if (source == source_rank.end()) {
            fail("RG " + info.rg_id + " uses a source absent from source order");
        }
        info.source_index = source->second;
        pending.push_back(std::make_tuple(info.source_index, row_number, info));
    }
    std::sort(pending.begin(), pending.end(),
              [](const auto& left, const auto& right) {
                  if (std::get<0>(left) != std::get<0>(right)) {
                      return std::get<0>(left) < std::get<0>(right);
                  }
                  return std::get<2>(left).rg_id < std::get<2>(right).rg_id;
              });
    RGManifest result;
    result.by_id.reserve(pending.size());
    for (const auto& item : pending) {
        const RGInfo info = std::get<2>(item);
        if (info.rg_id.empty()) fail("blank rg_id in RG metadata");
        if (result.by_id.count(info.rg_id)) fail("duplicate rg_id: " + info.rg_id);
        if (result.rows.size() >= 64) fail("library has more than 64 RGs");
        result.by_id[info.rg_id] = static_cast<uint32_t>(result.rows.size());
        result.rows.push_back(info);
    }
    if (result.rows.empty()) fail("no RG metadata rows matched library " + library);
    return result;
}

std::unordered_set<std::string> header_rg_ids(sam_hdr_t* header) {
    std::unordered_set<std::string> result;
    const char* text = sam_hdr_str(header);
    if (!text) return result;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        if (line.compare(0, 4, "@RG\t") != 0) continue;
        for (const std::string& field : split_tab(line)) {
            if (field.compare(0, 3, "ID:") == 0 && field.size() > 3) {
                result.insert(field.substr(3));
            }
        }
    }
    return result;
}

struct ClassInfo {
    std::string species;
    std::string reference_class;
    bool mitochondrial = false;
    bool rrna = false;
};

struct ClassManifest {
    bool supplied = false;
    std::unordered_map<std::string, ClassInfo> contigs;
    std::unordered_map<std::string, ClassInfo> features;
};

bool parse_bool(const std::string& value, const std::string& field) {
    if (value.empty() || value == "0" || value == "false" || value == "False" || value == "no") return false;
    if (value == "1" || value == "true" || value == "True" || value == "yes") return true;
    fail("classification field " + field + " is not boolean: " + value);
}

ClassManifest load_classes(const std::string& path) {
    ClassManifest result;
    if (path.empty()) return result;
    result.supplied = true;
    std::ifstream handle(path.c_str());
    if (!handle) fail("could not open classification manifest: " + path);
    std::string line;
    if (!std::getline(handle, line)) fail("empty classification manifest: " + path);
    const std::vector<std::string> header = split_tab(line);
    std::unordered_map<std::string, std::size_t> column;
    for (std::size_t i = 0; i < header.size(); ++i) column[header[i]] = i;
    const bool normalized = column.count("entity_type") && column.count("identifier");
    const bool legacy = column.count("contig") || column.count("feature_id") || column.count("GX");
    if (!normalized && !legacy) {
        fail("classification manifest requires entity_type+identifier, contig, feature_id, or GX columns");
    }
    while (std::getline(handle, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields = split_tab(line);
        fields.resize(header.size());
        ClassInfo info;
        if (column.count("species")) info.species = fields[column["species"]];
        if (column.count("reference_class")) info.reference_class = fields[column["reference_class"]];
        if (column.count("mitochondrial")) info.mitochondrial = parse_bool(fields[column["mitochondrial"]], "mitochondrial");
        if (column.count("rrna")) info.rrna = parse_bool(fields[column["rrna"]], "rrna");
        std::string entity;
        std::string identifier;
        if (normalized) {
            entity = fields[column["entity_type"]];
            identifier = fields[column["identifier"]];
        } else if (column.count("contig") && !fields[column["contig"]].empty()) {
            entity = "contig";
            identifier = fields[column["contig"]];
        } else {
            entity = "feature";
            if (column.count("feature_id")) {
                identifier = fields[column["feature_id"]];
            } else if (column.count("GX")) {
                identifier = fields[column["GX"]];
            } else {
                fail("blank contig classification identifier");
            }
        }
        if (identifier.empty()) fail("blank classification identifier");
        std::unordered_map<std::string, ClassInfo>* destination = NULL;
        if (entity == "contig") destination = &result.contigs;
        else if (entity == "feature" || entity == "GX" || entity == "feature_id") destination = &result.features;
        else fail("unknown classification entity_type: " + entity);
        if (!destination->insert(std::make_pair(identifier, info)).second) {
            fail("duplicate classification identifier: " + entity + ":" + identifier);
        }
    }
    if (result.contigs.empty() && result.features.empty()) {
        fail("classification manifest has no classification rows");
    }
    return result;
}

enum class AuxStatus { MISSING, VALID, MALFORMED };

AuxStatus aux_string(const bam1_t* record, const char tag[3], std::string& value) {
    value.clear();
    uint8_t* data = bam_aux_get(record, tag);
    if (!data) return AuxStatus::MISSING;
    if (*data != 'Z') return AuxStatus::MALFORMED;
    const char* decoded = bam_aux2Z(data);
    if (!decoded) return AuxStatus::MALFORMED;
    value.assign(decoded);
    return AuxStatus::VALID;
}

bool integer_aux_type(uint8_t type) {
    return type == 'c' || type == 'C' || type == 's' || type == 'S' ||
           type == 'i' || type == 'I';
}

AuxStatus aux_integer(const bam1_t* record, const char tag[3], int64_t& value) {
    uint8_t* data = bam_aux_get(record, tag);
    if (!data) return AuxStatus::MISSING;
    if (!integer_aux_type(*data)) return AuxStatus::MALFORMED;
    value = bam_aux2i(data);
    return AuxStatus::VALID;
}

bool valid_identifier(const std::string& value) {
    return !value.empty() && value != "-" && value != "0";
}

bool unambiguous_gene(const std::string& value) {
    return valid_identifier(value) && value.find(';') == std::string::npos &&
           value.find(',') == std::string::npos;
}

struct Metrics {
    uint64_t all_records = 0;
    uint64_t primary_records = 0;
    uint64_t primary_mapped_reads = 0;
    uint64_t secondary_records = 0;
    uint64_t supplementary_records = 0;
    uint64_t qcfail_records = 0;
    uint64_t bam_duplicate_flag_reads = 0;
    uint64_t nh1_reads = 0;
    uint64_t nh_gt1_reads = 0;
    uint64_t unambiguous_gx_reads = 0;
    uint64_t unique_gene_tagged_reads = 0;
    uint64_t candidate_countedU_reads = 0;
    uint64_t nh_gt1_unique_gene_countedU_reads = 0;
    uint64_t current_cell_associated_reads = 0;
    uint64_t raw_nonfiltered_droplet_reads = 0;
    uint64_t biological_classified_reads = 0;
    uint64_t classification_by_contig_reads = 0;
    uint64_t classification_by_feature_reads = 0;
    uint64_t mitochondrial_reads = 0;
    uint64_t rrna_reads = 0;
    uint64_t mapq_0 = 0, mapq_1_9 = 0, mapq_10_29 = 0, mapq_30_59 = 0;
    uint64_t mapq_60_plus = 0, mapq_255_unavailable = 0;
    int64_t nm_sum = 0;
    uint64_t nm_count = 0;
    int64_t as_sum = 0;
    uint64_t as_count = 0;
    uint64_t aligned_bases = 0, soft_clipped_bases = 0;
    uint64_t cr_equals_cb_reads = 0, cr_differs_cb_reads = 0;
    uint64_t gx_not_feature_manifest_reads = 0;
    uint64_t missing_rg = 0, missing_cb = 0, missing_cr = 0, missing_gx = 0;
    uint64_t missing_gn = 0, missing_ub = 0, missing_ur = 0, missing_nh = 0;
    uint64_t missing_as = 0, missing_nm = 0;
    uint64_t malformed_rg = 0, malformed_cb = 0, malformed_cr = 0, malformed_gx = 0;
    uint64_t malformed_gn = 0, malformed_ub = 0, malformed_ur = 0, malformed_nh = 0;
    uint64_t malformed_as = 0, malformed_nm = 0;
    uint64_t candidate_molecules = 0;
};

struct CoreMetrics {
    uint64_t all_records = 0, primary_records = 0, primary_mapped_reads = 0;
    uint64_t secondary_records = 0, supplementary_records = 0, qcfail_records = 0;
    uint64_t bam_duplicate_flag_reads = 0, nh1_reads = 0, nh_gt1_reads = 0;
    uint64_t unambiguous_gx_reads = 0, unique_gene_tagged_reads = 0;
    uint64_t candidate_countedU_reads = 0, nh_gt1_unique_gene_countedU_reads = 0;
    uint64_t current_cell_associated_reads = 0, raw_nonfiltered_droplet_reads = 0;
    uint64_t biological_classified_reads = 0, classification_by_contig_reads = 0;
    uint64_t classification_by_feature_reads = 0, mitochondrial_reads = 0, rrna_reads = 0;
    uint64_t candidate_molecules = 0;
};

struct Facts {
    bool primary = false;
    bool mapped = false;
    bool nh1_logical_representative = false;
    bool nh1 = false;
    bool nh_multi = false;
    bool gene_valid = false;
    bool gene_in_features = false;
    bool valid_rg = false;
    bool valid_cb = false;
    bool valid_ub = false;
    bool current_raw = false;
    bool current_filtered = false;
    const ClassInfo* class_info = NULL;
    bool class_by_feature = false;
};

struct RecordContribution {
    CoreMetrics core;
    bool has_mapped_detail = false;
    uint8_t mapq = 0;
    bool nm_valid = false;
    int64_t nm_value = 0;
    bool as_valid = false;
    int64_t as_value = 0;
    uint64_t aligned_bases = 0;
    uint64_t soft_clipped_bases = 0;
    // 0 is unavailable, 1 is CR==CB, and 2 is CR!=CB.
    uint8_t barcode_correction_relation = 0;
};

bool is_nh1_logical_representative(const bam1_t* record) {
    const uint16_t flag = record->core.flag;
    if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return false;
    if (!(flag & BAM_FPAIRED)) return true;
    // STAR emits the same whole-fragment uppercase GX on both mapped mates.
    // READ1 represents a mapped pair; READ2 represents the fragment only when
    // READ1 is unmapped and therefore cannot be selected above.
    if (flag & BAM_FREAD1) return true;
    return (flag & BAM_FREAD2) && (flag & BAM_FMUNMAP);
}

template <typename MetricType>
void add_nh_gt1_logical_read(MetricType& metrics, bool counted_u) {
    ++metrics.unique_gene_tagged_reads;
    if (!counted_u) return;
    ++metrics.candidate_countedU_reads;
    ++metrics.nh_gt1_unique_gene_countedU_reads;
}

void record_aux_status(Metrics& metrics, const char* tag, AuxStatus status) {
    if (status == AuxStatus::VALID) return;
    const bool missing = status == AuxStatus::MISSING;
#define TAG_STATUS(name, lower) if (std::strcmp(tag, name) == 0) { if (missing) ++metrics.missing_##lower; else ++metrics.malformed_##lower; return; }
    TAG_STATUS("RG", rg) TAG_STATUS("CB", cb) TAG_STATUS("CR", cr)
    TAG_STATUS("GX", gx) TAG_STATUS("GN", gn) TAG_STATUS("UB", ub)
    TAG_STATUS("UR", ur) TAG_STATUS("NH", nh) TAG_STATUS("AS", as)
    TAG_STATUS("nM", nm)
#undef TAG_STATUS
}

void add_core(CoreMetrics& m, const bam1_t* record, const Facts& f) {
    ++m.all_records;
    if (record->core.flag & BAM_FSECONDARY) ++m.secondary_records;
    if (record->core.flag & BAM_FSUPPLEMENTARY) ++m.supplementary_records;
    if (record->core.flag & BAM_FQCFAIL) ++m.qcfail_records;
    if (record->core.flag & BAM_FDUP) ++m.bam_duplicate_flag_reads;
    if (!f.primary) return;
    ++m.primary_records;
    if (!f.mapped) return;
    ++m.primary_mapped_reads;
    if (f.nh1) ++m.nh1_reads;
    if (f.nh_multi) ++m.nh_gt1_reads;
    if (f.gene_valid) ++m.unambiguous_gx_reads;
    if (f.nh1_logical_representative && f.gene_valid && f.gene_in_features) {
        ++m.unique_gene_tagged_reads;
    }
    if (f.nh1_logical_representative && f.valid_rg && f.valid_cb &&
        f.gene_valid && f.gene_in_features) {
        ++m.candidate_countedU_reads;
    }
    if (f.current_filtered) ++m.current_cell_associated_reads;
    if (f.current_raw && !f.current_filtered) ++m.raw_nonfiltered_droplet_reads;
    if (f.class_info) {
        ++m.biological_classified_reads;
        if (f.class_by_feature) ++m.classification_by_feature_reads;
        else ++m.classification_by_contig_reads;
        if (f.class_info->mitochondrial) ++m.mitochondrial_reads;
        if (f.class_info->rrna) ++m.rrna_reads;
    }
}

void apply_core(CoreMetrics& target, const CoreMetrics& source) {
#define COPY_CORE(field) target.field += source.field
    COPY_CORE(all_records); COPY_CORE(primary_records); COPY_CORE(primary_mapped_reads);
    COPY_CORE(secondary_records); COPY_CORE(supplementary_records); COPY_CORE(qcfail_records);
    COPY_CORE(bam_duplicate_flag_reads); COPY_CORE(nh1_reads); COPY_CORE(nh_gt1_reads);
    COPY_CORE(unambiguous_gx_reads); COPY_CORE(unique_gene_tagged_reads);
    COPY_CORE(candidate_countedU_reads); COPY_CORE(nh_gt1_unique_gene_countedU_reads);
    COPY_CORE(current_cell_associated_reads); COPY_CORE(raw_nonfiltered_droplet_reads);
    COPY_CORE(biological_classified_reads); COPY_CORE(classification_by_contig_reads);
    COPY_CORE(classification_by_feature_reads); COPY_CORE(mitochondrial_reads);
    COPY_CORE(rrna_reads); COPY_CORE(candidate_molecules);
#undef COPY_CORE
}

RecordContribution make_contribution(
    const bam1_t* record, const Facts& facts,
    AuxStatus cb_status, const std::string& cb,
    AuxStatus cr_status, const std::string& cr,
    AuxStatus as_status, int64_t as_value,
    AuxStatus nm_status, int64_t nm_value) {
    RecordContribution result;
    add_core(result.core, record, facts);
    if (!facts.primary || !facts.mapped) return result;
    result.has_mapped_detail = true;
    result.mapq = record->core.qual;
    result.nm_valid = nm_status == AuxStatus::VALID;
    result.nm_value = nm_value;
    result.as_valid = as_status == AuxStatus::VALID;
    result.as_value = as_value;
    const uint32_t* cigar = bam_get_cigar(record);
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        const int operation = bam_cigar_op(cigar[i]);
        const uint64_t length = bam_cigar_oplen(cigar[i]);
        if (operation == BAM_CMATCH || operation == BAM_CEQUAL ||
            operation == BAM_CDIFF) {
            result.aligned_bases += length;
        } else if (operation == BAM_CSOFT_CLIP) {
            result.soft_clipped_bases += length;
        }
    }
    if (cb_status == AuxStatus::VALID && cr_status == AuxStatus::VALID &&
        valid_identifier(cb) && valid_identifier(cr)) {
        result.barcode_correction_relation = cb == cr ? 1 : 2;
    }
    return result;
}

void apply_metrics(Metrics& m, const RecordContribution& contribution) {
    const CoreMetrics& core = contribution.core;
#define COPY_CORE(field) m.field += core.field
    COPY_CORE(all_records); COPY_CORE(primary_records); COPY_CORE(primary_mapped_reads);
    COPY_CORE(secondary_records); COPY_CORE(supplementary_records); COPY_CORE(qcfail_records);
    COPY_CORE(bam_duplicate_flag_reads); COPY_CORE(nh1_reads); COPY_CORE(nh_gt1_reads);
    COPY_CORE(unambiguous_gx_reads); COPY_CORE(unique_gene_tagged_reads);
    COPY_CORE(candidate_countedU_reads); COPY_CORE(nh_gt1_unique_gene_countedU_reads);
    COPY_CORE(current_cell_associated_reads); COPY_CORE(raw_nonfiltered_droplet_reads);
    COPY_CORE(biological_classified_reads); COPY_CORE(classification_by_contig_reads);
    COPY_CORE(classification_by_feature_reads); COPY_CORE(mitochondrial_reads);
    COPY_CORE(rrna_reads);
#undef COPY_CORE
    if (!contribution.has_mapped_detail) return;
    const uint8_t mapq = contribution.mapq;
    if (mapq == 255) ++m.mapq_255_unavailable;
    else if (mapq == 0) ++m.mapq_0;
    else if (mapq < 10) ++m.mapq_1_9;
    else if (mapq < 30) ++m.mapq_10_29;
    else if (mapq < 60) ++m.mapq_30_59;
    else ++m.mapq_60_plus;
    if (contribution.nm_valid) {
        m.nm_sum += contribution.nm_value;
        ++m.nm_count;
    }
    if (contribution.as_valid) {
        m.as_sum += contribution.as_value;
        ++m.as_count;
    }
    m.aligned_bases += contribution.aligned_bases;
    m.soft_clipped_bases += contribution.soft_clipped_bases;
    if (contribution.barcode_correction_relation == 1) ++m.cr_equals_cb_reads;
    else if (contribution.barcode_correction_relation == 2) ++m.cr_differs_cb_reads;
}

std::string classification_status(uint64_t primary_mapped, uint64_t classified) {
    if (!classified) return "unavailable";
    return classified == primary_mapped ? "complete" : "partial";
}

void write_optional_class_counts(std::ostream& out, uint64_t classified,
                                 uint64_t mitochondrial, uint64_t rrna) {
    out << classified << '\t';
    if (!classified) out << "\t\tunavailable";
    else out << mitochondrial << '\t' << rrna << '\t'
             << classification_status(classified, classified);
}

const char* core_metric_header() {
    return "all_records\tprimary_records\tprimary_mapped_reads\tsecondary_records\t"
           "supplementary_records\tqcfail_records\tbam_duplicate_flag_reads\t"
           "nh1_reads\tnh_gt1_reads\tunambiguous_gx_reads\tunique_gene_tagged_reads\t"
           "candidate_countedU_reads\tnh_gt1_unique_gene_countedU_reads\t"
           "current_cell_associated_reads\traw_nonfiltered_droplet_reads\t"
           "biological_classified_reads\tclassification_by_contig_reads\t"
           "classification_by_feature_reads\tmitochondrial_reads\trrna_reads\t"
           "biological_classification_status\tcandidate_matrix_molecules";
}

void write_core(std::ostream& out, const CoreMetrics& m) {
    out << m.all_records << '\t' << m.primary_records << '\t' << m.primary_mapped_reads << '\t'
        << m.secondary_records << '\t' << m.supplementary_records << '\t' << m.qcfail_records << '\t'
        << m.bam_duplicate_flag_reads << '\t' << m.nh1_reads << '\t' << m.nh_gt1_reads << '\t'
        << m.unambiguous_gx_reads << '\t' << m.unique_gene_tagged_reads << '\t'
        << m.candidate_countedU_reads << '\t' << m.nh_gt1_unique_gene_countedU_reads << '\t'
        << m.current_cell_associated_reads << '\t' << m.raw_nonfiltered_droplet_reads << '\t'
        << m.biological_classified_reads << '\t' << m.classification_by_contig_reads << '\t'
        << m.classification_by_feature_reads << '\t';
    if (!m.biological_classified_reads) out << "\t\tunavailable\t";
    else out << m.mitochondrial_reads << '\t' << m.rrna_reads << '\t'
             << classification_status(m.primary_mapped_reads, m.biological_classified_reads) << '\t';
    out << m.candidate_molecules;
}

const char* metric_header() {
    return "all_records\tprimary_records\tprimary_mapped_reads\tsecondary_records\t"
           "supplementary_records\tqcfail_records\tbam_duplicate_flag_reads\t"
           "nh1_reads\tnh_gt1_reads\tunambiguous_gx_reads\tunique_gene_tagged_reads\t"
           "candidate_countedU_reads\tnh_gt1_unique_gene_countedU_reads\t"
           "current_cell_associated_reads\traw_nonfiltered_droplet_reads\t"
           "biological_classified_reads\tclassification_by_contig_reads\t"
           "classification_by_feature_reads\tmitochondrial_reads\trrna_reads\t"
           "biological_classification_status\tmapq_0\tmapq_1_9\tmapq_10_29\t"
           "mapq_30_59\tmapq_60_plus\tmapq_255_unavailable\tnM_sum\tnM_count\t"
           "AS_sum\tAS_count\taligned_bases\tsoft_clipped_bases\t"
           "cr_equals_cb_reads\tcr_differs_cb_reads\tGX_not_in_feature_manifest_reads\t"
           "missing_RG\tmissing_CB\tmissing_CR\tmissing_GX\tmissing_GN\tmissing_UB\t"
           "missing_UR\tmissing_NH\tmissing_AS\tmissing_nM\tmalformed_RG\tmalformed_CB\t"
           "malformed_CR\tmalformed_GX\tmalformed_GN\tmalformed_UB\tmalformed_UR\t"
           "malformed_NH\tmalformed_AS\tmalformed_nM\tcandidate_matrix_molecules";
}

void write_metrics(std::ostream& out, const Metrics& m) {
    out << m.all_records << '\t' << m.primary_records << '\t' << m.primary_mapped_reads << '\t'
        << m.secondary_records << '\t' << m.supplementary_records << '\t' << m.qcfail_records << '\t'
        << m.bam_duplicate_flag_reads << '\t' << m.nh1_reads << '\t' << m.nh_gt1_reads << '\t'
        << m.unambiguous_gx_reads << '\t' << m.unique_gene_tagged_reads << '\t'
        << m.candidate_countedU_reads << '\t' << m.nh_gt1_unique_gene_countedU_reads << '\t'
        << m.current_cell_associated_reads << '\t' << m.raw_nonfiltered_droplet_reads << '\t'
        << m.biological_classified_reads << '\t' << m.classification_by_contig_reads << '\t'
        << m.classification_by_feature_reads << '\t';
    if (!m.biological_classified_reads) out << "\t\tunavailable\t";
    else out << m.mitochondrial_reads << '\t' << m.rrna_reads << '\t'
             << classification_status(m.primary_mapped_reads, m.biological_classified_reads) << '\t';
    out << m.mapq_0 << '\t' << m.mapq_1_9 << '\t' << m.mapq_10_29 << '\t'
        << m.mapq_30_59 << '\t' << m.mapq_60_plus << '\t' << m.mapq_255_unavailable << '\t'
        << m.nm_sum << '\t' << m.nm_count << '\t' << m.as_sum << '\t' << m.as_count << '\t'
        << m.aligned_bases << '\t' << m.soft_clipped_bases << '\t'
        << m.cr_equals_cb_reads << '\t' << m.cr_differs_cb_reads << '\t'
        << m.gx_not_feature_manifest_reads << '\t'
        << m.missing_rg << '\t' << m.missing_cb << '\t' << m.missing_cr << '\t'
        << m.missing_gx << '\t' << m.missing_gn << '\t' << m.missing_ub << '\t'
        << m.missing_ur << '\t' << m.missing_nh << '\t' << m.missing_as << '\t'
        << m.missing_nm << '\t' << m.malformed_rg << '\t' << m.malformed_cb << '\t'
        << m.malformed_cr << '\t' << m.malformed_gx << '\t' << m.malformed_gn << '\t'
        << m.malformed_ub << '\t' << m.malformed_ur << '\t' << m.malformed_nh << '\t'
        << m.malformed_as << '\t' << m.malformed_nm << '\t' << m.candidate_molecules;
}

struct BarcodeRGSlot {
    uint64_t key_plus_one = 0;
    CoreMetrics metrics;
};

class BarcodeRGTable {
public:
    explicit BarcodeRGTable(MemoryBudget& budget) : budget_(budget) {}

    CoreMetrics& get(uint64_t key) {
        if (key == std::numeric_limits<uint64_t>::max()) fail("barcode-RG key overflow");
        if (slots_.empty()) allocate(1024);
        if ((size_ + 1) * 10ULL >= slots_.size() * 7ULL) rehash(slots_.size() * 2ULL);
        std::size_t index = static_cast<std::size_t>(mix64(key) & (slots_.size() - 1));
        while (true) {
            BarcodeRGSlot& slot = slots_[index];
            if (!slot.key_plus_one) {
                slot.key_plus_one = key + 1;
                ++size_;
                return slot.metrics;
            }
            if (slot.key_plus_one == key + 1) return slot.metrics;
            index = (index + 1) & (slots_.size() - 1);
        }
    }

    std::vector<BarcodeRGSlot>& compact() {
        std::size_t destination = 0;
        for (std::size_t i = 0; i < slots_.size(); ++i) {
            if (!slots_[i].key_plus_one) continue;
            if (destination != i) slots_[destination] = slots_[i];
            ++destination;
        }
        slots_.resize(destination);
        return slots_;
    }

    uint64_t size() const { return size_; }
    uint64_t allocated_bytes() const { return slots_.capacity() * sizeof(BarcodeRGSlot); }

private:
    void allocate(uint64_t count) {
        const uint64_t bytes = count * sizeof(BarcodeRGSlot);
        budget_.claim(bytes, "barcode_rg_table");
        try {
            slots_.assign(static_cast<std::size_t>(count), BarcodeRGSlot());
        } catch (...) {
            budget_.release(bytes, "barcode_rg_table");
            throw;
        }
    }

    void rehash(uint64_t count) {
        std::vector<BarcodeRGSlot> old;
        old.swap(slots_);
        const uint64_t old_bytes = old.capacity() * sizeof(BarcodeRGSlot);
        allocate(count);
        size_ = 0;
        for (const BarcodeRGSlot& item : old) {
            if (!item.key_plus_one) continue;
            const uint64_t key = item.key_plus_one - 1;
            std::size_t index = static_cast<std::size_t>(mix64(key) & (slots_.size() - 1));
            while (slots_[index].key_plus_one) index = (index + 1) & (slots_.size() - 1);
            slots_[index] = item;
            ++size_;
        }
        budget_.release(old_bytes, "barcode_rg_table");
    }

    MemoryBudget& budget_;
    uint64_t size_ = 0;
    std::vector<BarcodeRGSlot> slots_;
};

struct MoleculeKey {
    uint32_t barcode;
    uint32_t gene;
    uint64_t umi;
};

struct MoleculeSlot {
    uint32_t barcode_plus_one = 0;
    uint32_t gene_plus_one = 0;
    uint64_t umi = 0;
    uint64_t source_mask = 0;
    uint64_t rg_mask = 0;
    uint64_t min_qname_hash = std::numeric_limits<uint64_t>::max();
};

static_assert(sizeof(MoleculeSlot) == 40, "unexpected MoleculeSlot layout");

uint64_t molecule_hash(const MoleculeKey& key) {
    return mix64(mix64((static_cast<uint64_t>(key.barcode) << 32) | key.gene) ^
                 mix64(key.umi));
}

class MoleculeTable {
public:
    explicit MoleculeTable(MemoryBudget& budget) : budget_(budget) {}

    void reserve(uint64_t expected) {
        uint64_t slots = 1024;
        if (expected > (std::numeric_limits<uint64_t>::max() - 7ULL) / 10ULL) {
            fail("expected molecule count is too large");
        }
        const uint64_t needed = expected ? expected * 10ULL / 7ULL + 1ULL : slots;
        while (slots < needed) {
            if (slots > (std::numeric_limits<uint64_t>::max() >> 1)) fail("molecule table overflow");
            slots <<= 1;
        }
        allocate(slots);
    }

    void update(const MoleculeKey& key, uint64_t source_bit, uint64_t rg_bit,
                uint64_t qname_hash) {
        if (slots_.empty()) reserve(0);
        if ((size_ + 1) * 10ULL >= slots_.size() * 7ULL) rehash(slots_.size() * 2ULL);
        std::size_t index = static_cast<std::size_t>(molecule_hash(key) & (slots_.size() - 1));
        while (true) {
            MoleculeSlot& slot = slots_[index];
            if (!slot.barcode_plus_one) {
                slot.barcode_plus_one = key.barcode + 1;
                slot.gene_plus_one = key.gene + 1;
                slot.umi = key.umi;
                slot.source_mask = source_bit;
                slot.rg_mask = rg_bit;
                slot.min_qname_hash = qname_hash;
                ++size_;
                return;
            }
            if (slot.barcode_plus_one == key.barcode + 1 &&
                slot.gene_plus_one == key.gene + 1 && slot.umi == key.umi) {
                slot.source_mask |= source_bit;
                slot.rg_mask |= rg_bit;
                slot.min_qname_hash = std::min(slot.min_qname_hash, qname_hash);
                return;
            }
            index = (index + 1) & (slots_.size() - 1);
        }
    }

    std::vector<MoleculeSlot>& compact() {
        std::size_t destination = 0;
        for (std::size_t i = 0; i < slots_.size(); ++i) {
            if (!slots_[i].barcode_plus_one) continue;
            if (destination != i) slots_[destination] = slots_[i];
            ++destination;
        }
        slots_.resize(destination);
        return slots_;
    }

    const std::vector<MoleculeSlot>& slots() const { return slots_; }
    uint64_t size() const { return size_; }
    uint64_t allocated_bytes() const { return slots_.capacity() * sizeof(MoleculeSlot); }

private:
    void allocate(uint64_t count) {
        const uint64_t bytes = count * sizeof(MoleculeSlot);
        budget_.claim(bytes, "ordinary_molecule_table");
        try {
            slots_.assign(static_cast<std::size_t>(count), MoleculeSlot());
        } catch (...) {
            budget_.release(bytes, "ordinary_molecule_table");
            throw;
        }
    }

    void rehash(uint64_t count) {
        std::vector<MoleculeSlot> old;
        old.swap(slots_);
        const uint64_t old_bytes = old.capacity() * sizeof(MoleculeSlot);
        allocate(count);
        size_ = 0;
        for (const MoleculeSlot& slot : old) {
            if (!slot.barcode_plus_one) continue;
            MoleculeKey key{slot.barcode_plus_one - 1, slot.gene_plus_one - 1, slot.umi};
            std::size_t index = static_cast<std::size_t>(molecule_hash(key) & (slots_.size() - 1));
            while (slots_[index].barcode_plus_one) index = (index + 1) & (slots_.size() - 1);
            slots_[index] = slot;
            ++size_;
        }
        budget_.release(old_bytes, "ordinary_molecule_table");
    }

    MemoryBudget& budget_;
    uint64_t size_ = 0;
    std::vector<MoleculeSlot> slots_;
};

// NH>1 alignments can place the only uppercase singleton GX on a secondary
// record. Coordinate-sorted BAMs do not keep those records adjacent, so retain
// one exact state per (declared RG, full QNAME). The arena stores every name
// byte for collision-exact comparison; hashes are only probe accelerators.
class QnameArena {
public:
    explicit QnameArena(MemoryBudget& budget) : budget_(budget) {}

    uint64_t append(const char* value, uint32_t length) {
        if (length > std::numeric_limits<uint64_t>::max() - used_) {
            fail("NH>1 QNAME arena size overflow");
        }
        if (blocks_.empty() ||
            static_cast<uint64_t>(blocks_.back().used) + length >
                blocks_.back().capacity) {
            add_block(std::max<uint32_t>(1048576U, length));
        }
        Block& block = blocks_.back();
        const uint32_t offset = block.used;
        std::memcpy(block.data.get() + offset, value, length);
        block.used += length;
        used_ += length;
        return (static_cast<uint64_t>(blocks_.size() - 1) << 32) | offset;
    }

    bool equals(uint64_t location, const char* value, uint32_t length) const {
        const uint32_t block_index = static_cast<uint32_t>(location >> 32);
        const uint32_t offset = static_cast<uint32_t>(location);
        if (block_index >= blocks_.size()) fail("internal NH>1 QNAME block index error");
        const Block& block = blocks_[block_index];
        if (offset > block.used || length > block.used - offset) {
            fail("internal NH>1 QNAME arena range error");
        }
        return std::memcmp(block.data.get() + offset, value, length) == 0;
    }

    uint64_t used_bytes() const { return used_; }
    uint64_t allocated_bytes() const { return allocated_; }

    void clear() {
        const uint64_t descriptor_bytes =
            blocks_.capacity() * static_cast<uint64_t>(sizeof(Block));
        std::vector<Block>().swap(blocks_);
        if (allocated_) budget_.release(allocated_, "nh_gt1_qname_arena");
        if (descriptor_bytes) {
            budget_.release(descriptor_bytes, "nh_gt1_qname_arena_descriptors");
        }
        used_ = 0;
        allocated_ = 0;
    }

private:
    struct Block {
        std::unique_ptr<char[]> data;
        uint32_t used = 0;
        uint32_t capacity = 0;
    };

    void reserve_descriptors(std::size_t requested) {
        if (requested <= blocks_.capacity()) return;
        const std::size_t old_capacity = blocks_.capacity();
        std::size_t new_capacity = std::max<std::size_t>(16, old_capacity * 2);
        while (new_capacity < requested) new_capacity *= 2;
        const uint64_t new_bytes =
            new_capacity * static_cast<uint64_t>(sizeof(Block));
        const uint64_t old_bytes =
            old_capacity * static_cast<uint64_t>(sizeof(Block));
        budget_.claim(new_bytes, "nh_gt1_qname_arena_descriptors");
        try {
            blocks_.reserve(new_capacity);
        } catch (...) {
            budget_.release(new_bytes, "nh_gt1_qname_arena_descriptors");
            throw;
        }
        if (old_bytes) {
            budget_.release(old_bytes, "nh_gt1_qname_arena_descriptors");
        }
    }

    void add_block(uint32_t capacity) {
        if (blocks_.size() >= std::numeric_limits<uint32_t>::max()) {
            fail("NH>1 QNAME arena exceeded 32-bit block space");
        }
        reserve_descriptors(blocks_.size() + 1);
        budget_.claim(capacity, "nh_gt1_qname_arena");
        Block block;
        try {
            block.data.reset(new char[capacity]);
        } catch (...) {
            budget_.release(capacity, "nh_gt1_qname_arena");
            throw;
        }
        block.capacity = capacity;
        allocated_ += capacity;
        blocks_.push_back(std::move(block));
    }

    MemoryBudget& budget_;
    std::vector<Block> blocks_;
    uint64_t used_ = 0;
    uint64_t allocated_ = 0;
};

struct NHMultiReadSlot {
    uint64_t key_hash = 0;
    uint64_t qname_location = 0;
    uint32_t qname_length = 0;
    uint32_t rg_plus_one = 0;
    uint32_t barcode_plus_one = 0;
    uint32_t feature_plus_one = 0;
    uint32_t nh = 0;
    uint32_t padding = 0;
};

static_assert(sizeof(NHMultiReadSlot) == 40, "unexpected NHMultiReadSlot layout");

class NHMultiReadTable {
public:
    explicit NHMultiReadTable(MemoryBudget& budget)
        : budget_(budget), qnames_(budget) {}

    void observe(uint32_t rg, const char* qname, uint32_t qname_length,
                 uint32_t nh, bool valid_cb, uint32_t barcode,
                 bool gene_in_features, uint32_t feature) {
        if (!qname || !qname_length) fail("mapped NH>1 record has an empty QNAME");
        if (slots_.empty()) allocate(1024);
        if ((size_ + 1) * 10ULL >= slots_.size() * 7ULL) {
            rehash(slots_.size() * 2ULL);
        }
        const uint64_t key_hash = mix64(
            stable_qname_hash(qname, qname_length, 0) ^
            mix64(static_cast<uint64_t>(rg) + 0x9e3779b97f4a7c15ULL));
        std::size_t index =
            static_cast<std::size_t>(key_hash & (slots_.size() - 1));
        while (true) {
            NHMultiReadSlot& slot = slots_[index];
            if (!slot.rg_plus_one) {
                slot.key_hash = key_hash;
                slot.qname_location = qnames_.append(qname, qname_length);
                slot.qname_length = qname_length;
                slot.rg_plus_one = rg + 1;
                slot.nh = nh;
                if (valid_cb) slot.barcode_plus_one = barcode + 1;
                if (gene_in_features) slot.feature_plus_one = feature + 1;
                ++size_;
                ++observed_records_;
                return;
            }
            if (slot.key_hash == key_hash && slot.rg_plus_one == rg + 1 &&
                slot.qname_length == qname_length &&
                qnames_.equals(slot.qname_location, qname, qname_length)) {
                if (slot.nh != nh) {
                    fail("NH differs among mapped alignments for one (RG,QNAME)");
                }
                if (valid_cb) {
                    if (slot.barcode_plus_one &&
                        slot.barcode_plus_one != barcode + 1) {
                        fail("corrected CB differs among mapped alignments for one (RG,QNAME)");
                    }
                    slot.barcode_plus_one = barcode + 1;
                }
                if (gene_in_features) {
                    if (slot.feature_plus_one &&
                        slot.feature_plus_one != feature + 1) {
                        fail("uppercase singleton GX differs among mapped alignments for one (RG,QNAME)");
                    }
                    slot.feature_plus_one = feature + 1;
                }
                ++observed_records_;
                return;
            }
            index = (index + 1) & (slots_.size() - 1);
        }
    }

    const std::vector<NHMultiReadSlot>& slots() const { return slots_; }
    uint64_t size() const { return size_; }
    uint64_t observed_records() const { return observed_records_; }
    uint64_t table_allocated_bytes() const {
        return slots_.capacity() * static_cast<uint64_t>(sizeof(NHMultiReadSlot));
    }
    uint64_t qname_used_bytes() const { return qnames_.used_bytes(); }
    uint64_t qname_allocated_bytes() const { return qnames_.allocated_bytes(); }

    void clear() {
        const uint64_t bytes = table_allocated_bytes();
        std::vector<NHMultiReadSlot>().swap(slots_);
        if (bytes) budget_.release(bytes, "nh_gt1_logical_read_table");
        qnames_.clear();
        size_ = 0;
    }

private:
    void allocate(uint64_t count) {
        const uint64_t bytes = count * sizeof(NHMultiReadSlot);
        budget_.claim(bytes, "nh_gt1_logical_read_table");
        try {
            slots_.assign(static_cast<std::size_t>(count), NHMultiReadSlot());
        } catch (...) {
            budget_.release(bytes, "nh_gt1_logical_read_table");
            throw;
        }
    }

    void rehash(uint64_t count) {
        std::vector<NHMultiReadSlot> old;
        old.swap(slots_);
        const uint64_t old_bytes =
            old.capacity() * static_cast<uint64_t>(sizeof(NHMultiReadSlot));
        allocate(count);
        size_ = 0;
        for (const NHMultiReadSlot& item : old) {
            if (!item.rg_plus_one) continue;
            std::size_t index =
                static_cast<std::size_t>(item.key_hash & (slots_.size() - 1));
            while (slots_[index].rg_plus_one) {
                index = (index + 1) & (slots_.size() - 1);
            }
            slots_[index] = item;
            ++size_;
        }
        budget_.release(old_bytes, "nh_gt1_logical_read_table");
    }

    MemoryBudget& budget_;
    QnameArena qnames_;
    uint64_t size_ = 0;
    uint64_t observed_records_ = 0;
    std::vector<NHMultiReadSlot> slots_;
};

bool read_gzip_line(gzFile handle, std::string& line) {
    line.clear();
    char buffer[65536];
    while (true) {
        char* result = gzgets(handle, buffer, static_cast<int>(sizeof(buffer)));
        if (!result) return !line.empty();
        line.append(buffer);
        if (!line.empty() && line.back() == '\n') return true;
        if (gzeof(handle)) return true;
    }
}

struct MatrixValidation {
    uint64_t nonzero_coordinates = 0;
    uint64_t molecules = 0;
};

MatrixValidation validate_matrix_exact(
    const std::string& path,
    const std::string& label,
    uint64_t expected_features,
    uint64_t expected_barcodes,
    const Interner& feature_registry,
    const BarcodeRegistry& barcode_registry,
    const std::vector<uint32_t>* column_barcode_ids,
    const std::vector<uint8_t>& barcode_flags,
    uint8_t membership_flag,
    bool reject_unselected_molecules,
    uint64_t expected_molecules,
    std::vector<MoleculeSlot>& molecules) {
    // STARsolo writes its sparse MatrixMarket coordinates in cell-column,
    // feature-row order. Reordering the already compact molecule table in
    // place permits an exact merge without a raw-barcode vector or hash map.
    std::sort(
        molecules.begin(), molecules.end(),
        [](const MoleculeSlot& left, const MoleculeSlot& right) {
            if (left.barcode_plus_one != right.barcode_plus_one) {
                return left.barcode_plus_one < right.barcode_plus_one;
            }
            if (left.gene_plus_one != right.gene_plus_one) {
                return left.gene_plus_one < right.gene_plus_one;
            }
            return left.umi < right.umi;
        });

    gzFile handle = gzopen(path.c_str(), "rb");
    if (!handle) fail("could not open " + label + ": " + path);
    auto matrix_fail = [&](const std::string& message) {
        gzclose(handle);
        fail(label + " reconciliation failed: " + message);
    };
    std::string line;
    if (!read_gzip_line(handle, line) ||
        line.compare(0, 32, "%%MatrixMarket matrix coordinate") != 0) {
        matrix_fail("unsupported MatrixMarket banner in " + path);
    }
    do {
        if (!read_gzip_line(handle, line)) {
            matrix_fail("missing MatrixMarket dimensions in " + path);
        }
    } while (line.empty() || line[0] == '%');
    uint64_t feature_count = 0, barcode_count = 0, declared_nnz = 0;
    {
        std::istringstream dimensions(line);
        std::string extra;
        if (!(dimensions >> feature_count >> barcode_count >> declared_nnz) ||
            dimensions >> extra) {
            matrix_fail("malformed MatrixMarket dimensions in " + path);
        }
    }
    if (feature_count != expected_features || barcode_count != expected_barcodes) {
        std::ostringstream message;
        message << "dimension mismatch in " << path << ": observed "
                << feature_count << "x" << barcode_count << ", expected "
                << expected_features << "x" << expected_barcodes;
        matrix_fail(message.str());
    }

    MatrixValidation result;
    std::size_t molecule_index = 0;
    uint64_t previous_barcode = 0, previous_feature = 0;
    bool have_previous = false;
    while (read_gzip_line(handle, line)) {
        if (line.empty() || line[0] == '%') continue;
        uint64_t feature = 0, barcode = 0, count = 0;
        std::istringstream fields(line);
        std::string extra;
        if (!(fields >> feature >> barcode >> count) || fields >> extra) {
            matrix_fail("malformed coordinate row in " + path);
        }
        if (!feature || feature > feature_count || !barcode ||
            barcode > barcode_count || !count) {
            matrix_fail("out-of-range or nonpositive coordinate in " + path);
        }
        if (have_previous &&
            (barcode < previous_barcode ||
             (barcode == previous_barcode && feature <= previous_feature))) {
            matrix_fail(
                "coordinates are not strictly ordered by barcode then feature in "
                + path);
        }
        have_previous = true;
        previous_barcode = barcode;
        previous_feature = feature;
        ++result.nonzero_coordinates;
        if (result.molecules > std::numeric_limits<uint64_t>::max() - count) {
            matrix_fail("molecule total overflow in " + path);
        }
        result.molecules += count;

        while (molecule_index < molecules.size()) {
            const uint32_t candidate_barcode =
                molecules[molecule_index].barcode_plus_one - 1;
            if (candidate_barcode < barcode_flags.size() &&
                (barcode_flags[candidate_barcode] & membership_flag)) {
                break;
            }
            if (reject_unselected_molecules) {
                matrix_fail("BAM molecule barcode is absent from the matrix roster");
            }
            ++molecule_index;
        }
        const uint32_t expected_barcode_id = column_barcode_ids
            ? column_barcode_ids->at(static_cast<std::size_t>(barcode - 1))
            : static_cast<uint32_t>(barcode - 1);
        const std::string& expected_feature_name =
            feature_registry.value(static_cast<uint32_t>(feature - 1));
        const std::string& expected_barcode_name =
            barcode_registry.strings.value(expected_barcode_id);
        if (molecule_index >= molecules.size()) {
            std::ostringstream message;
            message << "matrix coordinate has no BAM molecule: feature_index="
                    << feature << ", feature_id=" << expected_feature_name
                    << ", barcode_index=" << barcode << ", CB="
                    << expected_barcode_name << ", matrix_count=" << count
                    << ", bam_count=0";
            matrix_fail(message.str());
        }
        const uint32_t barcode_id = molecules[molecule_index].barcode_plus_one - 1;
        const uint32_t feature_id = molecules[molecule_index].gene_plus_one - 1;
        if (barcode_id != expected_barcode_id || feature_id + 1 != feature) {
            std::ostringstream message;
            message << "matrix and BAM molecule coordinates differ: expected_feature_index="
                    << feature << ", expected_feature_id=" << expected_feature_name
                    << ", expected_barcode_index=" << barcode << ", expected_CB="
                    << expected_barcode_name << ", next_BAM_feature_index="
                    << (feature_id + 1) << ", next_BAM_feature_id="
                    << feature_registry.value(feature_id) << ", next_BAM_CB="
                    << barcode_registry.strings.value(barcode_id);
            matrix_fail(message.str());
        }
        std::size_t end = molecule_index + 1;
        while (end < molecules.size() &&
               molecules[end].barcode_plus_one == barcode_id + 1 &&
               molecules[end].gene_plus_one == feature_id + 1) {
            ++end;
        }
        if (end - molecule_index != count) {
            std::ostringstream message;
            message << "matrix and BAM molecule counts differ at coordinate: feature_index="
                    << feature << ", feature_id=" << expected_feature_name
                    << ", barcode_index=" << barcode << ", CB="
                    << expected_barcode_name << ", matrix_count=" << count
                    << ", bam_count=" << (end - molecule_index);
            matrix_fail(message.str());
        }
        molecule_index = end;
    }
    if (!gzeof(handle)) matrix_fail("compressed matrix read error in " + path);
    gzclose(handle);
    if (result.nonzero_coordinates != declared_nnz) {
        fail(label + " reconciliation failed: nnz declaration mismatch");
    }
    while (molecule_index < molecules.size()) {
        const uint32_t candidate_barcode =
            molecules[molecule_index].barcode_plus_one - 1;
        const bool selected = candidate_barcode < barcode_flags.size() &&
            (barcode_flags[candidate_barcode] & membership_flag);
        if (selected || reject_unselected_molecules) {
            fail(label + " reconciliation failed: BAM has an extra molecule coordinate");
        }
        ++molecule_index;
    }
    if (expected_molecules != std::numeric_limits<uint64_t>::max() &&
        result.molecules != expected_molecules) {
        fail(label + " reconciliation failed: total molecule mismatch");
    }
    return result;
}

struct CorrectionSlot {
    uint32_t rg_plus_one = 0;
    uint32_t raw = 0;
    uint32_t corrected = 0;
    uint32_t padding = 0;
    uint64_t count = 0;
};

uint64_t correction_hash(uint32_t rg, uint32_t raw, uint32_t corrected) {
    return mix64((static_cast<uint64_t>(rg) << 32) ^
                 (static_cast<uint64_t>(raw) << 1) ^ corrected);
}

class CorrectionTable {
public:
    explicit CorrectionTable(MemoryBudget& budget) : budget_(budget) {}

    void increment(uint32_t rg, uint32_t raw, uint32_t corrected) {
        if (slots_.empty()) allocate(1024);
        if ((size_ + 1) * 10ULL >= slots_.size() * 7ULL) rehash(slots_.size() * 2ULL);
        std::size_t index = static_cast<std::size_t>(
            correction_hash(rg, raw, corrected) & (slots_.size() - 1));
        while (true) {
            CorrectionSlot& slot = slots_[index];
            if (!slot.rg_plus_one) {
                slot.rg_plus_one = rg + 1;
                slot.raw = raw;
                slot.corrected = corrected;
                slot.count = 1;
                ++size_;
                return;
            }
            if (slot.rg_plus_one == rg + 1 && slot.raw == raw &&
                slot.corrected == corrected) {
                ++slot.count;
                return;
            }
            index = (index + 1) & (slots_.size() - 1);
        }
    }

    std::vector<CorrectionSlot>& compact() {
        std::size_t destination = 0;
        for (std::size_t i = 0; i < slots_.size(); ++i) {
            if (!slots_[i].rg_plus_one) continue;
            if (destination != i) slots_[destination] = slots_[i];
            ++destination;
        }
        slots_.resize(destination);
        return slots_;
    }

    uint64_t size() const { return size_; }
    uint64_t allocated_bytes() const { return slots_.capacity() * sizeof(CorrectionSlot); }

private:
    void allocate(uint64_t count) {
        const uint64_t bytes = count * sizeof(CorrectionSlot);
        budget_.claim(bytes, "correction_table");
        try {
            slots_.assign(static_cast<std::size_t>(count), CorrectionSlot());
        } catch (...) {
            budget_.release(bytes, "correction_table");
            throw;
        }
    }

    void rehash(uint64_t count) {
        std::vector<CorrectionSlot> old;
        old.swap(slots_);
        const uint64_t old_bytes = old.capacity() * sizeof(CorrectionSlot);
        allocate(count);
        size_ = 0;
        for (const CorrectionSlot& slot : old) {
            if (!slot.rg_plus_one) continue;
            std::size_t index = static_cast<std::size_t>(
                correction_hash(slot.rg_plus_one - 1, slot.raw, slot.corrected) &
                (slots_.size() - 1));
            while (slots_[index].rg_plus_one) index = (index + 1) & (slots_.size() - 1);
            slots_[index] = slot;
            ++size_;
        }
        budget_.release(old_bytes, "correction_table");
    }

    MemoryBudget& budget_;
    uint64_t size_ = 0;
    std::vector<CorrectionSlot> slots_;
};

uint64_t encode_umi(const std::string& umi, Interner& fallback) {
    if (!umi.empty() && umi.size() <= 28) {
        uint64_t encoded = static_cast<uint64_t>(umi.size()) << 56;
        for (std::size_t i = 0; i < umi.size(); ++i) {
            uint64_t base = 0;
            switch (umi[i]) {
                case 'A': base = 0; break;
                case 'C': base = 1; break;
                case 'G': base = 2; break;
                case 'T': base = 3; break;
                default: return (1ULL << 63) | fallback.get(umi).first;
            }
            encoded |= base << (2 * i);
        }
        return encoded;
    }
    return (1ULL << 63) | fallback.get(umi).first;
}

uint32_t nested_bin(uint64_t hash, uint32_t bins) {
    return static_cast<uint32_t>((static_cast<__uint128_t>(hash) * bins) >> 64);
}

std::ofstream output_file(const Options& options, const std::string& name) {
    const std::string path = options.output_dir + "/" + name;
    std::ofstream handle(path.c_str(), std::ios::out | std::ios::trunc);
    if (!handle) fail("could not create output: " + path);
    return handle;
}

std::string barcode_category(uint8_t flags) {
    if (flags & CURRENT_FILTERED) return "current filtered cell";
    if (flags & OLD_FILTERED) return "historical-only filtered cell";
    if (flags & CURRENT_RAW) return "raw-matrix nonfiltered droplet";
    return "corrected barcode observed only in alignment evidence";
}

int profile(const Options& options) {
    const auto started = std::chrono::steady_clock::now();
    const uint64_t process_limit = configure_process_memory_limit(options.max_memory_bytes);
    MemoryBudget budget(process_limit);

    BarcodeRegistry barcodes(budget);
    const uint64_t current_raw_barcode_count =
        roster_record_count(options.raw_barcodes);
    const uint64_t current_filtered_barcode_count =
        roster_record_count(options.filtered_barcodes);
    const uint64_t barcode_reserve =
        current_raw_barcode_count + current_filtered_barcode_count +
        roster_record_count(options.old_raw_barcodes) +
        roster_record_count(options.old_filtered_barcodes) + 1024;
    barcodes.strings.reserve(barcode_reserve);
    barcodes.flags.reserve(static_cast<std::size_t>(barcode_reserve));
    budget.claim(barcode_reserve * sizeof(uint8_t), "barcode_flags");
    load_barcode_roster(options.raw_barcodes, CURRENT_RAW, barcodes);
    if (barcodes.strings.size() != current_raw_barcode_count) {
        fail("ordinary raw barcode roster contains duplicate values");
    }
    std::vector<uint32_t> current_filtered_barcode_ids;
    current_filtered_barcode_ids.reserve(
        static_cast<std::size_t>(current_filtered_barcode_count));
    budget.claim(
        current_filtered_barcode_count * sizeof(uint32_t),
        "filtered_barcode_column_ids");
    load_barcode_roster(
        options.filtered_barcodes,
        CURRENT_FILTERED,
        barcodes,
        &current_filtered_barcode_ids);
    if (barcodes.strings.size() != current_raw_barcode_count) {
        fail("ordinary filtered barcode roster is not a subset of the raw roster");
    }
    if (current_filtered_barcode_ids.size() != current_filtered_barcode_count ||
        !std::is_sorted(
            current_filtered_barcode_ids.begin(),
            current_filtered_barcode_ids.end())) {
        fail("ordinary filtered barcode roster does not preserve raw roster order");
    }
    load_barcode_roster(options.old_raw_barcodes, OLD_RAW, barcodes);
    load_barcode_roster(options.old_filtered_barcodes, OLD_FILTERED, barcodes);

    Interner features(budget, "feature_registry");
    features.reserve(roster_record_count(options.features) + 1);
    load_feature_roster(options.features, features);
    Interner fallback_umis(budget, "nonpacked_umi_interner");
    fallback_umis.reserve(1024);
    Interner raw_barcodes(budget, "raw_barcode_interner");
    raw_barcodes.reserve(1024);

    const std::vector<std::string> source_order = load_source_order(options.source_order);
    budget.claim(
        source_order.size() * 192ULL,
        "source_order_live_storage");
    RGManifest rg_manifest = load_rg_manifest(
        options.rg_metadata, options.library, source_order);
    budget.claim(
        rg_manifest.rows.size() * 512ULL,
        "rg_manifest_live_storage");
    const ClassManifest classes = load_classes(options.class_manifest);
    budget.claim(
        (classes.contigs.size() + classes.features.size()) * 512ULL,
        "classification_manifest_live_storage");

    std::vector<const ClassInfo*> feature_classes(features.size(), NULL);
    budget.claim(feature_classes.size() * sizeof(ClassInfo*), "feature_class_lookup");
    for (uint32_t id = 0; id < features.size(); ++id) {
        const auto found = classes.features.find(features.value(id));
        if (found != classes.features.end()) feature_classes[id] = &found->second;
    }

    samFile* input = sam_open(options.bam.c_str(), "rb");
    if (!input) fail("htslib could not open BAM: " + options.bam);
    if (hts_set_threads(input, static_cast<int>(options.threads)) != 0) {
        sam_close(input);
        fail("htslib could not configure BGZF threads");
    }
    sam_hdr_t* header = sam_hdr_read(input);
    if (!header) {
        sam_close(input);
        fail("htslib could not read BAM header");
    }
    const char* header_text = sam_hdr_str(header);
    budget.claim(
        (header_text ? std::strlen(header_text) : 0ULL) + 4096ULL,
        "bam_header_live_storage_estimate");
    const std::unordered_set<std::string> bam_header_rgs = header_rg_ids(header);
    budget.claim(
        bam_header_rgs.size() * 192ULL,
        "bam_header_rg_lookup_live_storage");
    for (RGInfo& info : rg_manifest.rows) {
        info.in_bam_header = bam_header_rgs.count(info.rg_id) != 0;
    }
    const int32_t target_count = sam_hdr_nref(header);
    if (target_count < 0) fail("htslib reported a negative reference count");
    std::vector<const ClassInfo*> contig_classes(static_cast<std::size_t>(target_count), NULL);
    budget.claim(contig_classes.size() * sizeof(ClassInfo*), "contig_class_lookup");
    for (int32_t tid = 0; tid < target_count; ++tid) {
        const char* name = sam_hdr_tid2name(header, tid);
        if (!name) continue;
        const auto found = classes.contigs.find(name);
        if (found != classes.contigs.end()) contig_classes[tid] = &found->second;
    }

    const uint64_t contig_slots =
        static_cast<uint64_t>(rg_manifest.rows.size()) *
        static_cast<uint64_t>(target_count);
    budget.claim(contig_slots * (sizeof(Metrics) + sizeof(uint8_t)),
                 "dense_rg_contig_accumulators");
    std::vector<Metrics> contig_metrics(static_cast<std::size_t>(contig_slots));
    std::vector<uint8_t> contig_observed(static_cast<std::size_t>(contig_slots), 0);

    MoleculeTable molecules(budget);
    molecules.reserve(options.expected_molecules);
    NHMultiReadTable nh_multi_reads(budget);
    BarcodeRGTable barcode_rg_metrics(budget);
    CorrectionTable corrections(budget);
    // Deque growth avoids the full old+new reallocation peak of a vector for
    // this high-cardinality payload. IDs remain stable array-style indices.
    std::deque<Metrics> observed_metrics;
    std::vector<uint32_t> observed_barcodes;
    std::vector<uint32_t> metric_plus_one(barcodes.strings.size(), 0);
    budget.claim(metric_plus_one.capacity() * sizeof(uint32_t), "barcode_metric_index");
    std::vector<Metrics> rg_metrics(rg_manifest.rows.size());
    Metrics unknown_rg_metrics;
    Metrics library_metrics;
    uint64_t records_rg_not_header = 0;
    uint64_t records_rg_not_manifest = 0;
    uint64_t nh1_logical_representative_records = 0;

    auto ensure_metric_index = [&](uint32_t barcode_id) -> uint32_t {
        if (barcode_id >= metric_plus_one.size()) {
            const std::size_t old_capacity = metric_plus_one.capacity();
            metric_plus_one.resize(barcodes.strings.size(), 0);
            if (metric_plus_one.capacity() > old_capacity) {
                budget.claim(
                    (metric_plus_one.capacity() - old_capacity) * sizeof(uint32_t),
                    "barcode_metric_index");
            }
        }
        if (!metric_plus_one[barcode_id]) {
            budget.claim(sizeof(Metrics) + sizeof(uint32_t) + 32,
                         "observed_barcode_metrics");
            observed_metrics.push_back(Metrics());
            observed_barcodes.push_back(barcode_id);
            metric_plus_one[barcode_id] =
                static_cast<uint32_t>(observed_metrics.size());
        }
        return metric_plus_one[barcode_id] - 1;
    };

    bam1_t* record = bam_init1();
    if (!record) fail("could not allocate BAM record");
    int status = 0;
    std::string rg, cb, cr, gx, gn, ub, ur;
    while ((status = sam_read1(input, header, record)) >= 0) {
        int64_t nh = 0, as_value = 0, nm_value = 0;
        const AuxStatus rg_status = aux_string(record, "RG", rg);
        const AuxStatus cb_status = aux_string(record, "CB", cb);
        const AuxStatus cr_status = aux_string(record, "CR", cr);
        const AuxStatus gx_status = aux_string(record, "GX", gx);
        const AuxStatus gn_status = aux_string(record, "GN", gn);
        const AuxStatus ub_status = aux_string(record, "UB", ub);
        const AuxStatus ur_status = aux_string(record, "UR", ur);
        AuxStatus nh_status = aux_integer(record, "NH", nh);
        const AuxStatus as_status = aux_integer(record, "AS", as_value);
        AuxStatus nm_status = aux_integer(record, "nM", nm_value);
        if (nh_status == AuxStatus::VALID && nh < 1) nh_status = AuxStatus::MALFORMED;
        if (nm_status == AuxStatus::VALID && nm_value < 0) nm_status = AuxStatus::MALFORMED;

        int32_t rg_index = -1;
        if (rg_status == AuxStatus::VALID && valid_identifier(rg)) {
            const auto found = rg_manifest.by_id.find(rg);
            if (found == rg_manifest.by_id.end()) {
                ++records_rg_not_manifest;
                if (!bam_header_rgs.count(rg)) ++records_rg_not_header;
            } else {
                rg_index = static_cast<int32_t>(found->second);
                if (!rg_manifest.rows[rg_index].in_bam_header) ++records_rg_not_header;
            }
        } else if (rg_status == AuxStatus::VALID) {
            ++records_rg_not_manifest;
            ++records_rg_not_header;
        }
        Metrics& rg_target = rg_index >= 0 ? rg_metrics[rg_index] : unknown_rg_metrics;
        for (const std::pair<const char*, AuxStatus>& item : {
                 std::make_pair("RG", rg_status), std::make_pair("CB", cb_status),
                 std::make_pair("CR", cr_status), std::make_pair("GX", gx_status),
                 std::make_pair("GN", gn_status), std::make_pair("UB", ub_status),
                 std::make_pair("UR", ur_status), std::make_pair("NH", nh_status),
                 std::make_pair("AS", as_status), std::make_pair("nM", nm_status)}) {
            record_aux_status(library_metrics, item.first, item.second);
            record_aux_status(rg_target, item.first, item.second);
        }

        Facts facts;
        facts.primary = !(record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
        facts.mapped = !(record->core.flag & BAM_FUNMAP);
        facts.nh1 = nh_status == AuxStatus::VALID && nh == 1;
        facts.nh_multi = nh_status == AuxStatus::VALID && nh > 1;
        facts.nh1_logical_representative =
            facts.nh1 && is_nh1_logical_representative(record);
        if (facts.nh1_logical_representative) {
            ++nh1_logical_representative_records;
        }
        facts.valid_rg = rg_index >= 0;
        facts.valid_cb = cb_status == AuxStatus::VALID && valid_identifier(cb);
        facts.valid_ub = ub_status == AuxStatus::VALID && valid_identifier(ub);
        facts.gene_valid = gx_status == AuxStatus::VALID && unambiguous_gene(gx);

        uint32_t barcode_id = 0;
        uint32_t feature_id = 0;
        if (facts.valid_cb) {
            const std::pair<uint32_t, bool> item = barcodes.strings.get(cb);
            barcode_id = item.first;
            if (item.second) {
                const std::size_t old_capacity = barcodes.flags.capacity();
                barcodes.flags.push_back(0);
                if (barcodes.flags.capacity() > old_capacity) {
                    budget.claim(barcodes.flags.capacity() - old_capacity,
                                 "barcode_flags_growth");
                }
            }
            facts.current_raw = barcodes.flags[barcode_id] & CURRENT_RAW;
            facts.current_filtered = barcodes.flags[barcode_id] & CURRENT_FILTERED;
        }
        if (facts.gene_valid) {
            facts.gene_in_features = features.find(gx, feature_id);
        }
        const char* qname = bam_get_qname(record);
        const std::size_t qname_size = qname ? std::strlen(qname) : 0;
        const bool nonsupplementary_mapped =
            facts.mapped && !(record->core.flag & BAM_FSUPPLEMENTARY);
        if (nonsupplementary_mapped && facts.nh_multi && rg_index >= 0) {
            if (nh > std::numeric_limits<uint32_t>::max() ||
                qname_size > std::numeric_limits<uint32_t>::max()) {
                fail("NH or QNAME length exceeds the compact logical-read state range");
            }
            nh_multi_reads.observe(
                static_cast<uint32_t>(rg_index), qname,
                static_cast<uint32_t>(qname_size), static_cast<uint32_t>(nh),
                facts.valid_cb, barcode_id, facts.gene_in_features, feature_id);
        }
        const ClassInfo* contig_class = NULL;
        if (facts.mapped && record->core.tid >= 0 &&
            record->core.tid < target_count) {
            contig_class = contig_classes[record->core.tid];
        }
        if (facts.gene_in_features && feature_classes[feature_id]) {
            facts.class_info = feature_classes[feature_id];
            facts.class_by_feature = true;
        } else {
            facts.class_info = contig_class;
        }

        const RecordContribution contribution = make_contribution(
            record, facts, cb_status, cb, cr_status, cr,
            as_status, as_value, nm_status, nm_value);

        apply_metrics(library_metrics, contribution);
        apply_metrics(rg_target, contribution);

        uint32_t metric_index = 0;
        if (facts.valid_cb) {
            metric_index = ensure_metric_index(barcode_id);
            Metrics& barcode_metric = observed_metrics[metric_index];
            for (const std::pair<const char*, AuxStatus>& item : {
                     std::make_pair("RG", rg_status), std::make_pair("CB", cb_status),
                     std::make_pair("CR", cr_status), std::make_pair("GX", gx_status),
                     std::make_pair("GN", gn_status), std::make_pair("UB", ub_status),
                     std::make_pair("UR", ur_status), std::make_pair("NH", nh_status),
                     std::make_pair("AS", as_status), std::make_pair("nM", nm_status)}) {
                record_aux_status(barcode_metric, item.first, item.second);
            }
            apply_metrics(barcode_metric, contribution);
            if (rg_index >= 0) {
                const uint64_t key =
                    (static_cast<uint64_t>(barcode_id) << 32) |
                    static_cast<uint32_t>(rg_index);
                apply_core(barcode_rg_metrics.get(key), contribution.core);
            }
        }

        if (facts.primary && facts.mapped && rg_index >= 0 &&
            record->core.tid >= 0 && record->core.tid < target_count) {
            const uint64_t flat =
                static_cast<uint64_t>(rg_index) * target_count + record->core.tid;
            contig_observed[flat] = 1;
            apply_metrics(contig_metrics[flat], contribution);
        }

        const bool ordinary_candidate =
            facts.mapped && facts.valid_cb && facts.valid_ub &&
            facts.gene_valid && facts.gene_in_features;
        if (facts.primary && facts.mapped && facts.gene_valid &&
            !facts.gene_in_features) {
            ++library_metrics.gx_not_feature_manifest_reads;
            ++rg_target.gx_not_feature_manifest_reads;
            if (facts.valid_cb) ++observed_metrics[metric_index].gx_not_feature_manifest_reads;
        }
        if (ordinary_candidate) {
            const uint64_t umi_code = encode_umi(ub, fallback_umis);
            const uint64_t source_bit = rg_index >= 0
                ? 1ULL << rg_manifest.rows[rg_index].source_index
                : 0;
            const uint64_t rg_bit = rg_index >= 0
                ? 1ULL << static_cast<uint32_t>(rg_index)
                : 0;
            molecules.update(
                MoleculeKey{barcode_id, feature_id, umi_code}, source_bit, rg_bit,
                stable_qname_hash(qname, qname_size, options.hash_seed));
        }

        if (facts.primary && facts.mapped && rg_index >= 0 && facts.valid_cb &&
            cr_status == AuxStatus::VALID && valid_identifier(cr)) {
            const uint32_t raw_id = raw_barcodes.get(cr).first;
            corrections.increment(static_cast<uint32_t>(rg_index), raw_id, barcode_id);
        }
    }
    if (status < -1) fail("htslib reported a truncated or corrupt BAM record stream");
    bam_destroy1(record);
    sam_close(input);

    uint64_t nh_gt1_singleton_feature_logical_reads = 0;
    uint64_t nh_gt1_counted_u_logical_reads = 0;
    for (const NHMultiReadSlot& logical_read : nh_multi_reads.slots()) {
        if (!logical_read.rg_plus_one || !logical_read.feature_plus_one) continue;
        ++nh_gt1_singleton_feature_logical_reads;
        const uint32_t rg_index = logical_read.rg_plus_one - 1;
        const bool counted_u = logical_read.barcode_plus_one != 0;
        add_nh_gt1_logical_read(library_metrics, counted_u);
        add_nh_gt1_logical_read(rg_metrics.at(rg_index), counted_u);
        if (!counted_u) continue;
        ++nh_gt1_counted_u_logical_reads;
        const uint32_t barcode_id = logical_read.barcode_plus_one - 1;
        const uint32_t metric_index = metric_plus_one.at(barcode_id) - 1;
        add_nh_gt1_logical_read(observed_metrics.at(metric_index), true);
        const uint64_t key =
            (static_cast<uint64_t>(barcode_id) << 32) | rg_index;
        add_nh_gt1_logical_read(barcode_rg_metrics.get(key), true);
    }
    const uint64_t nh_gt1_logical_read_states = nh_multi_reads.size();
    const uint64_t nh_gt1_alignment_records_observed =
        nh_multi_reads.observed_records();
    const uint64_t nh_gt1_qname_state_table_allocated_bytes =
        nh_multi_reads.table_allocated_bytes();
    const uint64_t nh_gt1_qname_bytes_used = nh_multi_reads.qname_used_bytes();
    const uint64_t nh_gt1_qname_arena_allocated_bytes =
        nh_multi_reads.qname_allocated_bytes();
    // Logical-read state is no longer needed once exact counters have been
    // propagated, so release it before high-cardinality output sorting.
    nh_multi_reads.clear();

    for (const MoleculeSlot& molecule : molecules.slots()) {
        if (!molecule.barcode_plus_one) continue;
        const uint32_t barcode_id = molecule.barcode_plus_one - 1;
        const uint32_t metric_index = metric_plus_one.at(barcode_id) - 1;
        ++observed_metrics[metric_index].candidate_molecules;
        for (uint32_t rg_index = 0; rg_index < rg_manifest.rows.size(); ++rg_index) {
            if (!(molecule.rg_mask & (1ULL << rg_index))) continue;
            ++rg_metrics[rg_index].candidate_molecules;
            const uint64_t key =
                (static_cast<uint64_t>(barcode_id) << 32) | rg_index;
            ++barcode_rg_metrics.get(key).candidate_molecules;
        }
    }
    library_metrics.candidate_molecules = molecules.size();

    budget.claim(observed_barcodes.size() * sizeof(uint32_t),
                 "barcode_output_sort");
    std::sort(observed_barcodes.begin(), observed_barcodes.end(),
              [&](uint32_t left, uint32_t right) {
                  return barcodes.strings.value(left) < barcodes.strings.value(right);
              });
    std::vector<uint32_t> barcode_rank(barcodes.strings.size(),
                                      std::numeric_limits<uint32_t>::max());
    budget.claim(barcode_rank.capacity() * sizeof(uint32_t), "barcode_lexical_rank");
    for (uint32_t rank = 0; rank < observed_barcodes.size(); ++rank) {
        barcode_rank[observed_barcodes[rank]] = rank;
    }

    {
        std::ofstream out = output_file(options, "barcode_read_metrics.tsv");
        out << "CB\tcurrent_raw_member\tcurrent_filtered_member\thistorical_raw_member\t"
               "historical_filtered_member\tbarcode_category\t"
            << metric_header() << '\n';
        for (uint32_t id : observed_barcodes) {
            const uint8_t flags = barcodes.flags[id];
            out << barcodes.strings.value(id) << '\t'
                << static_cast<int>(flags & CURRENT_RAW ? 1 : 0) << '\t'
                << static_cast<int>(flags & CURRENT_FILTERED ? 1 : 0) << '\t'
                << static_cast<int>(flags & OLD_RAW ? 1 : 0) << '\t'
                << static_cast<int>(flags & OLD_FILTERED ? 1 : 0) << '\t'
                << barcode_category(flags) << '\t';
            write_metrics(out, observed_metrics[metric_plus_one[id] - 1]);
            out << '\n';
        }
    }

    std::vector<BarcodeRGSlot>& barcode_rg_rows = barcode_rg_metrics.compact();
    std::sort(barcode_rg_rows.begin(), barcode_rg_rows.end(),
              [&](const BarcodeRGSlot& left, const BarcodeRGSlot& right) {
                  const uint64_t left_key = left.key_plus_one - 1;
                  const uint64_t right_key = right.key_plus_one - 1;
                  const uint32_t left_barcode = static_cast<uint32_t>(left_key >> 32);
                  const uint32_t right_barcode = static_cast<uint32_t>(right_key >> 32);
                  if (barcode_rank[left_barcode] != barcode_rank[right_barcode]) {
                      return barcode_rank[left_barcode] < barcode_rank[right_barcode];
                  }
                  return static_cast<uint32_t>(left_key) <
                         static_cast<uint32_t>(right_key);
              });
    {
        std::ofstream out = output_file(options, "barcode_rg_metrics.tsv");
        out << "CB\tRG\tsource_id\tsource_index\t" << core_metric_header() << '\n';
        for (const BarcodeRGSlot& item : barcode_rg_rows) {
            const uint64_t key = item.key_plus_one - 1;
            const uint32_t barcode_id = static_cast<uint32_t>(key >> 32);
            const uint32_t rg_id = static_cast<uint32_t>(key);
            const RGInfo& rg = rg_manifest.rows.at(rg_id);
            out << barcodes.strings.value(barcode_id) << '\t' << rg.rg_id << '\t'
                << rg.source_id << '\t' << rg.source_index << '\t';
            write_core(out, item.metrics);
            out << '\n';
        }
    }

    std::vector<MoleculeSlot>& molecule_rows = molecules.compact();
    const MatrixValidation raw_matrix_validation = validate_matrix_exact(
        options.raw_matrix,
        "ordinary raw matrix",
        features.size(),
        current_raw_barcode_count,
        features,
        barcodes,
        NULL,
        barcodes.flags,
        CURRENT_RAW,
        true,
        options.expected_molecules,
        molecule_rows);
    const MatrixValidation filtered_matrix_validation = validate_matrix_exact(
        options.filtered_matrix,
        "ordinary filtered matrix",
        features.size(),
        current_filtered_barcode_ids.size(),
        features,
        barcodes,
        &current_filtered_barcode_ids,
        barcodes.flags,
        CURRENT_FILTERED,
        false,
        std::numeric_limits<uint64_t>::max(),
        molecule_rows);
    std::sort(molecule_rows.begin(), molecule_rows.end(),
              [&](const MoleculeSlot& left, const MoleculeSlot& right) {
                  const uint32_t left_barcode = left.barcode_plus_one - 1;
                  const uint32_t right_barcode = right.barcode_plus_one - 1;
                  if (barcode_rank[left_barcode] != barcode_rank[right_barcode]) {
                      return barcode_rank[left_barcode] < barcode_rank[right_barcode];
                  }
                  if (left.source_mask != right.source_mask) {
                      return left.source_mask < right.source_mask;
                  }
                  return nested_bin(left.min_qname_hash, options.hash_bins) <
                         nested_bin(right.min_qname_hash, options.hash_bins);
              });
    uint64_t hash_bin_rows = 0;
    {
        std::ofstream out = output_file(options, "molecule_source_hash_bins.tsv");
        out << "CB\tsource_mask_hex\tnested_min_hash_bin\tn_hash_bins\t"
               "candidate_matrix_molecules\n";
        std::size_t begin = 0;
        while (begin < molecule_rows.size()) {
            const uint32_t barcode_id = molecule_rows[begin].barcode_plus_one - 1;
            const uint64_t mask = molecule_rows[begin].source_mask;
            const uint32_t bin =
                nested_bin(molecule_rows[begin].min_qname_hash, options.hash_bins);
            std::size_t end = begin + 1;
            while (end < molecule_rows.size() &&
                   molecule_rows[end].barcode_plus_one == barcode_id + 1 &&
                   molecule_rows[end].source_mask == mask &&
                   nested_bin(molecule_rows[end].min_qname_hash, options.hash_bins) == bin) {
                ++end;
            }
            out << barcodes.strings.value(barcode_id) << "\t0x" << std::hex
                << std::setw(16) << std::setfill('0') << mask << std::dec
                << std::setfill(' ') << '\t' << bin << '\t' << options.hash_bins
                << '\t' << (end - begin) << '\n';
            ++hash_bin_rows;
            begin = end;
        }
    }

    {
        std::ofstream out = output_file(options, "rg_summary.tsv");
        out << "library\tRG\tsource_id\tsource_index\tRG_in_BAM_header\t"
            << metric_header() << '\n';
        for (std::size_t i = 0; i < rg_manifest.rows.size(); ++i) {
            const RGInfo& rg = rg_manifest.rows[i];
            out << options.library << '\t' << rg.rg_id << '\t' << rg.source_id
                << '\t' << rg.source_index << '\t' << rg.in_bam_header << '\t';
            write_metrics(out, rg_metrics[i]);
            out << '\n';
        }
        if (unknown_rg_metrics.all_records) {
            out << options.library
                << "\t__MISSING_OR_UNDECLARED__\t\t-1\t0\t";
            write_metrics(out, unknown_rg_metrics);
            out << '\n';
        }
    }

    uint64_t observed_contig_rows = 0;
    {
        std::ofstream out = output_file(options, "rg_contig_class_summary.tsv");
        out << "library\tRG\tsource_id\tsource_index\tcontig_index\tcontig\t"
               "contig_class_available\tspecies\tmitochondrial\trrna\t"
               "reference_class\t"
            << metric_header() << '\n';
        for (uint32_t rg_index = 0; rg_index < rg_manifest.rows.size(); ++rg_index) {
            for (int32_t tid = 0; tid < target_count; ++tid) {
                const uint64_t flat = static_cast<uint64_t>(rg_index) * target_count + tid;
                if (!contig_observed[flat]) continue;
                const RGInfo& rg = rg_manifest.rows[rg_index];
                const char* name = sam_hdr_tid2name(header, tid);
                const ClassInfo* info = contig_classes[tid];
                out << options.library << '\t' << rg.rg_id << '\t' << rg.source_id
                    << '\t' << rg.source_index << '\t' << tid << '\t'
                    << (name ? name : "__UNKNOWN_CONTIG__") << '\t'
                    << (info ? "1" : "0") << '\t';
                if (info) {
                    out << info->species << '\t' << info->mitochondrial << '\t'
                        << info->rrna << '\t' << info->reference_class << '\t';
                } else {
                    out << "\t\t\t\t";
                }
                write_metrics(out, contig_metrics[flat]);
                out << '\n';
                ++observed_contig_rows;
            }
        }
    }

    std::vector<CorrectionSlot>& correction_rows = corrections.compact();
    std::sort(correction_rows.begin(), correction_rows.end(),
              [&](const CorrectionSlot& left, const CorrectionSlot& right) {
                  const std::string& left_raw = raw_barcodes.value(left.raw);
                  const std::string& right_raw = raw_barcodes.value(right.raw);
                  if (left_raw != right_raw) return left_raw < right_raw;
                  const std::string& left_rg =
                      rg_manifest.rows[left.rg_plus_one - 1].rg_id;
                  const std::string& right_rg =
                      rg_manifest.rows[right.rg_plus_one - 1].rg_id;
                  if (left_rg != right_rg) return left_rg < right_rg;
                  return barcode_rank[left.corrected] < barcode_rank[right.corrected];
              });
    uint64_t within_rg_conflict_pairs = 0;
    uint64_t cross_rg_conflict_raw_barcodes = 0;
    uint64_t global_conflict_raw_barcodes = 0;
    {
        std::ofstream out = output_file(options, "raw_to_corrected_barcode_counts.tsv");
        out << "CR\tRG\tsource_id\tsource_index\tCB\tread_count\t"
               "within_rg_conflict\tcross_rg_conflict\tglobal_conflict\n";
        std::size_t raw_begin = 0;
        while (raw_begin < correction_rows.size()) {
            std::size_t raw_end = raw_begin + 1;
            while (raw_end < correction_rows.size() &&
                   correction_rows[raw_end].raw == correction_rows[raw_begin].raw) {
                ++raw_end;
            }
            bool global_conflict = false;
            bool cross_rg_conflict = false;
            const uint32_t first_corrected = correction_rows[raw_begin].corrected;
            for (std::size_t i = raw_begin; i < raw_end; ++i) {
                if (correction_rows[i].corrected != first_corrected) {
                    global_conflict = true;
                }
            }
            // Cross-RG conflict describes disagreement between exact RG-level
            // mappings only. An ambiguous RG remains a within-RG/global
            // conflict and is not allowed to contaminate an exact mapping in
            // another RG.
            bool have_exact_rg_mapping = false;
            uint32_t first_exact_rg_corrected = 0;
            std::size_t exact_scan = raw_begin;
            while (exact_scan < raw_end) {
                std::size_t exact_end = exact_scan + 1;
                while (exact_end < raw_end &&
                       correction_rows[exact_end].rg_plus_one ==
                           correction_rows[exact_scan].rg_plus_one) {
                    ++exact_end;
                }
                if (exact_end - exact_scan == 1) {
                    const uint32_t corrected = correction_rows[exact_scan].corrected;
                    if (!have_exact_rg_mapping) {
                        first_exact_rg_corrected = corrected;
                        have_exact_rg_mapping = true;
                    } else if (corrected != first_exact_rg_corrected) {
                        cross_rg_conflict = true;
                    }
                }
                exact_scan = exact_end;
            }
            global_conflict_raw_barcodes += global_conflict;
            cross_rg_conflict_raw_barcodes += cross_rg_conflict;
            std::size_t rg_begin = raw_begin;
            while (rg_begin < raw_end) {
                std::size_t rg_end = rg_begin + 1;
                while (rg_end < raw_end &&
                       correction_rows[rg_end].rg_plus_one ==
                           correction_rows[rg_begin].rg_plus_one) {
                    ++rg_end;
                }
                const bool within_rg_conflict = rg_end - rg_begin > 1;
                within_rg_conflict_pairs += within_rg_conflict;
                for (std::size_t i = rg_begin; i < rg_end; ++i) {
                    const CorrectionSlot& item = correction_rows[i];
                    const RGInfo& rg = rg_manifest.rows[item.rg_plus_one - 1];
                    out << raw_barcodes.value(item.raw) << '\t' << rg.rg_id << '\t'
                        << rg.source_id << '\t' << rg.source_index << '\t'
                        << barcodes.strings.value(item.corrected) << '\t'
                        << item.count << '\t' << within_rg_conflict << '\t'
                        << cross_rg_conflict << '\t' << global_conflict << '\n';
                }
                rg_begin = rg_end;
            }
            raw_begin = raw_end;
        }
    }

    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) fail("getrusage failed");
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    {
        std::ofstream out = output_file(options, "profiler_summary.tsv");
        out << "metric\tvalue\n"
            << "profiler_version\t" << VERSION << '\n'
            << "cplusplus_standard\t" << __cplusplus << '\n'
#ifdef __VERSION__
            << "compiler_version\t" << __VERSION__ << '\n'
#endif
            << "htslib_version\t" << hts_version() << '\n'
            << "zlib_version\t" << zlibVersion() << '\n'
            << "hash_algorithm\t" << HASH_ALGORITHM << '\n'
            << "hash_seed\t" << options.hash_seed << '\n'
            << "hash_bins\t" << options.hash_bins << '\n'
            << "starsolo_feature\t" << options.starsolo_feature << '\n'
            << "starsolo_umi_filtering\t" << options.starsolo_umi_filtering << '\n'
            << "starsolo_umi_dedup\t" << options.starsolo_umi_dedup << '\n'
            << "starsolo_multimappers\t" << options.starsolo_multimappers << '\n'
            << "ordinary_countedU_read_definition\tdistinct_RG_QNAME_declared_RG_valid_CB_singleton_feature_GX_NH_unrestricted_UB_not_required_STARsolo_pre_UMI_filter_countedU\n"
            << "ordinary_matrix_molecule_definition\tdistinct_CB_GX_valid_STARsolo_corrected_UB_after_1MM_CR_and_MultiGeneUMI_CR_NH_unrestricted\n"
            << "summary_unique_read_metric\tUnique Reads in Cells Mapped to GeneFull_Ex50pAS\n"
            << "ordinary_matrix_reconciliation_target\tordinary_STARsolo_raw_and_filtered_matrix.mtx_and_Summary.csv\n"
            << "ordinary_raw_matrix_exact_coordinate_reconciliation\tPASS\n"
            << "ordinary_raw_matrix_nonzero_coordinates\t"
            << raw_matrix_validation.nonzero_coordinates << '\n'
            << "ordinary_raw_matrix_molecules\t"
            << raw_matrix_validation.molecules << '\n'
            << "ordinary_filtered_matrix_exact_coordinate_reconciliation\tPASS\n"
            << "ordinary_filtered_matrix_nonzero_coordinates\t"
            << filtered_matrix_validation.nonzero_coordinates << '\n'
            << "ordinary_filtered_matrix_molecules\t"
            << filtered_matrix_validation.molecules << '\n'
            << "nh_gt1_unique_gene_countedU_definition\tsubset_of_ordinary_countedU_reads_with_NH_gt1_and_singleton_feature_GX\n"
            << "starsolo_EM_evidence_definition\tSTARsolo_multi_gene_EM_unavailable_from_standard_uppercase_GX_UB_BAM_tags_NH_is_not_EM_membership\n"
            << "starsolo_EM_evidence_availability\tunavailable_from_standard_uppercase_GX_UB_BAM_tags\n"
            << "header_rg_count\t" << bam_header_rgs.size() << '\n'
            << "manifest_rg_count\t" << rg_manifest.rows.size() << '\n'
            << "source_order_count\t" << source_order.size() << '\n'
            << "records_rg_not_header\t" << records_rg_not_header << '\n'
            << "records_rg_not_manifest\t" << records_rg_not_manifest << '\n'
            << "total_records\t" << library_metrics.all_records << '\n'
            << "primary_mapped_reads\t" << library_metrics.primary_mapped_reads << '\n'
            << "secondary_records\t" << library_metrics.secondary_records << '\n'
            << "supplementary_records\t" << library_metrics.supplementary_records << '\n'
            << "qcfail_records\t" << library_metrics.qcfail_records << '\n'
            << "bam_duplicate_flag_reads\t" << library_metrics.bam_duplicate_flag_reads << '\n'
            << "candidate_countedU_reads\t" << library_metrics.candidate_countedU_reads << '\n'
            << "nh_gt1_unique_gene_countedU_reads\t"
            << library_metrics.nh_gt1_unique_gene_countedU_reads << '\n'
            << "nh1_logical_representative_records\t"
            << nh1_logical_representative_records << '\n'
            << "nh1_logical_representative_definition\tunpaired_or_READ1_or_mapped_READ2_with_unmapped_mate\n"
            << "nh_gt1_logical_read_key\tdeclared_RG_plus_full_QNAME_collision_exact\n"
            << "nh_gt1_logical_read_state_method\tmemory_budget_accounted_open_addressed_slots_plus_full_qname_byte_arena\n"
            << "nh_gt1_logical_read_state_released_before_output_sort\t1\n"
            << "nh_gt1_logical_read_states\t" << nh_gt1_logical_read_states << '\n'
            << "nh_gt1_alignment_records_observed\t"
            << nh_gt1_alignment_records_observed << '\n'
            << "nh_gt1_singleton_feature_logical_reads\t"
            << nh_gt1_singleton_feature_logical_reads << '\n'
            << "nh_gt1_countedU_logical_reads\t"
            << nh_gt1_counted_u_logical_reads << '\n'
            << "candidate_matrix_molecules\t" << molecules.size() << '\n'
            << "missing_RG\t" << library_metrics.missing_rg << '\n'
            << "missing_CB\t" << library_metrics.missing_cb << '\n'
            << "missing_CR\t" << library_metrics.missing_cr << '\n'
            << "missing_GX\t" << library_metrics.missing_gx << '\n'
            << "missing_GN\t" << library_metrics.missing_gn << '\n'
            << "missing_UB\t" << library_metrics.missing_ub << '\n'
            << "missing_UR\t" << library_metrics.missing_ur << '\n'
            << "missing_NH\t" << library_metrics.missing_nh << '\n'
            << "missing_AS\t" << library_metrics.missing_as << '\n'
            << "missing_nM\t" << library_metrics.missing_nm << '\n'
            << "malformed_RG\t" << library_metrics.malformed_rg << '\n'
            << "malformed_CB\t" << library_metrics.malformed_cb << '\n'
            << "malformed_CR\t" << library_metrics.malformed_cr << '\n'
            << "malformed_GX\t" << library_metrics.malformed_gx << '\n'
            << "malformed_GN\t" << library_metrics.malformed_gn << '\n'
            << "malformed_UB\t" << library_metrics.malformed_ub << '\n'
            << "malformed_UR\t" << library_metrics.malformed_ur << '\n'
            << "malformed_NH\t" << library_metrics.malformed_nh << '\n'
            << "malformed_AS\t" << library_metrics.malformed_as << '\n'
            << "malformed_nM\t" << library_metrics.malformed_nm << '\n'
            << "GX_not_in_feature_manifest_reads\t" << library_metrics.gx_not_feature_manifest_reads << '\n'
            << "biological_classification_manifest_supplied\t" << classes.supplied << '\n'
            << "biological_contig_classes\t" << classes.contigs.size() << '\n'
            << "biological_feature_classes\t" << classes.features.size() << '\n'
            << "biological_classified_reads\t" << library_metrics.biological_classified_reads << '\n'
            << "biological_classification_status\t"
            << classification_status(library_metrics.primary_mapped_reads,
                                     library_metrics.biological_classified_reads) << '\n'
            << "barcode_roster_union_ids\t" << barcodes.strings.size() << '\n'
            << "observed_corrected_barcodes\t" << observed_metrics.size() << '\n'
            << "ordinary_feature_ids\t" << features.size() << '\n'
            << "ordinary_feature_numeric_ids\t" << features.size() << '\n'
            << "interned_nonpacked_umis\t" << fallback_umis.size() << '\n'
            << "interned_raw_barcodes\t" << raw_barcodes.size() << '\n'
            << "barcode_rg_metric_rows\t" << barcode_rg_metrics.size() << '\n'
            << "correction_rows\t" << corrections.size() << '\n'
            << "within_rg_conflict_definition\tsame_CR_same_RG_has_multiple_corrected_CB_values\n"
            << "cross_rg_conflict_definition\ttwo_or_more_exact_RG_level_mappings_for_same_CR_disagree\n"
            << "global_conflict_definition\tsame_CR_has_multiple_corrected_CB_values_anywhere\n"
            << "within_rg_conflict_pairs\t" << within_rg_conflict_pairs << '\n'
            << "cross_rg_conflict_raw_barcodes\t" << cross_rg_conflict_raw_barcodes << '\n'
            << "global_conflict_raw_barcodes\t" << global_conflict_raw_barcodes << '\n'
            << "dense_rg_contig_slots\t" << contig_slots << '\n'
            << "observed_rg_contig_rows\t" << observed_contig_rows << '\n'
            << "molecule_hash_bin_rows\t" << hash_bin_rows << '\n'
            << "molecule_slot_bytes\t" << sizeof(MoleculeSlot) << '\n'
            << "nh_gt1_qname_state_slot_bytes\t" << sizeof(NHMultiReadSlot) << '\n'
            << "nh_gt1_qname_state_table_allocated_bytes\t"
            << nh_gt1_qname_state_table_allocated_bytes << '\n'
            << "nh_gt1_qname_bytes_used\t" << nh_gt1_qname_bytes_used << '\n'
            << "nh_gt1_qname_arena_allocated_bytes\t"
            << nh_gt1_qname_arena_allocated_bytes << '\n'
            << "ordinary_molecule_table_allocated_bytes\t" << molecules.allocated_bytes() << '\n'
            << "barcode_rg_table_allocated_bytes\t" << barcode_rg_metrics.allocated_bytes() << '\n'
            << "correction_table_allocated_bytes\t" << corrections.allocated_bytes() << '\n'
            << "hash_bin_aggregation_method\tin_place_sorted_ordinary_molecule_slots\n"
            << "hash_bin_aggregation_additional_heap_bytes\t0\n"
            << "output_sorting_method\tin_place_high_cardinality_tables_plus_numeric_barcode_order_vectors\n"
            << "record_contribution_method\tsingle_tag_decode_and_single_cigar_walk_reused_across_accumulators\n"
            << "configured_total_process_memory_bytes\t" << options.max_memory_bytes << '\n'
            << "effective_RLIMIT_AS_bytes\t" << process_limit << '\n'
            << "conservative_memory_admission_limit_bytes\t" << budget.admission_limit() << '\n'
            << "runtime_allocator_htslib_reserve_bytes\t"
            << (process_limit - budget.admission_limit()) << '\n'
            << "accounted_current_bytes\t" << budget.charged() << '\n'
            << "accounted_peak_bytes\t" << budget.peak() << '\n'
            << "process_MaxRSS_kib\t" << static_cast<uint64_t>(usage.ru_maxrss) << '\n';
        for (const auto& item : budget.component_peak()) {
            out << "memory_component_peak_bytes." << item.first << '\t'
                << item.second << '\n';
        }
        out << "htslib_bgzf_threads\t" << options.threads << '\n'
            << "logical_record_loops\t1\n"
            << "elapsed_seconds\t" << std::fixed << std::setprecision(6)
            << elapsed << '\n';
    }

    sam_hdr_destroy(header);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return profile(parse_options(argc, argv));
    } catch (const std::bad_alloc&) {
        std::cerr << "ERROR: allocation failed inside the enforced total-process memory ceiling\n";
        return 1;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
