SHELL = bash
COMP = g++
CCOMP = gcc
PREFIX ?= /usr/local
HTSLIB_PREFIX ?= /nvme/software/packages/htslib/1.20
HTSLIB_INCLUDE ?= $(HTSLIB_PREFIX)/include
HTSLIB_LIB ?= $(HTSLIB_PREFIX)/lib
HTSLIB_RPATH_FLAG = -Wl,--disable-new-dtags,-rpath,$(HTSLIB_LIB)
CXXIFLAGS = -I$(HTSLIB_INCLUDE) -I$(PREFIX)/include -Iinclude
CIFLAGS = -I$(HTSLIB_INCLUDE) -I$(PREFIX)/include -Iinclude
LFLAGS = -L$(HTSLIB_LIB) $(HTSLIB_RPATH_FLAG) -L$(PREFIX)/lib -Llib
CXXFLAGS = -std=c++11 -fPIC -DBC_LENX2=$(BC_LENX2) -DKX2=$(KX2)
RNA_CXXFLAGS = -std=c++17 -fPIC -O3 -Wall -Wextra -Wno-unused-function


BC_LENX2 = 32
KX2 = 16
DEPS = lib/libhtswrapper.a
DEPS2 = -lz -lhts -lpthread
HTSWRAPPER_PATCH = patches/htswrapper_local_fixes.patch
HTSWRAPPER_PATCH_INPUT := $(wildcard $(HTSWRAPPER_PATCH))
HTSWRAPPER_BUILD_DIR = build/htswrapper_cluster
HTSWRAPPER_PATCH_STAMP = $(HTSWRAPPER_BUILD_DIR)/.local_patch_applied
HTSWRAPPER_SOURCE_FILES := $(shell find dependencies/htswrapper -type f -not -path '*/.git*' -not -path '*/build/*' -not -path '*/lib/*' 2>/dev/null)

BJP_PREFIX ?= /nvme/software/packages/align_pipelines/bjp

BINARIES = \
	atac_fq_preprocess \
	split_read_files \
	vcf_depth_filter \
	rna_bam_evidence

MAPPING_SCRIPTS := $(sort $(wildcard scripts/mapping/*.py))
WORKFLOWS := $(sort $(wildcard workflows/*.nf))
ROOT_PIPELINES = align_pipelines.nf make_ref.nf
ROOT_CONFIGS = nextflow.config nextflow.config.slurm
ROOT_HELPERS = segment_genome.py compile_star_stats.R plot_filt_stats.R

.PHONY: all install-bjp install-bjp-scripts test-rna-bam-evidence clean clean_deps

all: $(BINARIES)

atac_fq_preprocess: src/atac_fq_preprocess.cpp $(DEPS)
	$(COMP) $(CXXFLAGS) $(CXXIFLAGS) src/atac_fq_preprocess.cpp $(LFLAGS) $(DEPS) -o atac_fq_preprocess $(DEPS2)

split_read_files: src/split_read_files.cpp $(DEPS)
	$(COMP) $(CXXFLAGS) $(CXXIFLAGS) src/split_read_files.cpp $(LFLAGS) $(DEPS) -o split_read_files $(DEPS2)

vcf_depth_filter: src/vcf_depth_filter.cpp
	$(COMP) $(CXXFLAGS) $(CXXIFLAGS) src/vcf_depth_filter.cpp $(LFLAGS) -o vcf_depth_filter $(DEPS2)

rna_bam_evidence: src/rna_bam_evidence.cpp
	$(COMP) $(RNA_CXXFLAGS) $(CXXIFLAGS) src/rna_bam_evidence.cpp $(LFLAGS) -o rna_bam_evidence $(DEPS2)

# Install the complete runtime into the dedicated align_pipelines/bjp package.
# DESTDIR supports packaging into a staging root without changing BJP_PREFIX.
install-bjp: all install-bjp-scripts
	install -d $(DESTDIR)$(BJP_PREFIX)/bin
	install -m 0755 $(BINARIES) $(DESTDIR)$(BJP_PREFIX)/bin/

# Update scripts and workflows without rebuilding the compiled utilities.
install-bjp-scripts: $(MAPPING_SCRIPTS) $(WORKFLOWS) $(ROOT_PIPELINES) $(ROOT_CONFIGS) $(ROOT_HELPERS)
	install -d \
		$(DESTDIR)$(BJP_PREFIX) \
		$(DESTDIR)$(BJP_PREFIX)/bin/mapping \
		$(DESTDIR)$(BJP_PREFIX)/workflows
	install -m 0755 $(MAPPING_SCRIPTS) $(DESTDIR)$(BJP_PREFIX)/bin/mapping/
	install -m 0644 $(WORKFLOWS) $(DESTDIR)$(BJP_PREFIX)/workflows/
	install -m 0644 $(ROOT_PIPELINES) $(ROOT_CONFIGS) $(DESTDIR)$(BJP_PREFIX)/
	install -m 0755 $(ROOT_HELPERS) $(DESTDIR)$(BJP_PREFIX)/

test-rna-bam-evidence: rna_bam_evidence
	RNA_BAM_EVIDENCE_BIN=$(CURDIR)/rna_bam_evidence python3 -m unittest discover -s tests/rna_bam_evidence -v

$(HTSWRAPPER_PATCH_STAMP): $(HTSWRAPPER_SOURCE_FILES) $(HTSWRAPPER_PATCH_INPUT)
	rm -rf -- "$(HTSWRAPPER_BUILD_DIR)"
	mkdir -p "$(HTSWRAPPER_BUILD_DIR)"
	cp -a dependencies/htswrapper/. "$(HTSWRAPPER_BUILD_DIR)/"
	rm -rf -- "$(HTSWRAPPER_BUILD_DIR)/.git" "$(HTSWRAPPER_BUILD_DIR)/build" "$(HTSWRAPPER_BUILD_DIR)/lib"
	mkdir -p "$(HTSWRAPPER_BUILD_DIR)/build" "$(HTSWRAPPER_BUILD_DIR)/lib"
	if [[ -s "$(CURDIR)/$(HTSWRAPPER_PATCH)" ]]; then patch --forward --batch -d "$(HTSWRAPPER_BUILD_DIR)" -p1 < "$(CURDIR)/$(HTSWRAPPER_PATCH)"; fi
	touch "$@"

lib/libhtswrapper.a: $(HTSWRAPPER_PATCH_STAMP)
	$(MAKE) -C "$(HTSWRAPPER_BUILD_DIR)" PREFIX="$(CURDIR)" IFLAGS="-I$(HTSLIB_INCLUDE) -I$(CURDIR)/include" LFLAGS="-L$(HTSLIB_LIB) $(HTSLIB_RPATH_FLAG) -L$(CURDIR)/lib"
	$(MAKE) -C "$(HTSWRAPPER_BUILD_DIR)" install PREFIX="$(CURDIR)"

clean: clean_deps
	rm -f $(BINARIES)
	rm -f lib/libhtswrapper.a

clean_deps:
	rm -rf -- "$(HTSWRAPPER_BUILD_DIR)"
