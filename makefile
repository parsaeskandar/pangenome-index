# Path to sdsl-lite
SDSL_DIR ?= /Users/seeskand/Documents/sdsl-lite
include $(SDSL_DIR)/Make.helper

# Directories
BUILD_BIN = bin
BUILD_LIB = lib
BUILD_OBJ = obj
SOURCE_DIR = src

# Compiler
MY_CXX ?= g++

# Initial flags
CXX_FLAGS += $(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS)
# gbwtgraph headers must precede $(INC_DIR): older system installs lack GBZ v2 support.
GBWTGRAPH_DIR ?= Giraffe_server/deps/gbwtgraph
CXX_FLAGS += -Iinclude -I$(GBWTGRAPH_DIR)/include -I$(INC_DIR) -Ideps/vg -Ideps/grlBWT/include -UNDEBUG

# Link vendored libgbwtgraph.a when built (make gbwtgraph-lib); else system -lgbwtgraph from $(LIB_DIR).
ifeq ($(wildcard $(GBWTGRAPH_DIR)/lib/libgbwtgraph.a),)
GBWTGRAPH_LIBS = -lgbwtgraph
$(info gbwtgraph: linking from $(LIB_DIR). If GBZ load fails on v2 graphs, run: make gbwtgraph-lib)
else
GBWTGRAPH_LIBS = $(GBWTGRAPH_DIR)/lib/libgbwtgraph.a
$(info gbwtgraph: linking vendored $(GBWTGRAPH_DIR)/lib/libgbwtgraph.a)
endif

# Zstandard (required by libgbwt). Default ~/lib for user installs; set ZSTD_LIB_DIR= for system-only paths.
ZSTD_LIB_DIR ?= $(HOME)/lib
ifeq ($(strip $(ZSTD_LIB_DIR)),)
ZSTD_LDFLAGS :=
else
ZSTD_LDFLAGS := -L$(ZSTD_LIB_DIR) -Wl,-rpath,$(ZSTD_LIB_DIR)
endif
GBWT_ZSTD_LIBS = $(ZSTD_LDFLAGS) -lzstd

# Parallelization flags
PARALLEL_FLAGS = -fopenmp -pthread

# Libraries (gbwtgraph before gbwt for correct static resolution)
LIBS = -L$(LIB_DIR) -Ldeps/grlBWT/build $(GBWTGRAPH_LIBS) -lgbwt -lhandlegraph -lsdsl -lgrlbwt -lcrypto $(GBWT_ZSTD_LIBS)

# macOS-specific OpenMP & compiler handling
ifeq ($(shell uname -s), Darwin)
    MY_CXX := clang++
    HOMEBREW_PREFIX := /opt/homebrew

    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info Detected Apple Clang; configuring OpenMP with libomp...)
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp via Homebrew.)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            $(info Found keg-only libomp via Homebrew.)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LIBS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp via MacPorts.)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        else
            $(error Could not find libomp. Please install it using Homebrew or MacPorts.)
        endif

        LIBS += -lomp
    endif

    # OpenSSL (keg-only on Homebrew)
    LIBS += -L$(HOMEBREW_PREFIX)/opt/openssl/lib
    CXX_FLAGS += -I$(HOMEBREW_PREFIX)/opt/openssl/include

    CXX_FLAGS += -I$(HOMEBREW_PREFIX)/include
endif

# Headers and objects
HEADERS = $(wildcard include/pangenome_index/*.hpp)
LIBOBJS = $(addprefix $(BUILD_OBJ)/,r-index.o tag_arrays.o sampled_tag_array.o translation_tables.o)
LIBRARY = $(BUILD_LIB)/libpanindexer.a

PROGRAMS = $(addprefix $(BUILD_BIN)/,build_tags merge_tags build_rindex query_tags tags_check find_mems convert_tags print_stats build_sampled_tags query_sampled_tags coordinate_translation build_translation_tables)

# Targets
.PHONY: all clean directories grlbwt gbwtgraph-lib test

all: grlbwt directories $(LIBRARY) $(PROGRAMS)

grlbwt:
	mkdir -p deps/grlBWT/build
	cd deps/grlBWT/build && cmake .. && make

# Build static libgbwtgraph.a into $(GBWTGRAPH_DIR)/lib (requires SDSL_DIR, gbwt/handlegraph in LIB_DIR).
gbwtgraph-lib:
	$(MAKE) -C "$(GBWTGRAPH_DIR)" SDSL_DIR="$(SDSL_DIR)" all

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)

$(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ):
	mkdir -p $@

$(BUILD_OBJ)/%.o: $(SOURCE_DIR)/%.cpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXX_FLAGS) -c -o $@ $<

$(LIBRARY): $(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

$(BUILD_BIN)/%: $(BUILD_OBJ)/%.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test: $(LIBRARY)
	cd tests && $(MAKE) test

# ── Python extension (.so) ──────────────────────────────────────────────────
PYBIND11_INCLUDES = $(shell python3 -m pybind11 --includes)
PYTHON_EXT_SUFFIX = $(shell python3-config --extension-suffix)

# Object files compiled with -fPIC for the shared library
$(BUILD_OBJ)/pangenome_server.pic.o: $(SOURCE_DIR)/pangenome_server.cpp $(SOURCE_DIR)/pangenome_server.hpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXX_FLAGS) -fPIC -c -o $@ $<

$(BUILD_OBJ)/bindings.pic.o: $(SOURCE_DIR)/bindings.cpp $(SOURCE_DIR)/pangenome_server.hpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXX_FLAGS) -fPIC $(PYBIND11_INCLUDES) -c -o $@ $<

$(BUILD_OBJ)/coordinate_translation.pic.o: $(SOURCE_DIR)/coordinate_translation.cpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXX_FLAGS) -fPIC -c -o $@ $<

# Build fPIC versions of the library objects needed by the .so
LIBOBJS_PIC = $(addprefix $(BUILD_OBJ)/,r-index.pic.o tag_arrays.pic.o sampled_tag_array.pic.o translation_tables.pic.o)

$(BUILD_OBJ)/%.pic.o: $(SOURCE_DIR)/%.cpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXX_FLAGS) -fPIC -c -o $@ $<

PYBIND11_LDFLAGS =
ifeq ($(shell uname -s), Darwin)
    PYBIND11_LDFLAGS += -undefined dynamic_lookup
endif

liftover_ext$(PYTHON_EXT_SUFFIX): directories $(BUILD_OBJ)/pangenome_server.pic.o $(BUILD_OBJ)/bindings.pic.o $(BUILD_OBJ)/coordinate_translation.pic.o $(LIBOBJS_PIC)
	$(MY_CXX) $(LDFLAGS) $(CXX_FLAGS) -shared -fPIC $(PYBIND11_LDFLAGS) -o $@ \
		$(BUILD_OBJ)/pangenome_server.pic.o \
		$(BUILD_OBJ)/bindings.pic.o \
		$(BUILD_OBJ)/coordinate_translation.pic.o \
		$(LIBOBJS_PIC) \
		$(LIBS)

liftover_ext.so: liftover_ext$(PYTHON_EXT_SUFFIX)
	@if [ "liftover_ext$(PYTHON_EXT_SUFFIX)" != "liftover_ext.so" ]; then \
		ln -sf liftover_ext$(PYTHON_EXT_SUFFIX) liftover_ext.so; \
	fi

clean:
	rm -rf $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)
	rm -f *.o *.a $(OBSOLETE) liftover_ext$(PYTHON_EXT_SUFFIX) liftover_ext.so
