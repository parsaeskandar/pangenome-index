# Path to sdsl-lite
SDSL_DIR ?= /Users/seeskand/Documents/sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
SOURCE_DIR=src

# Compiler - default to g++ for Linux, clang++ for macOS
MY_CXX ?= g++

# Parallelization flags
PARALLEL_FLAGS=-fopenmp -pthread

# Match root makefile: prefer vendored gbwtgraph (GBZ v2).
GBWTGRAPH_DIR ?= ../Giraffe_server/deps/gbwtgraph

ifeq ($(wildcard $(GBWTGRAPH_DIR)/lib/libgbwtgraph.a),)
GBWTGRAPH_LIBS = -lgbwtgraph
else
GBWTGRAPH_LIBS = $(GBWTGRAPH_DIR)/lib/libgbwtgraph.a
endif

# Match root makefile: Zstd for libgbwt (default ~/lib).
ZSTD_LIB_DIR ?= $(HOME)/lib
ifeq ($(strip $(ZSTD_LIB_DIR)),)
ZSTD_LDFLAGS :=
else
ZSTD_LDFLAGS := -L$(ZSTD_LIB_DIR) -Wl,-rpath,$(ZSTD_LIB_DIR)
endif
GBWT_ZSTD_LIBS = $(ZSTD_LDFLAGS) -lzstd

# Libraries - base libs without paths (paths added per-platform)
BASE_LIBS=$(GBWTGRAPH_LIBS) -lgbwt -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64 -lgrlbwt $(GBWT_ZSTD_LIBS)
GTEST_LIBS=-lgtest -lgtest_main
# Default LIBS for non-Darwin (will be overridden on Darwin)
LIBS=-L$(LIB_DIR) -L../deps/grlBWT/build $(BASE_LIBS) $(GTEST_LIBS)

# macOS-specific OpenMP & compiler handling
ifeq ($(shell uname -s), Darwin)
    MY_CXX := clang++
    HOMEBREW_PREFIX := /opt/homebrew

    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info Detected Apple Clang; configuring OpenMP with libomp...)

        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        # Find libomp installed by Homebrew or MacPorts.
        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp via Homebrew.)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
            LIBOMP_PATH := $(HOMEBREW_PREFIX)/lib
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            $(info Found keg-only libomp via Homebrew.)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LIBOMP_PATH := $(HOMEBREW_PREFIX)/opt/libomp/lib
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp via MacPorts.)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBOMP_PATH := /opt/local/lib/libomp
        else
            $(error Could not find libomp. Please install it using Homebrew or MacPorts.)
        endif
    endif

    CXX_FLAGS += -I$(HOMEBREW_PREFIX)/opt/openssl/include
    CXX_FLAGS += -I$(HOMEBREW_PREFIX)/include
    
    # Override LIBS: Homebrew path FIRST (for gtest), then other libs
    LIBS := -L$(HOMEBREW_PREFIX)/lib -L$(LIB_DIR) -L../deps/grlBWT/build $(BASE_LIBS) -L$(LIBOMP_PATH) -lomp -L$(HOMEBREW_PREFIX)/opt/openssl/lib -lcrypto $(GTEST_LIBS)
endif

CXX_FLAGS += $(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I../include -I$(GBWTGRAPH_DIR)/include -I$(INC_DIR) -I../lib -I../deps/r-index/internal -I../deps/vg -I../deps/grlBWT/include

HEADERS=$(wildcard ../include/pangenome_index/*.hpp)

#LIBOBJS=$(addprefix $(BUILD_OBJ)/,tag_arrays.o r-index.o)
#LIBRARY=$(BUILD_LIB)/libpanindexer.a

PROGRAMS=test_rindex test_bplus_tree benchmark_query_tags

LIBRARY=../lib/libpanindexer.a

.PHONY: all clean test
all:$(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)


%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

$(PROGRAMS):%:%.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)


test: $(PROGRAMS)
	$(foreach PROG,$(PROGRAMS),./$(PROG) || exit 1;)

clean:
	rm -f $(PROGRAMS) *.o
