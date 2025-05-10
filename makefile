# Path to sdsl-lite
SDSL_DIR ?= /Users/seeskand/Documents/sdsl-lite
include $(SDSL_DIR)/Make.helper

# Directories
BUILD_BIN = bin
BUILD_LIB = lib
BUILD_OBJ = obj
SOURCE_DIR = src

# Default compiler
MY_CXX ?= g++

# Default parallelization flags
PARALLEL_FLAGS = -fopenmp -pthread

# Default library list
LIBS = -L$(LIB_DIR) -lgbwtgraph -lgbwt -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64 -lgrlbwt

# macOS-specific settings
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
endif

# Compiler flags
CXX_FLAGS = $(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -Iinclude -I$(INC_DIR) -Ideps/r-index/internal -Ideps/vg -Ideps/grlBWT/include -UNDEBUG

# macOS-specific includes
ifeq ($(shell uname -s), Darwin)
    CXX_FLAGS += -I$(HOMEBREW_PREFIX)/include
endif

# Header and object definitions
HEADERS = $(wildcard include/pangenome_index/*.hpp)
LIBOBJS = $(addprefix $(BUILD_OBJ)/,r-index.o tag_arrays.o)
LIBRARY = $(BUILD_LIB)/libpanindexer.a

PROGRAMS = $(addprefix $(BUILD_BIN)/,build_tags merge_tags build_rindex query_tags tags_check find_mems)

.PHONY: all clean directories test

all: directories $(LIBRARY) $(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)

$(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ):
	mkdir -p $@

$(BUILD_OBJ)/%.o: $(SOURCE_DIR)/%.cpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

$(LIBRARY): $(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

$(BUILD_BIN)/%: $(BUILD_OBJ)/%.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test: $(LIBRARY)
	cd tests && $(MAKE) test

clean:
	rm -rf $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)
	rm -f *.o *.a $(OBSOLETE)
	# cd tests && $(MAKE) clean
