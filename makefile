SDSL_DIR=/Users/seeskand/Documents/sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
SOURCE_DIR=src
MY_CXX := clang++
OMP_HOMEBREW_PREFIX := /opt/homebrew/Cellar/libomp/18.1.8

PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lgbwt -lgbwtgraph -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64


# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        # If HOMEBREW_PREFIX is specified, libomp probably cannot be found automatically.
        ifdef HOMEBREW_PREFIX
            ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
                # libomp moved to these directories, recently, because it is now keg-only to not fight GCC
                PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
                LIBS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
            else
                PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
                LIBS += -L$(HOMEBREW_PREFIX)/lib
            endif
        # Macports installs libomp to /opt/local/lib/libomp
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        endif

        # We also need to link it.
        LIBS += -lomp
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -Iinclude -I$(INC_DIR) -Ideps/r-index/internal -Ideps/vg


HEADERS=$(wildcard include/pangenome_index/*.hpp)

LIBOBJS=$(addprefix $(BUILD_OBJ)/,bplus_tree.o tag_arrays.o unique_kmer.o algorithm.o)
LIBRARY=$(BUILD_LIB)/libpanindexer.a


PROGRAMS=$(addprefix $(BUILD_BIN)/,build_tags query_tags)


.PHONY: all clean directories test
all: directories $(LIBRARY) $(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)

$(BUILD_BIN):
	mkdir -p $@

$(BUILD_LIB):
	mkdir -p $@

$(BUILD_OBJ):
	mkdir -p $@

$(BUILD_OBJ)/%.o:$(SOURCE_DIR)/%.cpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

$(LIBRARY):$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

$(BUILD_BIN)/%:$(BUILD_OBJ)/%.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test:$(LIBRARY)
	cd tests && $(MAKE) test

clean:
	rm -rf $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)
	rm -f *.o *.a $(OBSOLETE)
	#cd tests && $(MAKE) clean
