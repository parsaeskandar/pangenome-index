SDSL_DIR=/Users/seeskand/Documents/sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
SOURCE_DIR=src
MY_CXX := clang++
OMP_HOMEBREW_PREFIX := /opt/homebrew/Cellar/libomp/18.1.8

PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lgbwt -lgbwtgraph -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64 -lgrlbwt -lgtest -lgtest_main


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

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -Iinclude -I$(INC_DIR) -I../lib -I../deps/r-index/internal -I../deps/vg -I../deps/grlBWT/include


HEADERS=$(wildcard ../include/pangenome_index/*.hpp)

#LIBOBJS=$(addprefix $(BUILD_OBJ)/,tag_arrays.o r-index.o)
#LIBRARY=$(BUILD_LIB)/libpanindexer.a


PROGRAMS=test_rindex

LIBRARY=../lib/libpanindexer.a

.PHONY: all clean test
all:$(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)


%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

$(PROGRAMS):%:%.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)


test:$(PROGRAMS)
	$(foreach PROG,$(PROGRAMS),./$(PROG) || exit 1;)

clean:
	rm -f $(PROGRAMS) *.o
