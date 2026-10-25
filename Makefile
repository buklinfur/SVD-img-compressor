# Plain-make build for systems without CMake. CMake is the primary build;
# this file only needs a C11 compiler and GNU make (3.81+).
#
#   make            build build-make/svdc
#   make test       build and run the unit tests and a CLI round trip
#   make OPENMP=1   parallelize channels with OpenMP (gcc/clang with libgomp/libomp)
#   make clean

CC      ?= cc
CFLAGS  ?= -O3
WARN     = -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wstrict-prototypes \
           -Wmissing-prototypes -Wdouble-promotion -Wundef
STD      = -std=c11
INC      = -Isrc -isystem third_party/stb
LDLIBS   = -lm
BUILD    = build-make

ifeq ($(OPENMP),1)
CFLAGS  += -fopenmp
LDFLAGS += -fopenmp
endif

CORE     = src/matrix.c src/svd.c src/image.c src/codec.c src/metrics.c
CORE_OBJ = $(CORE:src/%.c=$(BUILD)/obj/%.o) $(BUILD)/obj/stb_impl.o

all: $(BUILD)/svdc

$(BUILD)/obj/%.o: src/%.c src/*.h
	@mkdir -p $(BUILD)/obj
	$(CC) $(STD) $(CFLAGS) $(WARN) $(INC) -c $< -o $@

# third-party code: compiled without the project's warning flags
$(BUILD)/obj/stb_impl.o: src/stb_impl.c
	@mkdir -p $(BUILD)/obj
	$(CC) $(STD) $(CFLAGS) -w $(INC) -c $< -o $@

$(BUILD)/svdc: $(BUILD)/obj/main.o $(CORE_OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BUILD)/%: tests/%.c $(CORE_OBJ) tests/test_util.h
	$(CC) $(STD) $(CFLAGS) $(WARN) $(INC) -DEXAMPLES_DIR='"examples"' $< $(CORE_OBJ) $(LDFLAGS) $(LDLIBS) -o $@

test: $(BUILD)/svdc $(BUILD)/test_svd $(BUILD)/test_codec
	$(BUILD)/test_svd
	$(BUILD)/test_codec
	$(BUILD)/svdc compress examples/test2.png $(BUILD)/roundtrip.svdc -k 20
	$(BUILD)/svdc decompress $(BUILD)/roundtrip.svdc $(BUILD)/roundtrip.png

clean:
	rm -rf $(BUILD)/obj $(BUILD)/svdc $(BUILD)/test_svd $(BUILD)/test_codec \
	       $(BUILD)/roundtrip.svdc $(BUILD)/roundtrip.png

.PHONY: all test clean
