# Pangenome-Index

A comprehensive toolkit for building and querying pangenome indices using tag arrays and r-index structures. This project provides efficient algorithms for indexing pangenome graphs and performing various queries including MEM (Maximal Exact Match) finding and tag-based queries.

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Large Graph Processing](#large-graph-processing)
- [Executables](#executables)
  - [build_tags](#build_tags)
  - [build_rindex](#build_rindex)
  - [merge_tags](#merge_tags)
  - [find_mems](#find_mems)
  - [query_tags](#query_tags)
- [File Formats](#file-formats)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

## Dependencies

### Required Dependencies

- **[SDSL](https://github.com/vgteam/sdsl-lite)** (vgteam fork) - For low-level data structures
- **[GBWT](https://github.com/jltsiren/gbwt)** - For the backend
- **[GBWTGraph](https://github.com/jltsiren/gbwtgraph)** - For the graph-based backend
- **[HandleGraph](https://github.com/vgteam/libhandlegraph)** - For the handle graph interface
- **[grlBWT](https://github.com/ddiazdom/grlBWT)** - For the run-length BWT interface (automatically built during compilation)

Please refer to each dependency's repository for installation instructions.

### System Dependencies

- **C++ Compiler**: GCC 7+ or Clang 6+
- **OpenMP**: For parallel processing
- **CMake**: For building grlBWT
- **Make**: For building the project

**macOS Users**: Install OpenMP via Homebrew:
```bash
brew install libomp
```

## Installation

### Step 1: Install Dependencies

First, install all required dependencies in accessible locations. The makefile looks for dependencies in `$(LIB_DIR)` which should be set to the location where you installed the dependencies. Follow the installation instructions provided in the [Dependencies](#dependencies) section above.

### Step 2: Clone and Build Pangenome-Index

```bash
# Clone the repository with submodules
git clone --recursive https://github.com/parsaeskandar/pangenome-index.git
cd pangenome-index

# Edit the Makefile to set the correct path to your SDSL installation
# Change this line in the Makefile:
# SDSL_DIR ?= /Users/seeskand/Documents/sdsl-lite
# To point to your SDSL installation location

# Build the project
make -j 8
```

The build process will:
1. Automatically build grlBWT dependency
2. Create the `bin/` directory with all executables
3. Create the `lib/` directory with the library
4. Create the `obj/` directory with object files

## Quick Start

For small to medium-sized graphs, you can run the entire pipeline using only the `build_tags` binary:

### Step 1: Prepare Your Graph

You need a graph file in `.gbz` format. If you don't have one, you can create it from your graph using GBWTGraph tools.

### Step 2: Extract Sequences and Create RL-BWT

```bash
# Extract sequences from the graph (requires gbwtgraph installation)
gbz_extract -t 8 -b your_graph.gbz > graph_info

# Create the run-length BWT file using grlbwt
grlbwt-cli -t 8 graph_info
# This creates graph_info.rl_bwt
```

### Step 3: Build Tag Arrays

```bash
# Build tag arrays using the graph and RL-BWT files
./bin/build_tags your_graph.gbz graph_info.rl_bwt output.tags
```

## Large Graph Processing

For very large graphs (e.g., whole human genome), the pipeline should be run per-chromosome and then merged. This approach reduces memory usage and enables parallel processing.

### Step 1: Process Each Chromosome

Create a script to process each chromosome independently:

```bash
#!/bin/bash
# Process a single chromosome
CHR=$1
THREADS=16

# Create chromosome-specific directory
mkdir -p "${CHR}"
cd "${CHR}"

# Extract sequences from chromosome graph
gbz_extract -t ${THREADS} -b -p ${CHR}.gbz > ${CHR}_info

# Create run-length BWT for chromosome
grlbwt-cli -t ${THREADS} -T ${PWD}/../tmp -f 0 ${CHR}_info

# Build tag arrays for chromosome
./bin/build_tags ${CHR}.gbz ${CHR}_info.rl_bwt ${CHR}.tags

cd ..
```

Run this for each chromosome (e.g., chr1, chr2, ..., chrX, chrY).

### Step 2: Create Whole-Genome R-Index

```bash
# Extract sequences from whole-genome graph
gbz_extract -t 32 -b -p whole_genome.gbz > whole_genome_info

# Create whole-genome run-length BWT
grlbwt-cli -t 32 -T ${PWD}/tmp whole_genome_info

# Build whole-genome r-index
./bin/build_rindex whole_genome_info.rl_bwt > whole.ri
```

### Step 3: Merge Chromosome Tag Arrays

```bash
# Merge all chromosome tag arrays into whole-genome index
./bin/merge_tags whole_genome.gbz whole.ri tags/
```

**Directory Structure After Processing:**
```
project/
├── whole_genome.gbz
├── whole.ri
├── tags/
│   ├── chr1.tags
│   ├── chr2.tags
│   ├── ...
│   └── chrY.tags
└── whole_genome_tag_array_compressed.tags
```

## Executables

### build_tags

**Purpose**: Constructs tag arrays over the pangenome graph.

**Input Files Required**:
- `.gbz` file: The pangenome graph in GBZ format
- `.rl_bwt` file: Run-length BWT file created by grlbwt

**How to Create Input Files**:
```bash
# Create .gbz file (if you don't have one)
# This depends on your graph format - convert to GBZ using GBWTGraph tools

# Create .rl_bwt file from .gbz
gbz_extract -t 8 -b your_graph.gbz > graph_info
grlbwt-cli -t 8 graph_info
```

**Usage**:
```bash
./bin/build_tags <graph.gbz> <graph_info.rl_bwt> <output.tags>
```

**Example**:
```bash
./bin/build_tags test_data/x.giraffe.gbz test_data/x.rl_bwt output.tags
```

**Output**: 
- `output.tags`: Binary file containing the tag arrays index

**What it does**:
1. Loads the pangenome graph from the GBZ file
2. Computes unique k-mers in the graph (default k=31)
3. Builds a B+ tree structure for efficient k-mer lookup
4. Creates tag arrays that map k-mers to their positions in the graph
5. Serializes the tag arrays to the output file

---

### build_rindex

**Purpose**: Builds the r-index over the pangenome from an RL-BWT file.

**Input Files Required**:
- `.rl_bwt` file: Run-length BWT file created by grlbwt

**Usage**:
```bash
./bin/build_rindex <graph_info.rl_bwt>
```

**Example**:
```bash
./bin/build_rindex test_data/x.rl_bwt
```

**Output**: 
- R-index data printed to stdout (can be redirected to a file)

**What it does**:
1. Loads the run-length BWT file
2. Constructs the r-index data structure
3. Serializes the r-index to stdout

---

### merge_tags

**Purpose**: Merges multiple tag arrays from different chromosomes into a whole-genome tag array index.

**Input Files Required**:
- Whole-genome r-index file
- Directory containing chromosome-specific `.rl_bwt` files
- Multiple `.tags` files (one per chromosome)

**Usage**:
```bash
./bin/merge_tags <whole_genome_rindex> <rl_bwt_directory> <output_merged.tags>
```

**Example**:
```bash
./bin/merge_tags whole_genome.ri /path/to/chromosome_rlbwt_files/ merged_output.tags
```

**Output**: 
- `merged_output.tags`: Merged tag arrays file

**What it does**:
1. Loads the whole-genome r-index
2. Reads multiple chromosome-specific RL-BWT files
3. Merges tag arrays from different chromosomes
4. Creates a unified tag array index for the whole genome

---

### find_mems

**Purpose**: Finds Maximal Exact Matches (MEMs) between query sequences and the pangenome, reporting unique tags of those MEMs.

**Input Files Required**:
- R-index file (`.ri`)
- Tag arrays file (`.tags`)
- Reads file (text file with one sequence per line)

**Usage**:
```bash
./bin/find_mems <r_index.ri> <tag_arrays.tags> <reads.txt> <min_mem_length> <min_occurrences>
```

**Example**:
```bash
./bin/find_mems test_data/x.giraffe.ri test_data/x.tags.index test_data/small_test_nl.txt 10 1
```

**Output**: 
- MEMs found in the pangenome with their positions and tag information

**What it does**:
1. Loads the r-index and tag arrays
2. Reads query sequences from the input file
3. Finds all MEMs between queries and the pangenome
4. Reports MEM positions and associated tag information
5. Filters results by minimum MEM length and occurrence count

---

### query_tags

**Purpose**: Queries the pangenome using tag arrays for k-mer lookups. This is for evaluating the performance of finding the k-mers on the tag arrays.

**Input Files Required**:
- R-index file (`.ri`)
- Tag arrays file (`.tags`)
- Sequence file (text file with one sequence per line)

**Usage**:
```bash
./bin/query_tags <r_index.ri> <tag_arrays.tags> <k-mer_length> <sequences.txt>
```

**Example**:
```bash
./bin/query_tags test_data/x.giraffe.ri test_data/x.tags.index 31 test_data/small_test_nl.txt
```

**Output**: 
- K-mer lookup results showing positions in the pangenome

**What it does**:
1. Loads the r-index and tag arrays
2. Reads query sequences
3. Extracts k-mers from the sequences (Randomly selected)
4. Performs lookups in the tag arrays
5. Reports k-mer positions in the pangenome

---


## File Formats

### Input Files

1. **`.gbz` files**: Pangenome graphs in GBZ format (GBWTGraph format)
2. **`.rl_bwt` files**: Run-length BWT files created by grlbwt
3. **`.ri` files**: R-index files (serialized r-index data)
4. **`.tags` files**: Tag arrays index files (binary format)
5. **Text files**: Plain text files with one sequence per line

### Output Files

1. **`.tags` files**: Binary tag arrays index files
2. **R-index data**: Serialized r-index output (can be saved to `.ri` files)
3. **Console output**: MEM results, query results, validation reports

## Examples

### Complete Pipeline Example

```bash
# 1. Prepare your graph (assuming you have a graph in GBZ format)
GRAPH_FILE="your_graph.gbz"

# 2. Extract sequences and create RL-BWT
gbz_extract -t 8 -b $GRAPH_FILE > graph_info
grlbwt-cli -t 8 graph_info

# 3. Build tag arrays
./bin/build_tags $GRAPH_FILE graph_info.rl_bwt output.tags

# 4. Build r-index
./bin/build_rindex graph_info.rl_bwt > output.ri

# 5. Query with reads
./bin/find_mems output.ri output.tags your_reads.txt 10 1
```

### Working with Test Data

```bash
# Use the provided test data
./bin/build_tags test_data/x.giraffe.gbz test_data/x.rl_bwt test_output.tags
./bin/build_rindex test_data/x.rl_bwt > test_output.ri
./bin/find_mems test_output.ri test_output.tags test_data/small_test_nl.txt 5 1
```

## Troubleshooting

### Common Issues

1. **Dependency not found errors**:
   - Ensure all dependencies are installed and built
   - Check that `SDSL_DIR` in the Makefile points to the correct location
   - Verify that library files are in the expected locations

2. **OpenMP errors on macOS**:
   - Install libomp: `brew install libomp`
   - The Makefile should automatically detect and configure OpenMP

3. **Memory issues with large graphs**:
   - For very large graphs, consider processing per-chromosome
   - Use the `merge_tags` tool to combine chromosome-specific results

4. **File format errors**:
   - Ensure `.gbz` files are valid GBWTGraph format
   - Verify `.rl_bwt` files are created correctly with grlbwt
   - Check that input text files have proper line endings

### Getting Help

If you encounter issues:

1. Check that all dependencies are properly installed
2. Verify file formats and paths
3. Try with the provided test data first
4. Check the console output for specific error messages

## Citation

If you use this software in your research, please cite:

Parsa Eskandar, Benedict Paten, and Jouni Sirén: Lossless Pangenome Indexing Using Tag Arrays. bioRxiv 2025.05.12.653561, 2025. DOI: 10.1101/2025.05.12.653561

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.





