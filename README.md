# pangenome-index



# Dependencies

- [SDSL](https://github.com/vgteam/sdsl-lite) (vgteam fork) for low-level data structures. 
- [GBWT](https://github.com/jltsiren/gbwt) for the backend.
- [GBWTGraph](https://github.com/jltsiren/gbwtgraph) for the graph-based backend.
- [HandleGraph](https://github.com/vgteam/libhandlegraph) for the handle graph interface.
- [grlBWT](https://github.com/ddiazdom/grlBWT) for the run-length BWT interface.

# Compiling pangenome-index

Clone the repository:
```sh
git clone --recursive https://github.com/parsaeskandar/pangenome-index.git
cd pangenome-index
```
Edit the `Makefile` to set the correct path to your `sdsl-lite` installation:
```sh
SDSL_DIR=/path/to/sdsl-lite
```
Build the Project:
```sh
make
```

# Usage

The toolkit provides several command-line programs:

build_tags: Constructs tag arrays over the pangenome graph.

build_rindex: Builds the r-index over the pangenome.

merge_tags: Merges multiple tag arrays.

find_mems: Finds Maximal Exact Matches between query sequences and the pangenome. Also report unique tags of those MEMs.

query_tags: Queries the pangenome using tag arrays. 





