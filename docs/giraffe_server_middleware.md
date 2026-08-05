# Giraffe Server Middleware

This repo now includes middleware for a long-lived `vg giraffe-server` process:

- `middleware/giraffe_server_middleware.py`: starts/stops `vg giraffe-server`, sends query batches, receives framed results.
- `middleware/pangenome_middleware.py`: unified API for coordinate translation + read mapping.

## Minimal mapping example

```python
from middleware.giraffe_server_middleware import GiraffeServerConfig, GiraffeServerMiddleware

cfg = GiraffeServerConfig(
    vg_binary="/private/groups/cgl/seeskand/1.server/Giraffe_server/bin/vg",
    gbz_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/graph.gbz",
    minimizer_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/index.shortread.withzip.min",
    distance_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/index.dist",
    zipcode_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/index.shortread.zipcodes",
    threads=8,
    max_multimaps=1,
    batch_size=256,
)

mw = GiraffeServerMiddleware(cfg)
result = mw.map_reads([("r1", "ACTAGAGAGA", "IIIIIIIIII")])
print(result)
mw.stop()
```

The middleware sends an internal `FLUSH_NOW` control line after each request batch so
`vg giraffe-server` processes immediately even when `--batch-size` is larger than the request size.

## Unified middleware example

```python
from middleware.giraffe_server_middleware import GiraffeServerConfig
from middleware.pangenome_middleware import CoordinateIndexPaths, PangenomeMiddleware

coord = CoordinateIndexPaths(
    gbz_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/graph.gbz",
    ri_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/rlbwt_rindex.ri",
    tags_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/sampled.tags",
    gbwt_ri_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/gbwt_fastlocate.ri",
    table1_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/output.t1",
    table2_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/output.t2",
)

gcfg = GiraffeServerConfig(
    vg_binary="/private/groups/cgl/seeskand/1.server/Giraffe_server/bin/vg",
    gbz_path=coord.gbz_path,
    minimizer_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/index.shortread.withzip.min",
    distance_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/index.dist",
    zipcode_path="/private/groups/cgl/seeskand/graph_index/11.server/test_files/chrM/index.shortread.zipcodes",
)

svc = PangenomeMiddleware(coord, gcfg)
print(svc.map_sequences(["ACTAGAGAGA", "GGCCAGTGCCCTCCTAGTTGGGGGGTAGGGGC"]))
print(svc.translate("CHM13#0", 100, 120, "HG002#1"))
svc.close()
```

