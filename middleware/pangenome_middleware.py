#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Tuple

import liftover_ext

from .giraffe_server_middleware import (
    FastqRead,
    GiraffeAlignment,
    GiraffeServerConfig,
    GiraffeServerMiddleware,
)


@dataclass
class CoordinateIndexPaths:
    gbz_path: str
    ri_path: str
    tags_path: str
    gbwt_ri_path: str
    table1_path: str
    table2_path: str


class PangenomeMiddleware:
    """
    Unified middleware:
    - coordinate translation via liftover_ext.Index (in-process)
    - read mapping via vg giraffe-server (long-lived subprocess)
    """

    def __init__(self, coord_paths: CoordinateIndexPaths, giraffe_cfg: GiraffeServerConfig) -> None:
        self._coord = liftover_ext.Index()
        self._coord.load(
            coord_paths.gbz_path,
            coord_paths.ri_path,
            coord_paths.tags_path,
            coord_paths.gbwt_ri_path,
            coord_paths.table1_path,
            coord_paths.table2_path,
        )
        self._giraffe = GiraffeServerMiddleware(giraffe_cfg)
        self._giraffe.start()

    def close(self) -> None:
        self._giraffe.stop()

    def translate(self, src_haplotype: str, start: int, end: int, tgt_haplotype: str):
        return self._coord.translate(src_haplotype, start, end, tgt_haplotype)

    def map_reads(self, reads: Iterable[FastqRead]) -> List[List[str]]:
        return self._giraffe.map_reads(reads)

    def map_reads_with_graph_paths(self, reads: Iterable[FastqRead]) -> List[List[GiraffeAlignment]]:
        return self._giraffe.map_reads_with_graph_paths(reads)

    def map_sequences(self, sequences: Iterable[str]) -> List[List[str]]:
        return self._giraffe.map_sequences(sequences)

    def map_sequences_with_graph_paths(self, sequences: Iterable[str]) -> List[List[GiraffeAlignment]]:
        return self._giraffe.map_sequences_with_graph_paths(sequences)

    def accept_graph_path(self, graph_path: List[Tuple[int, bool, int]]) -> None:
        visits = [liftover_ext.GraphPathVisit(int(n), bool(rev), int(fl)) for n, rev, fl in graph_path]
        self._coord.accept_graph_path(visits)
