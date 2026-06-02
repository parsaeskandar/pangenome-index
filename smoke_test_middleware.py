#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys

from middleware.giraffe_server_middleware import GiraffeServerConfig, GiraffeServerMiddleware

DEFAULT_SURJECT_TARGET = "GRCh38"

HELP_TEXT = """\
Commands:
  <sequence>                                  map (no surjection, or default --surject-target if given)
  <sequence> <target>                         map + surject to <target>     (one giraffe-server call)
  <sequence> <target> anchor                  map → build anchors → surject-with-anchors
                                              (3-step pipeline: giraffe-server, coord-trans, giraffe-server)
  translate <src> <start> <end> <tgt>         coordinate translation
  haplotypes                                  list all haplotype names (requires --ri/--tags/--gbwt-ri/--t1/--t2)
  help                                        show this message
  quit / exit / Ctrl-D                        exit
"""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Smoke test / interactive session for the pangenome middleware",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--vg", required=True, help="Path to vg binary")
    p.add_argument("--gbz", required=True, help="Path to graph.gbz")
    p.add_argument("--minimizer", required=True, help="Path to minimizer index")
    p.add_argument("--dist", required=True, help="Path to distance index")
    p.add_argument("--zipcodes", required=True, help="Path to zipcode index")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--max-multimaps", type=int, default=0,
                   help="Override max mappings per read. 0 (default) = use the "
                        "engine's built-in BLAT default of 100.")
    p.add_argument("--surject-target", default="",
                   help="Pre-index this path/haplotype name for surjection (e.g. GRCh38#0#chr1). "
                        "If omitted, sequences are mapped to the graph only unless you specify "
                        "a target inline in the REPL.")

    # Interactive mode flag
    p.add_argument("--interactive", action="store_true",
                   help="Start an interactive REPL session instead of running the smoke test")

    # Anchor pipeline flag (non-interactive smoke test only).
    p.add_argument("--anchors", action="store_true",
                   help="Route the smoke-test reads through the 3-step anchor pipeline "
                        "(graph-map → build anchors → surject-with-anchors) instead of "
                        "the regular single-call map+surject. Requires --surject-target "
                        "and the coordinate index paths (--ri/--tags/--gbwt-ri/--t1/--t2). "
                        "Useful for diffing anchor-driven vs ReferencePathOverlay-driven "
                        "surjection output side-by-side.")

    # Optional coordinate-index paths (only used in --interactive with translate support)
    coord = p.add_argument_group("coordinate index (optional, enables 'translate' command)")
    coord.add_argument("--ri",       help="RLBWT r-index (.ri)")
    coord.add_argument("--tags",     help="Sampled tag array (.tags)")
    coord.add_argument("--gbwt-ri",  help="GBWT FastLocate r-index (.ri)")
    coord.add_argument("--t1",       help="Translation table 1 (.t1)")
    coord.add_argument("--t2",       help="Translation table 2 (.t2)")

    return p.parse_args()


_SMOKE_READS = [
    ("r1", "ACTAGAGAGA", "IIIIIIIIII"),
    ("r2", "GGCCAGTGCCCTCCTAGTTGGGGGGTAGGGGC", "I" * 32),
]


def _run_smoke_test(mw: GiraffeServerMiddleware,
                    surject_target: str = "") -> int:
    """Regular smoke test: one map+surject call per read. With `surject_target`
    set, every read gets surjection results appended as GAF tags
    (sj:Z:..., sn:Z:..., sp:i:..., etc.) — that's what we compare against
    the --anchors output."""
    print(f"# mode=regular target={surject_target or '(none)'}")
    reads: list = []
    for name, seq, qual in _SMOKE_READS:
        if surject_target:
            reads.append((name, seq, qual, surject_target))
        else:
            reads.append((name, seq, qual))
    try:
        out = mw.map_reads(reads)
    finally:
        mw.stop()

    if len(out) != len(reads):
        print(f"FAIL: expected {len(reads)} result groups, got {len(out)}", file=sys.stderr)
        return 1

    non_empty = 0
    for i, read_out in enumerate(out):
        print(f"[read {i}] mappings={len(read_out)}")
        for line in read_out:
            print(line)
        if read_out:
            non_empty += 1

    if non_empty == 0:
        print("FAIL: received no mappings for any read", file=sys.stderr)
        return 1

    print("PASS: middleware returned framed mapping output.")
    return 0


def _run_anchor_smoke_test(mw: GiraffeServerMiddleware,
                           coord_index,
                           surject_target: str) -> int:
    """Anchor-mode smoke test: same reads as `_run_smoke_test`, but routed
    through the 3-step anchor pipeline (graph-map → build anchors →
    surject-with-anchors) instead of a single map+surject call. The output
    GAF tags follow the same schema as the regular path so the two runs
    can be diffed line-by-line.
    """
    if coord_index is None:
        print("FAIL: --anchors requires --ri/--tags/--gbwt-ri/--t1/--t2 "
              "(coordinate index not loaded)", file=sys.stderr)
        return 1
    if not surject_target:
        print("FAIL: --anchors requires --surject-target", file=sys.stderr)
        return 1

    print(f"# mode=anchors target={surject_target}")
    non_empty = 0
    try:
        for read_idx, (name, seq, qual) in enumerate(_SMOKE_READS):
            print(f"\n=== Anchor pipeline: {name} -> {surject_target} ===")
            # Step 1 — graph alignment only.
            out = mw.map_reads([(name, seq, qual)])
            graph_alignments = out[0] if out else []
            print(f"[read {read_idx}] step1 graph_alignments={len(graph_alignments)}")
            for gl in graph_alignments:
                print(gl)
            if not graph_alignments:
                continue

            # Step 2 — anchors per alignment.
            use_full = hasattr(coord_index, "build_surject_anchors_full")
            anchor_records = []
            for gaf in graph_alignments:
                if use_full:
                    res = coord_index.build_surject_anchors_full(gaf, surject_target)
                    anchor_records.append((list(res.anchors), res.target_path_length))
                else:
                    anchor_records.append(
                        (coord_index.build_surject_anchors(gaf, surject_target), 0)
                    )

            # Step 3 — surject each using the anchors from step 2.
            surjected_lines = 0
            for gaf, (anchors, path_len) in zip(graph_alignments, anchor_records):
                if not anchors:
                    # No anchors → we can't drive the anchor surjector. Emit
                    # the graph alignment unchanged to keep parity with the
                    # regular smoke test's output count.
                    print(gaf + "\tsj:Z:no_anchors")
                    surjected_lines += 1
                    continue
                surjected = mw.surject_with_anchors(
                    gaf, anchors, surject_target, target_path_length=path_len,
                )
                for line in surjected:
                    print(line)
                    surjected_lines += 1
            if surjected_lines:
                non_empty += 1
    finally:
        mw.stop()

    if non_empty == 0:
        print("FAIL: no anchor-pipeline surjections completed for any read",
              file=sys.stderr)
        return 1

    print("PASS: anchor-pipeline smoke test completed.")
    return 0


def _run_anchor_path(mw: GiraffeServerMiddleware, coord_index, seq: str, target: str) -> None:
    """REPL handler for the 3-step anchor pipeline:
      step 1: giraffe-server map (no surjection)         — map_reads
      step 2: coord-trans build_surject_anchors_full     — liftover_ext.Index
      step 3: giraffe-server surject-with-anchors        — SURJECT_WITH_ANCHORS

    Step 3 talks to giraffe-server via the SURJECT_WITH_ANCHORS stdin command,
    so the surjection results returned here come from the AnchorBackedPosition-
    Graph driven Surjector — not from a fallback map+surject call. To compare
    against the ReferencePathOverlay path, use the same sequence and target
    without the `anchor` keyword (or run the non-interactive smoke test with
    and without --anchors).

    The hasattr-gated TODO branches below remain only as a safety net for an
    incomplete build (e.g. liftover_ext.so missing the new bindings); they
    are never hit when the project is built normally.
    """
    name = f"q_{seq[:8]}"
    qual = "I" * len(seq)

    print(f"=== Anchor pipeline: {name} -> {target} ===")

    # Step 1 — giraffe-server: map without surjection. Already wired:
    # passing no target produces a graph alignment, no surjection step.
    print(f"\n[step 1] giraffe-server: map (no surjection)")
    try:
        out = mw.map_reads([(name, seq, qual)])
    except Exception as exc:
        print(f"  ERROR mapping: {exc}", file=sys.stderr)
        return

    graph_alignments = out[0] if out else []
    print(f"  got {len(graph_alignments)} graph alignment(s):")
    for gaf in graph_alignments:
        print(f"    {gaf}")
    if not graph_alignments:
        print("  (no graph alignments returned — aborting)")
        return

    # Step 2 — coord-trans: build anchors against the requested target. This
    # needs a Python binding for build_surject_anchors_for_path. If the
    # binding isn't there, say what's missing and stop with the partial
    # results so the user can still inspect step 1's output.
    print(f"\n[step 2] coord-trans: build anchors for target='{target}'")
    if coord_index is None:
        print("  ERROR: coordinate index not loaded; pass --ri/--tags/--gbwt-ri/--t1/--t2",
              file=sys.stderr)
        return
    if not hasattr(coord_index, "build_surject_anchors"):
        print("  TODO: coord_index.build_surject_anchors(...) is not yet exposed.")
        print("  Need to add a Python binding wrapping panindexer::build_surject_anchors_for_path.")
        print("  Expected signature (approx.):")
        print("    coord_index.build_surject_anchors(graph_alignment_gaf: str, target_haplotype: str)")
        print("      -> list[Anchor(step_begin_packed, step_end_packed,")
        print("                    path_offset_step_begin, path_offset_step_end,")
        print("                    read_begin, read_end,")
        print("                    gbwt_edge_begin, gbwt_edge_end)]")
        return

    # Each list entry is (anchors_list, target_path_length).
    # We use build_surject_anchors_full when available so step 3 can pass the
    # cached path length to the server (otherwise the server has to re-walk
    # the GBWT path, which is O(path_length)).
    all_anchors = []
    try:
        use_full = hasattr(coord_index, "build_surject_anchors_full")
        for gaf in graph_alignments:
            if use_full:
                res = coord_index.build_surject_anchors_full(gaf, target)
                anchors = list(res.anchors)
                path_len = res.target_path_length
                status = res.status
                print(f"  status='{status}', got {len(anchors)} anchor(s), "
                      f"target_path_length={path_len}:")
            else:
                anchors = coord_index.build_surject_anchors(gaf, target)
                path_len = 0
                print(f"  got {len(anchors)} anchor(s) for one graph alignment:")
            for a in anchors:
                print(f"    {a}")
            all_anchors.append((anchors, path_len))
    except Exception as exc:
        print(f"  ERROR building anchors: {exc}", file=sys.stderr)
        return

    if not any(a for a, _ in all_anchors):
        print("  (no anchors built — aborting before step 3)")
        return

    # Step 3 — giraffe-server: surject using pre-computed anchors. This needs
    # a SURJECT_WITH_ANCHORS stdin command in giraffe-server plus a Python
    # helper that sends the alignment + anchors and reads the response.
    print(f"\n[step 3] giraffe-server: surject with anchors")
    if not hasattr(mw, "surject_with_anchors"):
        print("  TODO: GiraffeServerMiddleware.surject_with_anchors(...) is not yet wired.")
        print("  Need: (a) SURJECT_WITH_ANCHORS stdin command on giraffe-server side,")
        print("        (b) corresponding Python helper that frames the alignment + anchors")
        print("            payload, sends it, and reads the framed GAF response.")
        print("  Expected signature (approx.):")
        print("    mw.surject_with_anchors(graph_alignment_gaf: str, anchors: list,")
        print("                            target_haplotype: str, target_path_length: int)")
        print("      -> list[str]   # GAF lines of the surjected alignment")
        return

    try:
        for gaf, (anchors, path_len) in zip(graph_alignments, all_anchors):
            if not anchors:
                continue
            surjected = mw.surject_with_anchors(
                gaf, anchors, target, target_path_length=path_len
            )
            print(f"  got {len(surjected)} surjected line(s):")
            for line in surjected:
                print(f"    {line}")
    except Exception as exc:
        print(f"  ERROR surjecting with anchors: {exc}", file=sys.stderr)
        return


def _run_interactive(mw: GiraffeServerMiddleware, coord_index=None, default_target: str = "") -> int:
    """REPL: keep indexes open, accept mapping and translate commands until EOF/quit."""
    print("Pangenome interactive session. Type 'help' for commands.", file=sys.stderr)
    if default_target:
        print(f"Default surjection target: {default_target}", file=sys.stderr)
    else:
        print("No default surjection target — graph-only mapping unless you specify one inline.", file=sys.stderr)
    if coord_index is None:
        print("Coordinate index not loaded — 'translate' and 'haplotypes' unavailable.", file=sys.stderr)
    print(file=sys.stderr)

    prompt = "> " if sys.stdin.isatty() else ""

    while True:
        try:
            if prompt:
                print(prompt, end="", flush=True, file=sys.stderr)
            line = sys.stdin.readline()
        except KeyboardInterrupt:
            print("\nInterrupted.", file=sys.stderr)
            break

        if not line:          # EOF (Ctrl-D or pipe closed)
            break

        line = line.strip()
        if not line:
            continue

        parts = line.split()
        cmd = parts[0].lower()

        if cmd in ("quit", "exit"):
            break

        if cmd == "help":
            print(HELP_TEXT)
            continue

        if cmd == "haplotypes":
            if coord_index is None:
                print("ERROR: coordinate index not loaded (pass --ri/--tags/--gbwt-ri/--t1/--t2)", file=sys.stderr)
            else:
                names = coord_index.get_haplotype_names()
                print(f"{len(names)} haplotypes:")
                for n in names:
                    print(f"  {n}")
            continue

        if cmd == "translate":
            # translate <src> <start> <end> <tgt>
            if coord_index is None:
                print("ERROR: coordinate index not loaded", file=sys.stderr)
                continue
            if len(parts) != 5:
                print("Usage: translate <src_haplotype> <start> <end> <tgt_haplotype>", file=sys.stderr)
                continue
            try:
                src, start, end, tgt = parts[1], int(parts[2]), int(parts[3]), parts[4]
                results = coord_index.translate(src, start, end, tgt)
                if not results:
                    print("(no translation results)")
                for r in results:
                    print(r)
            except Exception as exc:
                print(f"ERROR: {exc}", file=sys.stderr)
            continue

        # Otherwise treat as one of:
        #   <sequence>                          -> map only (or default target if given)
        #   <sequence> <target>                 -> map + surject (single giraffe-server call)
        #   <sequence> <target> anchor          -> anchor pipeline (3 steps)
        seq = parts[0].upper()
        anchor_mode = (len(parts) >= 2 and parts[-1].lower() == "anchor")
        if anchor_mode:
            # Strip the "anchor" keyword from the token list for target parsing.
            non_anchor_parts = parts[:-1]
            target = non_anchor_parts[1] if len(non_anchor_parts) >= 2 else default_target
            if not target:
                print("ERROR: anchor mode requires a target haplotype",
                      file=sys.stderr)
                continue
            _run_anchor_path(mw, coord_index, seq, target)
            continue

        # Regular path: map + (optional) surject in a single giraffe-server call.
        target = parts[1] if len(parts) >= 2 else default_target
        # Only include surjection target field when one is actually set;
        # an absent 4th field means graph-only mapping (avoids server hang
        # on unknown/unindexed targets).
        if target:
            read = (f"q_{seq[:8]}", seq, "I" * len(seq), target)
        else:
            read = (f"q_{seq[:8]}", seq, "I" * len(seq))
        try:
            out = mw.map_reads([read])
            mappings = out[0]
            print(f"mappings={len(mappings)} target={target}")
            for m in mappings:
                print(m)
            if not mappings:
                print("(no mappings returned)")
        except Exception as exc:
            print(f"ERROR: {exc}", file=sys.stderr)

    return 0


def main() -> int:
    args = parse_args()

    surject_target = getattr(args, "surject_target", "")
    cfg = GiraffeServerConfig(
        vg_binary=args.vg,
        gbz_path=args.gbz,
        minimizer_path=args.minimizer,
        distance_path=args.dist,
        zipcode_path=args.zipcodes,
        threads=args.threads,
        max_multimaps=args.max_multimaps,
        batch_size=args.batch_size,
        surject_target_paths=[surject_target] if surject_target else [],
    )

    mw = GiraffeServerMiddleware(cfg)
    mw.start()

    if args.interactive or args.anchors:
        # Both interactive and --anchors do per-read work and need the server
        # to be ready before the first call. (The plain smoke test doesn't
        # bother — it just waits on the first map_reads call.)
        print("Waiting for giraffe-server to finish loading indexes…", file=sys.stderr)
        mw.wait_until_ready()
        print("giraffe-server is ready.", file=sys.stderr)

    # Try to load the coordinate index if all paths are present. Needed for
    # interactive 'translate'/'haplotypes' commands and for --anchors.
    coord_index = None
    coord_args = (args.ri, args.tags, getattr(args, "gbwt_ri"), args.t1, args.t2)
    if all(coord_args):
        try:
            from middleware.pangenome_middleware import CoordinateIndexPaths, PangenomeMiddleware
            coord_paths = CoordinateIndexPaths(
                gbz_path=args.gbz,
                ri_path=args.ri,
                tags_path=args.tags,
                gbwt_ri_path=getattr(args, "gbwt_ri"),
                table1_path=args.t1,
                table2_path=args.t2,
            )
            import liftover_ext
            coord_index = liftover_ext.Index()
            coord_index.load(
                coord_paths.gbz_path,
                coord_paths.ri_path,
                coord_paths.tags_path,
                coord_paths.gbwt_ri_path,
                coord_paths.table1_path,
                coord_paths.table2_path,
            )
            print("Coordinate index loaded.", file=sys.stderr)
        except Exception as exc:
            print(f"WARNING: could not load coordinate index: {exc}", file=sys.stderr)
            coord_index = None
    elif any(coord_args):
        print("WARNING: provide all of --ri/--tags/--gbwt-ri/--t1/--t2 for coordinate translation.",
              file=sys.stderr)

    # Dispatch:
    #   --anchors        → anchor-pipeline smoke test (3-step)
    #   --interactive    → REPL session
    #   (otherwise)      → regular single-call map+surject smoke test
    if args.anchors:
        # _run_anchor_smoke_test owns mw.stop() in its finally clause.
        return _run_anchor_smoke_test(mw, coord_index, surject_target)

    if not args.interactive:
        # _run_smoke_test owns mw.stop() in its finally clause.
        return _run_smoke_test(mw, surject_target=surject_target)

    try:
        return _run_interactive(mw, coord_index, default_target=surject_target)
    finally:
        mw.stop()


if __name__ == "__main__":
    raise SystemExit(main())
