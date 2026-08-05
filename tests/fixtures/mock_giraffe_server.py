#!/usr/bin/env python3
"""A mock `vg giraffe-server` for testing the middleware multiplexer.

Speaks just enough of the framed-output protocol to exercise the dispatcher +
demuxer without needing a real vg or any indexes:

  * a read line ``name<TAB>seq[<TAB>qual[<TAB>target]]`` is buffered;
  * ``PROCESS_BATCH`` / ``FLUSH_NOW`` emits, for each buffered read, a frame
    ``READ<TAB>name<TAB>1`` followed by ``name<TAB>gaf:<seq>`` — echoing the
    sequence so a test can verify each frame reached the right caller;
  * ``SURJECT_WITH_ANCHORS<TAB>name<TAB>target<TAB>path_len<TAB>n`` flushes any
    pending reads (like the real server), consumes ``n`` anchor lines + one GAF
    line, and emits ``READ<TAB>name<TAB>1`` + ``name<TAB>surj:<gaf>``.

All command-line arguments are ignored. Reads stdin until EOF.
"""
import sys


def main() -> None:
    out = sys.stdout
    stdin = sys.stdin
    buffered = []  # list of (name, seq)

    def flush_batch() -> None:
        for name, seq in buffered:
            out.write(f"READ\t{name}\t1\n")
            out.write(f"{name}\tgaf:{seq}\n")
        buffered.clear()
        out.flush()

    while True:
        line = stdin.readline()
        if not line:                      # EOF: stdin closed by the client
            break
        line = line.rstrip("\n")

        if line in ("PROCESS_BATCH", "FLUSH_NOW"):
            flush_batch()
            continue

        if line.startswith("SURJECT_WITH_ANCHORS\t"):
            flush_batch()                 # keep output order, like the real server
            fields = line.split("\t")
            name = fields[1] if len(fields) > 1 else "?"
            n_anchors = int(fields[4]) if len(fields) > 4 else 0
            for _ in range(n_anchors):
                stdin.readline()          # discard anchor lines
            gaf = stdin.readline().rstrip("\n")
            out.write(f"READ\t{name}\t1\n")
            out.write(f"{name}\tsurj:{gaf}\n")
            out.flush()
            continue

        # Otherwise a read line.
        parts = line.split("\t")
        if len(parts) >= 2:
            buffered.append((parts[0], parts[1]))


if __name__ == "__main__":
    main()
