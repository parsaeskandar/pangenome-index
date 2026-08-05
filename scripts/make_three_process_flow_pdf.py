#!/usr/bin/env python3
"""
Generate an editable vector PDF of the three-process flow diagram.

The PDF keeps text as text (TrueType, fonttype=42) and shapes as vectors,
so it opens in Illustrator with every box, arrow, and label individually
selectable and editable.

Run:
    python3 scripts/make_three_process_flow_pdf.py
Output:
    three_process_flow.pdf  (in the repo root)
"""
from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Critical for Illustrator editability:
#   fonttype = 42 -> embed TrueType, text stays as text glyphs
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["pdf.use14corefonts"] = False
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]

# Colors (soft pastels, good contrast for slides)
MIDDLEWARE_FACE = "#FFF1D6"
MIDDLEWARE_EDGE = "#B8860B"
GIRAFFE_FACE    = "#DCEEFB"
GIRAFFE_EDGE    = "#1F6FA8"
COORD_FACE      = "#DCF1E0"
COORD_EDGE      = "#2E7D32"
ARROW_REQ       = "#1F6FA8"   # request direction
ARROW_REP       = "#A8281F"   # reply direction


def add_box(ax, x, y, w, h, title, subtitle, bullets, face, edge):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.08,rounding_size=0.18",
        linewidth=1.8, edgecolor=edge, facecolor=face,
    )
    ax.add_patch(box)
    cx = x + w / 2
    ax.text(cx, y + h - 0.45, title, ha="center", va="center",
            fontsize=14, fontweight="bold", color=edge)
    if subtitle:
        ax.text(cx, y + h - 0.90, subtitle, ha="center", va="center",
                fontsize=10, style="italic", color="#444444")
    line_y = y + h - 1.45
    for line in bullets:
        ax.text(x + 0.30, line_y, line, ha="left", va="center",
                fontsize=10.5, color="#202020")
        line_y -= 0.38


def add_arrow(ax, x1, y1, x2, y2, color, label,
              label_pos=None, label_color=None):
    """label_pos: (x, y) absolute coordinates for the label. If None, midpoint."""
    arr = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>", mutation_scale=18,
        linewidth=1.6, color=color, shrinkA=0, shrinkB=0,
    )
    ax.add_patch(arr)
    if label_pos is None:
        lx = (x1 + x2) / 2
        ly = (y1 + y2) / 2
    else:
        lx, ly = label_pos
    ax.text(lx, ly, label, ha="center", va="center", fontsize=9.5,
            color=(label_color or color),
            bbox=dict(facecolor="white", edgecolor="none", pad=1.5))


def build_figure():
    fig, ax = plt.subplots(figsize=(14.5, 9.0))
    ax.set_xlim(0, 14.5)
    ax.set_ylim(0, 9)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- Middleware (top center) ----
    mw_x, mw_y, mw_w, mw_h = 4.5, 7.0, 5.5, 1.6
    add_box(
        ax, mw_x, mw_y, mw_w, mw_h,
        title="middleware (Python)",
        subtitle="orchestrates the pipeline",
        bullets=[],
        face=MIDDLEWARE_FACE, edge=MIDDLEWARE_EDGE,
    )

    # ---- giraffe-server (bottom left) ----
    gs_x, gs_y, gs_w, gs_h = 0.4, 1.0, 5.6, 4.6
    add_box(
        ax, gs_x, gs_y, gs_w, gs_h,
        title="giraffe-server",
        subtitle="(always open)",
        bullets=[
            "• GBZ",
            "• mapper",
            "• surjector",
            "      ─ map()",
            "      ─ surject_with_anchors()",
        ],
        face=GIRAFFE_FACE, edge=GIRAFFE_EDGE,
    )

    # ---- coord-trans-server (bottom right) ----
    ct_x, ct_y, ct_w, ct_h = 8.5, 1.0, 5.6, 4.6
    add_box(
        ax, ct_x, ct_y, ct_w, ct_h,
        title="coord-trans-server",
        subtitle="(always open)",
        bullets=[
            "• GBZ",
            "• FastLocate",
            "• r-index",
            "• tag array",
            "• translation tables",
            "• build_anchors()",
        ],
        face=COORD_FACE, edge=COORD_EDGE,
    )

    # ---- Arrows ----
    # Middleware bottom edge: y = mw_y = 7.0
    # Giraffe top edge: y = gs_y + gs_h = 5.6
    # Coord top edge:   y = ct_y + ct_h = 5.6

    # Arrows between middleware and giraffe-server: 4 arrows fan across the
    # box's top edge. Labels are stacked at four DIFFERENT y positions so they
    # never collide regardless of how close the arrow x positions are.
    g_top_y = gs_y + gs_h
    mw_bot_y = mw_y
    # vertical band between the two boxes: 5.6 → 7.0 (1.4 units of space)
    y_label_4 = 6.78  # near middleware (top)
    y_label_3 = 6.49
    y_label_2 = 6.20
    y_label_1 = 5.91  # near giraffe (bottom)

    # Step 1 (request): middleware → giraffe
    add_arrow(ax, 1.2, mw_bot_y, 1.2, g_top_y + 0.05,
              ARROW_REQ, "1. read", label_pos=(1.2, y_label_4))
    # Step 2 (reply): giraffe → middleware
    add_arrow(ax, 2.4, g_top_y + 0.05, 2.4, mw_bot_y,
              ARROW_REP, "2. alignment", label_pos=(2.4, y_label_3))
    # Step 5 (request): middleware → giraffe (alignment + anchors)
    add_arrow(ax, 4.0, mw_bot_y, 4.0, g_top_y + 0.05,
              ARROW_REQ, "5. aln + anchors", label_pos=(4.0, y_label_2))
    # Step 6 (reply): giraffe → middleware (surjected)
    add_arrow(ax, 5.2, g_top_y + 0.05, 5.2, mw_bot_y,
              ARROW_REP, "6. surjected", label_pos=(5.2, y_label_1))

    # Arrows between middleware and coord-trans: 2 arrows, plenty of room.
    c_top_y = ct_y + ct_h
    add_arrow(ax, 10.6, mw_bot_y, 10.6, c_top_y + 0.05,
              ARROW_REQ, "3. aln + target", label_pos=(10.6, y_label_3))
    add_arrow(ax, 11.9, c_top_y + 0.05, 11.9, mw_bot_y,
              ARROW_REP, "4. anchors", label_pos=(11.9, y_label_2))

    # ---- Title ----
    ax.text(7.25, 8.85,
            "Three-process flow — middleware orchestrates giraffe-server + coord-trans-server",
            ha="center", va="center", fontsize=14, fontweight="bold",
            color="#202020")

    # ---- Legend ----
    legend_x, legend_y = 0.4, 0.30
    ax.add_patch(FancyBboxPatch((legend_x, legend_y), 6.5, 0.55,
                                 boxstyle="round,pad=0.04",
                                 linewidth=0.8, edgecolor="#999999",
                                 facecolor="white"))
    # request swatch
    ax.add_patch(FancyArrowPatch((legend_x + 0.2, legend_y + 0.28),
                                  (legend_x + 0.85, legend_y + 0.28),
                                  arrowstyle="-|>", mutation_scale=14,
                                  linewidth=1.6, color=ARROW_REQ))
    ax.text(legend_x + 1.0, legend_y + 0.28, "request",
            ha="left", va="center", fontsize=10, color="#202020")
    # reply swatch
    ax.add_patch(FancyArrowPatch((legend_x + 2.4, legend_y + 0.28),
                                  (legend_x + 3.05, legend_y + 0.28),
                                  arrowstyle="-|>", mutation_scale=14,
                                  linewidth=1.6, color=ARROW_REP))
    ax.text(legend_x + 3.2, legend_y + 0.28, "reply",
            ha="left", va="center", fontsize=10, color="#202020")
    ax.text(legend_x + 4.3, legend_y + 0.28,
            "numbers = sequential step order",
            ha="left", va="center", fontsize=9.5, style="italic",
            color="#666666")

    return fig


def main() -> None:
    fig = build_figure()
    out = "three_process_flow.pdf"
    fig.savefig(out, format="pdf", bbox_inches="tight")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
