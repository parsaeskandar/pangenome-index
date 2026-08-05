# Coordinate Translation Algorithm - Toy Example

This document explains the coordinate translation algorithm using a simple toy graph example.

---

## 1. Problem Statement

**Goal:** Given a position interval `[start, end]` on a **source haplotype**, find the corresponding positions on a **target haplotype** by identifying shared nodes (anchors) and tracing paths using GBWT LF-mapping.

---

## 2. Toy Graph Structure

Consider a simple pangenome graph with 6 nodes:

```
                    ┌─────┐
                    │  3  │
                    │ CAT │
                    └──┬──┘
                       │
    ┌─────┐    ┌─────┐ │ ┌─────┐    ┌─────┐
    │  1  │───▶│  2  │─┼▶│  4  │───▶│  6  │
    │ ATG │    │ CGA │ │ │ TTA │    │ GGC │
    └─────┘    └─────┘ │ └─────┘    └─────┘
                       │
                    ┌──┴──┐
                    │  5  │
                    │ AAA │
                    └─────┘
```

**Node sequences and lengths:**
| Node ID | Sequence | Length (bp) |
|---------|----------|-------------|
| 1       | ATG      | 3           |
| 2       | CGA      | 3           |
| 3       | CAT      | 3           |
| 4       | TTA      | 3           |
| 5       | AAA      | 3           |
| 6       | GGC      | 3           |

---

## 3. Haplotype Paths

Two haplotypes traverse this graph differently:

### Source Haplotype (seq_id = 0)
```
Path:  1 → 2 → 3 → 4 → 6
Seq:   ATG-CGA-CAT-TTA-GGC
       └─┬─┘└─┬─┘└─┬─┘└─┬─┘└─┬─┘
Base:  0-2  3-5  6-8  9-11 12-14
```

### Target Haplotype (seq_id = 1)
```
Path:  1 → 2 → 5 → 4 → 6
Seq:   ATG-CGA-AAA-TTA-GGC
       └─┬─┘└─┬─┘└─┬─┘└─┬─┘└─┬─┘
Base:  0-2  3-5  6-8  9-11 12-14
```

**Key observation:** Nodes 1, 2, 4, and 6 are **common** to both haplotypes. Node 3 is unique to source, and node 5 is unique to target.

---

## 4. Algorithm Walkthrough

### Input
- **Source haplotype:** seq_id = 0
- **Source interval:** bases [3, 11] (covers nodes 2, 3, 4)
- **Target haplotype:** seq_id = 1

### Step 1: Find Tags in Source Interval

Using the RLBWT r-index and sampled tag array, find all nodes visited by the source haplotype in the interval `[3, 11]`:

```
Tags found in source interval:
┌───────────┬──────────────┬───────────────┐
│ Tag Code  │ Node ID      │ Source Offset │
├───────────┼──────────────┼───────────────┤
│ 3         │ 2            │ 3             │
│ 5         │ 3            │ 6             │
│ 7         │ 4            │ 9             │
└───────────┴──────────────┴───────────────┘
```

**Tag encoding:** `tag_code = 1 + ((node_id - 1) << 1) | is_reverse`

### Step 2: Find First and Last Common Nodes

Sort tags by tag_code and search for nodes that exist in **both** source and target:

```
Checking tag_code=3 (node 2): ✓ Found in both!
  → FIRST common node
  → Source: base_offset=3, node_offset=1
  → Target: base_offset=3, node_offset=1

Checking tag_code=5 (node 3): ✗ Not in target

Checking tag_code=7 (node 4): ✓ Found in both!
  → LAST common node
  → Source: base_offset=9, node_offset=3
  → Target: base_offset=9, node_offset=3
```

### Step 3: Build Offset Mapping using GBWT LF

Starting from the **anchor** (first common node), traverse both paths forward using GBWT LF-mapping:

```
GBWT LF Traversal (forward along path):

Source Path Traversal:
  Anchor (node 2) → base_offset=3
  LF step → node 3 → base_offset=6  [NOT in target - skip]
  LF step → node 4 → base_offset=9  [COMMON - map!]
  LF step → node 6 → base_offset=12 [beyond interval - stop]

Target Path Traversal:
  Anchor (node 2) → base_offset=3
  LF step → node 5 → base_offset=6  [NOT in source]
  LF step → node 4 → base_offset=9  [COMMON]
  LF step → node 6 → base_offset=12
```

### Step 4: Create Position Mappings

Only positions at **common nodes** can be directly mapped:

```
┌─────────────────┬─────────────────┬────────────┐
│ Source Offset   │ Target Offset   │ Common Node│
├─────────────────┼─────────────────┼────────────┤
│ 3               │ 3               │ Node 2     │
│ 9               │ 9               │ Node 4     │
└─────────────────┴─────────────────┴────────────┘
```

**Note:** Source offset 6 (node 3) has **no direct mapping** because node 3 is not visited by the target haplotype.

---

## 5. Visual Summary

```
Source (seq_id=0):
  [0-2]   [3-5]   [6-8]   [9-11]  [12-14]
    1   →   2   →   3   →   4   →   6
    │       │       ╳       │       │
    │       │   no mapping  │       │
    ▼       ▼               ▼       ▼
    1   →   2   →   5   →   4   →   6
  [0-2]   [3-5]   [6-8]   [9-11]  [12-14]
Target (seq_id=1):

Query interval: [3, 11]

Mapped positions:
  Source 3  ──────────────▶  Target 3   (at node 2)
  Source 6  ──────────────▶  NO MAPPING (node 3 unique to source)
  Source 9  ──────────────▶  Target 9   (at node 4)
```

---

## 6. Key Concepts

### Tag Encoding
```
tag_code = 1 + ((node_id - 1) << 1) | is_reverse
```
- Tag 0 is reserved for gaps
- Even bits encode forward orientation, odd bits encode reverse

### GBWT LF-Mapping
- Unlike standard FM-index LF (moves backward), GBWT's LF moves **forward** along the path
- Each node's record stores outgoing edges (next nodes)
- `LF(pos)` returns the next node position along the path

### Finding Common Nodes
1. For each tag in source interval, use `decompressSA(node)` to find all paths visiting that node
2. Check if target sequence appears in the results
3. First common node = anchor point for alignment
4. Last common node = stopping point for traversal

### Offset Types
- **Base offset (RLBWT):** Position in bases from start of haplotype sequence
- **Node offset (GBWT):** Position in nodes from start of path (larger = earlier)

---

## 7. Algorithm Complexity

| Operation | Complexity |
|-----------|------------|
| Find tags in interval | O(interval_size) |
| Find common nodes | O(tags × path_occurrences) |
| GBWT LF traversal | O(nodes_in_interval) |
| Total | O(interval_size + common_nodes) |

---

## 8. Edge Cases

1. **No common nodes:** Source and target paths don't share any nodes in the interval
2. **Multiple occurrences:** A node may appear multiple times in a path (cycles/loops)
3. **Diverge and reconverge:** Paths may share a node, diverge, then share another node later
4. **Fragments:** Target may be split across multiple sequence IDs

---

## 9. Data Structures Used

| Structure | Purpose |
|-----------|---------|
| RLBWT R-index | Find text positions and BWT positions for base-level coordinates |
| GBWT Index | Store haplotype paths as node sequences |
| GBWT FastLocate | Efficiently decompress SA values for specific nodes |
| Sampled Tag Array | Map BWT positions to graph node tags |
| GBWTGraph | Get node lengths for base offset calculations |

---

## 10. Example Output

```
Coordinate Translation Results:
Source haplotype ID: 0
Source interval: 3..11
Target haplotype ID: 1
First common node (anchor): tag_code=3
  Source: node_offset=1, base_offset=3
  Target: node_offset=1, base_offset=3
Last common node (stop): tag_code=7
  Source: node_offset=3, base_offset=9

Translations:
Source_Offset	Target_Seq_ID	Target_Offset
3	1	3
6	1	0   ← No mapping (node 3 not in target)
9	1	9
```
