# Coordinate Translation Algorithm for Pangenome Index

## Overview

This document describes the algorithm for translating coordinates between two haplotypes (source and target) in a pangenome graph using the GBWT (Graph Burrows-Wheeler Transform) index, FastLocate r-index, and RLBWT r-index. The algorithm identifies common nodes between haplotypes and performs forward traversal to map base-level sequence positions.

## Problem Statement

**Given:**
- A pangenome graph indexed with GBWT
- A source haplotype with sequence ID `s_id` and an interval `[s_start, s_end]` (base offsets)
- A target haplotype with sequence ID `t_id`
- RLBWT r-index for finding tags and base offsets
- GBWT FastLocate r-index for finding common nodes
- Sampled tag array containing node information
- GBWTGraph for node lengths

**Find:** A mapping from each base offset in the source interval `[s_start, s_end]` to the corresponding base offset in the target haplotype.

## Key Concepts

### Coordinate Systems

1. **Node Offsets (GBWT)**: Number of nodes from the start of a path (0-indexed)
   - GBWT stores paths as sequences of nodes
   - `gbwt::FastLocate::seqOffset()` returns node offsets
   - Used for GBWT operations (LF mapping, decompressSA)

2. **Base Offsets (RLBWT)**: Number of bases from the start of a sequence
   - User input and output are in base offsets
   - Obtained from RLBWT r-index
   - Used for mapping between sequences

3. **Tag Code**: Encoded representation of a node
   - Format: `tag_code = ((node_id - 1) << 1) | is_reverse) + 1`
   - Uniquely identifies a node and its orientation
   - Used for sorting and comparison

### GBWT LF Mapping

Unlike standard FM-index LF mapping which moves backward, GBWT's LF mapping moves **forward** along paths:
- Each node's record stores outgoing edges (next nodes)
- `LF(pos)` returns the next node position along the path
- One LF operation advances one node forward (not one base)
- Efficient because: one operation = one node, good memory locality

## Algorithm Steps

### Step 1: Find Tags in Source Interval

**Input:**
- RLBWT r-index (`rlbwt_rindex`)
- Sampled tag array (`sampled`)
- Source sequence ID (`source_seq_id`)
- Source interval `[seq_start, seq_end]` (base offsets)

**Process:**
1. Convert sequence interval to packed text positions:
   ```
   text_pos_i = r_index.pack(source_seq_id, seq_start)
   text_pos_j = r_index.pack(source_seq_id, seq_end)
   ```

2. Find successor position using `last_successor(text_pos_j)` to get starting BWT position

3. Iterate backwards through the interval `[text_pos_i, text_pos_j]`:
   - For each position, get tag from sampled tag array
   - If tag exists (tag_val != 0), unpack position to get `(seq_id, offset)`
   - Store tag information: `tag_code`, `source_offsets` (base offsets from RLBWT)

4. Return vector of `TagInfo` structures containing:
   - `tag_code`: Encoded node identifier
   - `source_offsets`: List of base offsets where this tag appears in source
   - `source_bwt_positions`: BWT positions
   - `source_packed_positions`: Packed text positions

**Output:** `vector<TagInfo>` - All tags (nodes) visited by source haplotype in the interval

**Complexity:** O(k) where k is the number of positions in the interval

---

### Step 2: Find First and Last Common Nodes

**Input:**
- GBWT FastLocate r-index (`gbwt_fast_locate`)
- RLBWT r-index (`rlbwt_rindex`)
- Sampled tag array (`sampled`)
- Source tags from Step 1
- Source sequence ID (`source_seq_id`)
- Target sequence ID (`target_seq_id`)

**Process:**

#### 2.1 Sort Tags by Tag Code
```
sorted_tags = sort(source_tags, by tag_code ascending)
```

#### 2.2 Find First Common Node (from beginning)
Iterate through `sorted_tags` from start until first common node found:

For each tag `tag_info`:
1. Decode tag to get GBWT node:
   ```
   (node_id, is_rev) = decode_tag(tag_info.tag_code)
   node = gbwt::Node::encode(node_id, is_rev)
   ```

2. Use GBWT FastLocate to find all occurrences:
   ```
   sa_values = gbwt_fast_locate.decompressSA(node)
   ```
   - Returns all suffix array values (occurrences) for this node
   - Size equals number of path occurrences

3. Filter occurrences for source and target sequences:
   ```
   for each sa_value in sa_values:
       seq_id = gbwt_fast_locate.seqId(sa_value)
       seq_offset = gbwt_fast_locate.seqOffset(sa_value)  // node offset
       
       if seq_id == source_seq_id:
           source_visits.append((sa_value, seq_offset))
       if seq_id == target_seq_id:
           target_visits.append((sa_value, seq_offset))
   ```

4. If both `source_visits` and `target_visits` are non-empty:
   - **Found common node!**
   - Sort visits by offset (larger offset = earlier in GBWT path)
   - Select earliest occurrence: `source_visits[0]`, `target_visits[0]`
   - Get base offsets from RLBWT:
     - Source base: `tag_info.source_offsets[0]` (already from RLBWT)
     - Target base: Look up tag in RLBWT using `find_sequences_for_tag()`
   - Store as first common node and **break**

#### 2.3 Find Last Common Node (from end)
Iterate through `sorted_tags` from end until last common node found:

Same process as Step 2.2, but iterate in reverse order (largest tag_code first).

**Output:** `CommonNodes` structure containing:
- `first_source_offset`: GBWT node offset
- `first_target_offset`: GBWT node offset
- `first_source_base`: RLBWT base offset
- `first_target_base`: RLBWT base offset
- `first_tag_code`: Tag code
- `last_source_offset`: GBWT node offset
- `last_target_offset`: GBWT node offset
- `last_source_base`: RLBWT base offset
- `last_target_base`: RLBWT base offset
- `last_tag_code`: Tag code

**Complexity:** O(k × m) where k is number of tags checked and m is average occurrences per node. Optimized to only check until first/last found.

**Key Optimization:** Only checks tags until first common node is found (from start) and last common node is found (from end), avoiding checking all tags in between.

---

### Step 3: Forward Traversal Using LF Mapping

**Input:**
- GBWT index (`gbwt_index`)
- GBWT FastLocate (`gbwt_fast_locate`)
- GBWTGraph (`graph`) - for node lengths
- Source sequence ID (`source_seq_id`)
- Source interval `[source_start, source_end]` (base offsets)
- Target sequence ID (`target_seq_id`)
- Anchor node offsets (from Step 2)
- Last common node tag_code (from Step 2)

**Process:**

#### 3.1 Initialize Positions

1. Decode anchor node:
   ```
   (anchor_node_id, anchor_is_rev) = decode_tag(anchor_tag_code)
   anchor_node = gbwt::Node::encode(anchor_node_id, anchor_is_rev)
   ```

2. Get all occurrences of anchor node:
   ```
   sa_values = gbwt_fast_locate.decompressSA(anchor_node)
   ```

3. Find source and target positions in SA:
   ```
   for i in range(len(sa_values)):
       seq_id = gbwt_fast_locate.seqId(sa_values[i])
       seq_offset = gbwt_fast_locate.seqOffset(sa_values[i])
       
       if seq_id == source_seq_id && seq_offset == anchor_source_offset:
           source_offset_in_node = i
       if seq_id == target_seq_id && seq_offset == anchor_target_offset:
           target_offset_in_node = i
   ```

4. Create GBWT edge positions:
   ```
   source_pos = (anchor_node, source_offset_in_node)
   target_pos = (anchor_node, target_offset_in_node)
   ```

5. Initialize offsets:
   ```
   source_path_offset = anchor_source_offset  // node offset
   target_path_offset = anchor_target_offset  // node offset
   source_base_offset = anchor_source_base    // base offset from RLBWT
   target_base_offset = anchor_target_base    // base offset from RLBWT
   ```

6. Map anchor point:
   ```
   if anchor_source_base in [source_start, source_end]:
       offset_map[anchor_source_base] = anchor_target_base
   ```

#### 3.2 Forward Traversal Loop

**Important**: Paths can **diverge and converge** again. For example:
- Source path: `1 → 2 → 3`
- Target path: `1 → 4 → 5 → 6 → 2 → 3`

We traverse each path **completely independently**, store all nodes visited, then find common nodes to map positions. This ensures correct coordinate translation even when paths take different routes.

**Data Structure:**
```
struct PathNode {
    uint64_t tag_code;        // Tag code for this node
    gbwt::node_type node;     // GBWT node
    size_t node_offset;       // Node offset from start of path
    size_t base_offset;       // Base offset from start of path
    gbwt::edge_type gbwt_position;  // GBWT edge position for LF mapping
}
```

**Algorithm:**

```
// Step 1: Initialize paths with anchor node
source_path = [anchor_node]
target_path = [anchor_node]

// Step 2: Traverse source path independently
src_pos = anchor_source_position
src_node_offset = anchor_source_offset
src_base_offset = anchor_source_base

while true:
    next_src_pos = gbwt_index.LF(src_pos)
    
    if next_src_pos.first == ENDMARKER:
        break
    
    // Update base offset by adding length of node we just left
    src_base_offset += graph.get_length(src_pos.node)
    src_node_offset++
    
    // Get tag code for next node
    next_src_tag_code = encode_tag(next_src_pos.first)
    
    // Check stop conditions
    if next_src_tag_code > last_common_tag_code:
        break
    if src_base_offset > source_end:
        break
    
    // Store this node
    source_path.append(PathNode{
        tag_code: next_src_tag_code,
        node: next_src_pos.first,
        node_offset: src_node_offset,
        base_offset: src_base_offset,
        gbwt_position: next_src_pos
    })
    
    src_pos = next_src_pos

// Step 3: Traverse target path independently
tgt_pos = anchor_target_position
tgt_node_offset = anchor_target_offset
tgt_base_offset = anchor_target_base

while true:
    next_tgt_pos = gbwt_index.LF(tgt_pos)
    
    if next_tgt_pos.first == ENDMARKER:
        break
    
    // Update base offset by adding length of node we just left
    tgt_base_offset += graph.get_length(tgt_pos.node)
    tgt_node_offset++
    
    // Get tag code for next node
    next_tgt_tag_code = encode_tag(next_tgt_pos.first)
    
    // Check stop condition
    if next_tgt_tag_code > last_common_tag_code:
        break
    
    // Store this node
    target_path.append(PathNode{
        tag_code: next_tgt_tag_code,
        node: next_tgt_pos.first,
        node_offset: tgt_node_offset,
        base_offset: tgt_base_offset,
        gbwt_position: next_tgt_pos
    })
    
    tgt_pos = next_tgt_pos

// Step 4: Find common nodes and map positions
// Create index: tag_code -> list of indices in target_path
target_nodes_by_tag = {}
for i, node in enumerate(target_path):
    target_nodes_by_tag[node.tag_code].append(i)

// Iterate through source path and find matches
for src_node in source_path:
    if src_node.base_offset not in [source_start, source_end]:
        continue
    
    if src_node.tag_code in target_nodes_by_tag:
        // Found common node!
        tgt_idx = target_nodes_by_tag[src_node.tag_code][0]  // Use first occurrence
        tgt_node = target_path[tgt_idx]
        
        offset_map[src_node.base_offset] = tgt_node.base_offset
```

**Stop Conditions:**

**Source path stops when:**
1. Reaches ENDMARKER
2. Current node's `tag_code > last_common_tag_code` (passed last common node)
3. `source_base_offset > source_end` (exceeded source interval)

**Target path stops when:**
1. Reaches ENDMARKER
2. Current node's `tag_code > last_common_tag_code` (passed last common node)

**Key Insight:** By storing paths separately, we can handle divergence correctly. When paths converge again (e.g., both reach node 2), we correctly map the positions even though they took different routes to get there.

**Output:** `unordered_map<size_t, size_t>` mapping source base offsets to target base offsets

**Complexity:** O(n) where n is number of nodes between first and last common nodes

---

### Step 4: Generate Translation Results

**Input:**
- Offset map from Step 3
- Source interval `[source_start, source_end]`

**Process:**
```
for src_off in [source_start, source_end]:
    result.source_offset = src_off
    result.target_seq_id = target_seq_id
    
    if src_off in offset_map:
        result.target_offset = offset_map[src_off]
    else:
        result.target_offset = 0
    
    result.tag_code = anchor_tag_code
    translations.append(result)
```

**Output:** `vector<TranslationResult>` - Complete translation results for all positions in source interval

---

## Complete Algorithm Flow

```
1. Load all required indices:
   - RLBWT r-index
   - GBWT index
   - GBWT FastLocate r-index
   - GBWTGraph (from GBZ file)
   - Sampled tag array

2. Find tags in source interval [seq_start, seq_end]
   → Returns: vector<TagInfo> with base offsets from RLBWT

3. Find first and last common nodes:
   a. Sort tags by tag_code
   b. Search from beginning → find first common node
   c. Search from end → find last common node
   → Returns: CommonNodes with both node offsets (GBWT) and base offsets (RLBWT)

4. Forward traversal:
   a. Initialize at anchor node using GBWT positions
   b. Start with base offsets from RLBWT
   c. Traverse forward using LF mapping
   d. Accumulate node lengths to track base offsets
   e. Map positions within [source_start, source_end]
   f. Stop when passing last common node (by tag_code)
   → Returns: offset_map[source_base] = target_base

5. Generate results:
   → Returns: vector<TranslationResult> for all positions in interval
```

## Key Implementation Details

### Tag Code Encoding/Decoding

**Encoding:**
```cpp
uint64_t decoded = ((node_id - 1) << 1) | (is_rev ? 1 : 0);
uint64_t tag_code = decoded + 1;
```

**Decoding:**
```cpp
uint64_t decoded = tag_code - 1;
int64_t node_id = (decoded >> 1) + 1;
bool is_rev = (decoded & 1) != 0;
```

### Base Offset Tracking

During traversal:
- Start with base offsets from RLBWT at anchor node
- For each node traversed:
  ```
  source_base_offset += graph.get_length(current_node_handle)
  target_base_offset += graph.get_length(current_node_handle)
  ```
- This accumulates node lengths to convert node offsets to base offsets

### Stop Condition

The algorithm stops when:
1. **Primary**: Either `source_tag_code > last_common_tag_code` OR `target_tag_code > last_common_tag_code`
   - Checks BOTH source and target paths independently
   - Uses tag_code comparison (reliable, node-independent)
   - Does NOT use node_offset comparison (unreliable, path-dependent)
   - Stops if EITHER path passes the last common node

2. **Secondary**: `source_base_offset > source_end`
   - Safety check to avoid exceeding input interval

3. **Path end**: Either path reaches ENDMARKER
   - Checked independently for each path

## Example Execution

**Input:**
- Source sequence ID: 0
- Source interval: [10, 50] (base offsets)
- Target sequence ID: 2

**Step 1:** Find tags in interval [10, 50]
- Finds tags: [5, 7, 9, 15, 23, 29] (tag_codes)

**Step 2:** Find common nodes
- First common node: tag_code=5, source_base=10, target_base=12
- Last common node: tag_code=29, source_base=45, target_base=47

**Step 3:** Traverse from anchor
- Start: node_offset=143 (source), base_offset=10
- Step 1: node_offset=144, base_offset=13 → map(13, 15)
- Step 2: node_offset=145, base_offset=18 → map(18, 20)
- ...
- Step N: node_offset=150, base_offset=45 → map(45, 47)
- Step N+1: tag_code=31 > 29 → **stop**

**Output:**
```
Source_Offset  Target_Seq_ID  Target_Offset
10             2              12
13             2              15
18             2              20
...
45             2              47
```

## Complexity Analysis

| Step | Operation | Complexity |
|------|-----------|------------|
| 1. Find tags | Iterate interval | O(k) where k = interval size |
| 2. Find common nodes | Check tags until found | O(k' × m) where k' = tags checked, m = avg occurrences |
| 3. Traversal | LF mapping forward | O(n) where n = nodes between first/last |
| 4. Generate results | Iterate interval | O(k) |

**Overall:** O(k + k'×m + n) where typically k' << k (optimization reduces checks)

## Notes

- **Node offsets vs Base offsets**: GBWT uses node offsets internally, but user input/output uses base offsets. Conversion happens via node lengths from GBWTGraph.

- **Tag code ordering**: Tags are sorted by tag_code, which corresponds to node ID ordering. This allows efficient binary search-like behavior.

- **RLBWT for base offsets**: Base offsets come from RLBWT r-index, not calculated from path start. This is more efficient and accurate.

- **Stop condition**: Uses tag_code comparison, not node_offset, because node offsets can differ between paths for the same logical position.

