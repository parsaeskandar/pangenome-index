//
// Created by seeskand on 9/16/24.
//

#ifndef PANGENOME_INDEX_WANG_HASH_HPP
#define PANGENOME_INDEX_WANG_HASH_HPP

/// Thomas Wang's integer hash function. In many implementations, std::hash
/// is identity function for integers, which leads to performance issues.
inline size_t wang_hash_64(size_t key) {
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

#endif //PANGENOME_INDEX_WANG_HASH_HPP
