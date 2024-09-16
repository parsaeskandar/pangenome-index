//
// Created by seeskand on 9/16/24.
//

#ifndef PANGENOME_INDEX_UNIQUE_KMER_HPP
#define PANGENOME_INDEX_UNIQUE_KMER_HPP

#include <gbwtgraph/index.h>
#include <omp.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/utils.h>
#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/minimizer.h>
#include <hash_map.hpp>


namespace gbwtgraph {
    using namespace std;

// Function to extract unique kmers from the graph in parallel
    template<class KeyType>
    void
    unique_kmers_parallel(const GBWTGraph &graph,
                          hash_map <gbwtgraph::Key64::value_type, gbwtgraph::Position> &index,
                          size_t k);


}

#endif //PANGENOME_INDEX_UNIQUE_KMER_HPP
