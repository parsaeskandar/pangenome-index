
//
// Created by seeskand on 5/22/24.
/*
   This implementation of the B+ tree is capable handling the Run Length encoding of a string.
   The Runs having a starting position as the keys and a graph position as the value. The B+ tree insert function gets
   Run and Run length as the input and inserts the run into the tree and it also add a "gap" run which is used to keep
   the length of the actual run. For example adding a Run(1, G) with length 5 will add two runs to the tree, one being
    Run(1, G) and the other being Run(6, 0) where 0 is the graph position of the gap run (a special value). The tree
    is capable of merging, splitting and handling the insertion of the runs.

*/
#ifndef PANGENOME_INDEX_BPLUS_TREE_HPP
#define PANGENOME_INDEX_BPLUS_TREE_HPP

#include <iostream>
#include <gbwtgraph/gbwtgraph.h>
//#include "../gbwtgraph_helper.hpp"
#include <gbwtgraph/index.h>
#include <vector>

// Create a Run struct to store the run data structure for the BPlusTree
struct Run {
    size_t start_position;
    gbwtgraph::Position graph_position;

    bool operator<(const Run &other) const;
    bool operator==(const Run &other) const;
    Run &operator=(const Run &other);
    Run &operator=(int zero);
    friend std::ostream &operator<<(std::ostream &os, const Run &run);
};

template<typename T>
class bpNode {
public:
    bpNode(std::size_t degree, bool leaf = false);
    ~bpNode();

    void set_next(bpNode<T> *next_node);
    void set_prev(bpNode<T> *prev_node);
    void set_parent(bpNode<T> *parent_node);
    void add_child(bpNode<T> *child);
    void add_item(T item);
    void add_item(T item, size_t i);
    void add_child(bpNode<T> *child, size_t i);
    void change_item(T new_item, size_t i);
    size_t search(T &data);
    bpNode<T> *get_parent();
    int find_index_in_parent();
    bool is_underflowing() const;
    void remove_item(size_t i);
    void remove_child(size_t i);
    bool is_full() const;
    bool is_leaf() const;
    bool is_root() const;
    bool is_empty() const;
    std::vector<T> get_items();
    std::vector<bpNode<T> *> get_children();
    bpNode<T> *get_next();
    bpNode<T> *get_prev();
    void remove_first_item();
    void remove_last_item();
    void remove_last_child();
    int get_size();
    bpNode<T> *get_child(size_t i);
    T get_item(size_t i);
    void replace_items_leaf(std::vector<T> new_items);
    void replace_children(std::vector<bpNode<T> *> new_children);
    std::vector<T> insert(const T &data, size_t run_length);
    std::vector<T> insert_success(const T &data, size_t run_length, bool &success);
    int search_non_leaf(const T &data);
    void print();
    bool is_gap(T run);

private:
    bool leaf;
    size_t degree;
    std::vector<T> items;
    std::vector<bpNode<T> *> children;
    bpNode<T> *parent;
    bpNode<T> *next;
    bpNode<T> *prev;
    T create_gap(size_t gap_start_position);
    std::vector<T> run_insert(std::vector<T> new_items, T new_run, int run_length, size_t &inserted_index);
    bool merge_item_next(std::vector<T> &new_items);
    bool merge_item_prev(std::vector<T> &new_items);
};

template<typename T>
class BplusTree {
public:
    BplusTree(std::size_t _degree);
    ~BplusTree();
    bpNode<T> *get_root();
    bpNode<T> *leaf_search(bpNode<T> *node, T data);
    size_t *search_index(T data);
    T search(size_t position);
    void insert(const T &data, size_t run_length);
    bool insert_success(const T &data, size_t run_length);
    void print(bpNode<T> *node);
    void print_whole();
    size_t get_bpt_size();
    size_t get_bpt_run_count();

    class Iterator {
        bpNode<T> *node;
        int index;

    public:
        Iterator(bpNode<T> *node, int index);
        T operator*() const;
        Iterator &operator++();
        Iterator &operator--();
        bool operator!=(const Iterator &other) const;
        bool operator==(const Iterator &other) const;
    };

    Iterator begin();
    Iterator end();

private:
    bpNode<T> *root;
    std::size_t degree;
    void parent_insert(bpNode<T> *node, T data, bpNode<T> *child);
    void remove_from_parent(bpNode<T> *node);
    void handle_parent_underflow(bpNode<T> *node);
    void leaf_underflow(bpNode<T> *node);
    size_t get_size(bpNode<T> *cursor);
    size_t get_run_count(bpNode<T> *cursor);
};

#endif // PANGENOME_INDEX_BPLUS_TREE_HPP
