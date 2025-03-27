
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
//
//#include <iostream>
//#include <gbwtgraph/gbwtgraph.h>
////#include "../gbwtgraph_helper.hpp"
//#include <gbwtgraph/index.h>
//#include <vector>
//
//// Create a Run struct to store the run data structure for the BPlusTree
//struct Run {
//    size_t start_position;
//    gbwtgraph::Position graph_position;
//
//    bool operator<(const Run &other) const;
//    bool operator==(const Run &other) const;
//    Run &operator=(const Run &other);
//    Run &operator=(int zero);
//    friend std::ostream &operator<<(std::ostream &os, const Run &run);
//};
//
//template<typename T>
//class bpNode {
//public:
//    bpNode(std::size_t degree, bool leaf = false);
//    ~bpNode();
//
//    void set_next(bpNode<T> *next_node);
//    void set_prev(bpNode<T> *prev_node);
//    void set_parent(bpNode<T> *parent_node);
//    void add_child(bpNode<T> *child);
//    void add_item(T item);
//    void add_item(T item, size_t i);
//    void add_child(bpNode<T> *child, size_t i);
//    void change_item(T new_item, size_t i);
//    size_t search(T &data);
//    bpNode<T> *get_parent();
//    int find_index_in_parent();
//    bool is_underflowing() const;
//    void remove_item(size_t i);
//    void remove_child(size_t i);
//    bool is_full() const;
//    bool is_leaf() const;
//    bool is_root() const;
//    bool is_empty() const;
//    std::vector<T> get_items();
//    std::vector<bpNode<T> *> get_children();
//    bpNode<T> *get_next();
//    bpNode<T> *get_prev();
//    void remove_first_item();
//    void remove_last_item();
//    void remove_last_child();
//    int get_size();
//    bpNode<T> *get_child(size_t i);
//    T get_item(size_t i);
//    void replace_items_leaf(std::vector<T> new_items);
//    void replace_children(std::vector<bpNode<T> *> new_children);
//    std::vector<T> insert(const T &data, size_t run_length);
//    std::vector<T> insert_success(const T &data, size_t run_length, bool &success);
//    int search_non_leaf(const T &data);
//    void print();
//    bool is_gap(T run);
//
//private:
//    bool leaf;
//    size_t degree;
//    std::vector<T> items;
//    std::vector<bpNode<T> *> children;
//    bpNode<T> *parent;
//    bpNode<T> *next;
//    bpNode<T> *prev;
//    T create_gap(size_t gap_start_position);
//    std::vector<T> run_insert(std::vector<T> new_items, T new_run, int run_length, size_t &inserted_index);
//    bool merge_item_next(std::vector<T> &new_items);
//    bool merge_item_prev(std::vector<T> &new_items);
//};
//
//template<typename T>
//class BplusTree {
//public:
//    BplusTree(std::size_t _degree);
//    ~BplusTree();
//    bpNode<T> *get_root();
//    bpNode<T> *leaf_search(bpNode<T> *node, T data);
//    size_t *search_index(T data);
//    T search(size_t position);
//    void insert(const T &data, size_t run_length);
//    bool insert_success(const T &data, size_t run_length);
//    void print(bpNode<T> *node);
//    void print_whole();
////    size_t get_bpt_size();
//    size_t get_bpt_run_count();
//
//    class Iterator {
//        bpNode<T> *node;
//        int index;
//
//    public:
//        Iterator(bpNode<T> *node, int index);
//        T operator*() const;
//        Iterator &operator++();
//        Iterator &operator--();
//        bool operator!=(const Iterator &other) const;
//        bool operator==(const Iterator &other) const;
//    };
//
//    Iterator begin();
//    Iterator end();
//
//    size_t get_bpt_size() {
//        return get_size(this->root);
//    }
//
//    size_t get_size(bpNode <T> *cursor) {
//        size_t size = 0;
//        if (cursor != NULL) {
//            if (cursor->is_leaf()) {
//                size += cursor->get_size();
//            } else {
//                for (int i = 0; i < cursor->get_size() + 1; ++i) {
//                    size += get_size(cursor->get_child(i));
//                }
//            }
//        }
//        return size;
//    }
//
//private:
//    bpNode<T> *root;
//    std::size_t degree;
//    void parent_insert(bpNode<T> *node, T data, bpNode<T> *child);
//    void remove_from_parent(bpNode<T> *node);
//    void handle_parent_underflow(bpNode<T> *node);
//    void leaf_underflow(bpNode<T> *node);
////    size_t get_size(bpNode<T> *cursor);
//    size_t get_run_count(bpNode<T> *cursor);
//
//};


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

#include <iostream>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/index.h>

using namespace std;

namespace panindexer {


// Create a Run struct to store the run data structure for the BPlusTree
    struct Run {
        size_t start_position;
        gbwtgraph::Position graph_position;

        // Operators for the struct
        bool operator<(const Run &other) const {
            return start_position < other.start_position;
        }

//    bool operator==(const Run &other) const {
//        return start_position == other.start_position;
//    }

        bool operator==(const Run &other) const {
            return (start_position == other.start_position && graph_position.value == other.graph_position.value);
        }


        Run &operator=(const Run &other) {
            if (this == &other) {
                return *this;
            }
            start_position = other.start_position;
            graph_position = other.graph_position;
//        run_char = other.run_char;
            return *this;
        }

        // Set the Run struct to a zero one when it is assigned to 0
        Run &operator=(int zero) {
            if (zero == 0) {  // Ensures it only responds to 0
                start_position = 0;
                graph_position = gbwtgraph::Position::no_value();
            }
            return *this;
        }

        // print the Run struct
        friend std::ostream &operator<<(std::ostream &os, const Run &run) {
            os << "start_position: " << run.start_position
               << ", graph_position: " << run.graph_position.value;

            return os;
        }

    };

    template<typename T>
    class bpNode {
    public:
        bpNode(std::size_t degree, bool leaf = false) : leaf(leaf), degree(degree) {
            items.reserve(degree);
            children.reserve(degree + 1);
            parent = nullptr;
            next = nullptr;
            prev = nullptr;
//        current_size = 0;
        }

        // destructor
        ~bpNode() {
            children.clear();
            items.clear();
        }

        bool operator==(const bpNode<T>& other) const {
            return (degree == other.degree) &&
                   (leaf == other.leaf) &&
                   (items == other.items) &&
                   (children == other.children) &&
                   (parent == other.parent) &&
                   (next == other.next) &&
                   (prev == other.prev);
        }

        // set next
        void set_next(bpNode<T> *next_node) {
            next = next_node;
        }

        // set prev
        void set_prev(bpNode<T> *prev_node) {
            prev = prev_node;
        }

        // set parent
        void set_parent(bpNode<T> *parent_node) {
            parent = parent_node;
        }

        // add child
        void add_child(bpNode<T> *child) {
            children.push_back(child);
        }

        // add item
        void add_item(T item) {
            items.push_back(item);
        }

        // add item in ith position
        void add_item(T item, size_t i) {
            items.insert(items.begin() + i, item);
        }

        // add child in ith position
        void add_child(bpNode<T> *child, size_t i) {
            children.insert(children.begin() + i, child);
        }

        // change ith item to new item
        void change_item(T new_item, size_t i) {
            items[i] = new_item;

//        if (parent != nullptr && i == 0) {
//            auto index_parent = this->find_index_in_parent();
//            cout << "index parent " << index_parent << " " << items[0] << endl;
//            if (index_parent != 0) {
//                parent->change_item(items[0], index_parent - 1);
//            }
//        }


        }


        // find the index of the item in the node
        size_t search(T &data) {
            for (size_t i = 0; i < items.size(); i++) {
                if (data < items[i]) {
                    return i;
                }
                if (i == items.size() - 1) {
                    return i + 1;
                }
            }
            return -1;
        }

        // get parent node
        bpNode<T> *get_parent() {
            return parent;
        }

        // find the index of this node in the parent
        int find_index_in_parent() {
            assert(parent != nullptr);

            for (int i = 0; i < parent->get_size() + 1; i++) {
                if (parent->get_child(i) == this) {
                    return i;
                }
            }

            throw std::runtime_error("Node not found in its parent.");
        }

        // check if it is underflowing
        bool is_underflowing() const {
            if (leaf) {
                return items.size() < ((degree + 1) / 2 - 1);
            } else {
                return children.size() < ((degree + 1) / 2 - 1);
            }
        }

        // remove the item at ith position
        void remove_item(size_t i) {
            items.erase(items.begin() + i);
        }

        // remove the child at ith position
        void remove_child(size_t i) {
            children.erase(children.begin() + i);
        }


        // check if the node is full
        bool is_full() const {
            return items.size() == items.capacity();
        }

        bool is_leaf() const {
            return leaf;
        }

        bool is_root() const {
            return parent == nullptr;
        }

        bool is_empty() const {
            return items.empty();
        }

        // get items
        std::vector <T> get_items() {
            return items;
        }

        // get children
        std::vector<bpNode<T> *> get_children() {
            return children;
        }

        // get next node
        bpNode<T> *get_next() {
            return next;
        }

        // get previous node
        bpNode<T> *get_prev() {
            return prev;
        }

        void remove_first_item() {
            items.erase(items.begin());

            // change the parent corresponding item
//        if (parent != nullptr) {
//            auto index_parent = this->find_index_in_parent();
//            cout << "index parent " << index_parent << " " << items[0] << endl;
//            if (index_parent != 0) {
//                parent->change_item(items[0], index_parent - 1);
//            }
//        }
        }

        void remove_last_item() {
            items.pop_back();
        }

        // remove last children
        void remove_last_child() {
            children.pop_back();
        }

        int get_size() {
            return items.size();
        }

        // get ith children
        bpNode<T> *get_child(size_t i) {
            bpNode<T> *child = children[i];
            return child;
        }

        T get_item(size_t i) {
            return items[i];
        }


        // replace the items in the leaf node with new items
        void replace_items_leaf(vector <T> new_items) {
//        if (DEEBUG) cout << "replacing items in leaf node" << endl;
            items.assign(new_items.begin(), new_items.end());

        }

        // replace children
        void replace_children(vector<bpNode<T> *> new_children) {
            children.assign(new_children.begin(), new_children.end());
        }


        /*
           This function handle the insertion of a new run into the node. It will handle the splitting and merging of the
             nodes if needed.

        */
        vector <T> insert(const T &data, size_t run_length) {
//        if (DEEBUG) cout << "inserting into node " << data << " " << run_length << endl;
            assert(leaf);
            size_t inserted_index;
            std::vector <T> new_items;
            new_items.reserve(degree + 2);
            for (size_t i = 0; i < items.size(); i++) {
                new_items.push_back(items[i]);
            }

            new_items = run_insert(new_items, data, run_length, inserted_index);
            if (inserted_index == std::numeric_limits<size_t>::max()) {
                return items;
            }


//        if (DEEBUG) cout << "OOOOOO" << inserted_index << " " << new_items.size() << endl;
            // if the new item is added at the end of the node
            if (inserted_index == new_items.size()) {
                if (next != nullptr) {
                    bool res = merge_item_next(new_items);
                    if (!res) {
                        return items;
                    }
                }
            }

            if (inserted_index < 1 ) {
                if (prev != nullptr) {
                    bool res_prev = merge_item_prev(new_items);
                    if (!res_prev) {
                        return items;
                    }
                }
            }

            return new_items;


        }

        // This is the same as the insert function, however this function change the success boolean to true if the insertion
        // was successful and false otherwise.
        vector <T> insert_success(const T &data, size_t run_length, bool &success) {
//        if (DEEBUG) cout << "inserting into node " << data << " " << run_length << endl;
            success = true;
            assert(leaf);
            size_t inserted_index;
            std::vector <T> new_items;
            new_items.reserve(degree + 2);
            for (size_t i = 0; i < items.size(); i++) {
                new_items.push_back(items[i]);
            }


            new_items = run_insert(new_items, data, run_length, inserted_index);
            if (inserted_index == std::numeric_limits<size_t>::max()) {
//                std::cerr << "Here for data " << data << endl;
                success = false;
                return new_items;
            }


//            cout << "OOOOOO" << inserted_index << " " << new_items.size() << endl;
            // if the new item is added at the end of the node
            if ((inserted_index >= new_items.size() || inserted_index == 0) && new_items.size() > 0) {
                if (next != nullptr) {
                    bool res = merge_item_next(new_items);
                    if (!res) {
//                        std::cerr << "MERGE NEXT ERROR for data " << data << endl;
                        success = false;
                        return new_items;
//                        return items;
                    }
                }
            }

            if (inserted_index <= 1 && new_items.size() > 0) {
                if (prev != nullptr) {

                    bool res_prev = merge_item_prev(new_items);
                    if (!res_prev) {
                        success = false;
                        return items;
                    }
                }
            }

//        for (auto &item: new_items) {
//            cout << item << " ";
//        }
//        cout << "Done inserting " << endl;

            return new_items;


        }


        // search an item in the node and return the index of the item
        int search_non_leaf(const T &data) {
            for (int i = 0; i < items.size(); i++) {
                if (data < items[i]) {
                    return i;
                }
                if (i == items.size() - 1) {
                    return i + 1;
                }
            }
            return -1;
        }

        // print function that prints the node items and the size
        void print() {
            std::cout << "Node: ";
            for (auto &item: items) {
                std::cout << item << " ";
            }
            std::cout << "Size: " << items.size() << endl;
        }

        bool is_gap(T run) {
            return run.graph_position.value == 0;
        }

    private:
        bool leaf;
        size_t degree;
//    size_t current_size;
        std::vector <T> items;
        std::vector<bpNode<T> *> children;
        bpNode<T> *parent;
        bpNode<T> *next;
        bpNode<T> *prev;


        T create_gap(size_t gap_start_position) {
            T gap = {gap_start_position, 0};
            return gap;
        }


        /*
         * This function insert a new run into an array of runs without handling any of the merging or splitting and
         * just insert the new run into the array.
         */
        vector <T> run_insert(vector <T> new_items, T new_run, int run_length, size_t &inserted_index) {

            size_t len = new_items.size();
            int index = std::upper_bound(new_items.begin(), new_items.end(), new_run) - new_items.begin();
            int gap_index =
                    std::lower_bound(new_items.begin(), new_items.end(),
                                     create_gap(new_run.start_position + run_length)) -
                    new_items.begin();
//            cerr << "index: " << index << "len " << len << " gap index " << gap_index << endl;

            if (gap_index - index > 0) {
//            if (DEEBUG)
//                cout << "Error: the new run is not inserted correctly. The new run cannot inserted in another run!!"
//                     << endl;
//                cout << "================111111" << endl;
                inserted_index = std::numeric_limits<size_t>::max();
                return new_items;
            }




            if (new_items.size() == 0) {
//            if (DEEBUG) cout << "there is nothing in the node!" << endl;
                new_items.push_back(new_run);
                new_items.push_back(create_gap(new_run.start_position + run_length));
                inserted_index = 0;
                return new_items;
            }



            if (index == 0) {
                // This case is when the new data start position is less than the first item start position

                if (prev != nullptr && !is_gap(prev->get_item(prev->get_size() - 1))) {
                    // this is an error case and should not happen
//                if (DEEBUG)
//                    cerr << "Error: the new run is not inserted correctly. The new run cannot inserted in another run!"
//                         << endl;
//                    cout << "================222222" << endl;
                    inserted_index = std::numeric_limits<size_t>::max();
                    return new_items;

                }
            } else {
                if (!is_gap(new_items[index - 1])) {
//                    cout << "================3333333" << endl;
                    inserted_index = std::numeric_limits<size_t>::max();
                    return new_items;
                }
                if (index == len) {
                    if (next != nullptr && new_run.start_position + run_length > next->get_item(0).start_position) {
//                        cout << "================444444444" << endl;
                        inserted_index = std::numeric_limits<size_t>::max();
                        return new_items;
                    }
                }
            }

//        cout << "index: " << index << "len " << len << new_items[index] << endl;
//        if (index == 1) {
//            cout << "index is 1 " << new_items[index - 1] << " " << new_items[index] << endl;
//        }
            // first case is when the new_run doesn't hit previous run or the next run
            if ((index == 0 || new_run.start_position > new_items[index - 1].start_position) &&
                (index == len || new_run.start_position + run_length < new_items[index].start_position)) {
//                std::cout << "insert run case 1" << std::endl;
                if (index == 0) {
                    new_items.insert(new_items.begin(), new_run);
                    new_items.insert(new_items.begin() + 1, create_gap(new_run.start_position + run_length));
                    inserted_index = 0;
                    return new_items;
                } else if (index == len) {
                    new_items.push_back(new_run);
                    new_items.push_back(create_gap(new_run.start_position + run_length));
                    inserted_index = new_items.size();
                    return new_items;
                } else {
                    new_items.insert(new_items.begin() + index, new_run);
                    new_items.insert(new_items.begin() + index + 1, create_gap(new_run.start_position + run_length));
                    inserted_index = index;
                    return new_items;
                }
            }


                // second case is when we hit the previous run end point and not the next run starting point
            else if (index != 0 && new_run.start_position == new_items[index - 1].start_position &&
                     ((index == len || new_run.start_position + run_length < new_items[index].start_position))) {
//                std::cout << "insert run case 2" << std::endl;
                // the case that index = 1
                if (index == 1) {
                    // if adding after a gap run
                    new_items[index - 1] = new_run;
                    new_items.insert(new_items.begin() + index, create_gap(new_run.start_position + run_length));
                    inserted_index = index - 1;
                    return new_items;
                }
//            else if (index == len){
//                if (DEEBUG) std::cout << "insert run case 2.0" << std::endl;
//                // in this case have to remove the gap run and insert the new run in its place and then add a gap for the new run
//                new_items[index - 1] = new_run;
//
//                new_items.insert(new_items.begin() + index, create_gap(new_run.start_position + run_length));
//                inserted_index = new_items.size();
//                return new_items;
//            }

                    // two cases here, first being the new run graph position not be the same as the previous run graph position
                else if (new_run.graph_position.value != new_items[index - 2].graph_position.value) {
//                if (DEEBUG) std::cout << "insert run case 2.1" << std::endl;
                    // in this case have to remove the gap run and insert the new run in its place and then add a gap for the new run
                    new_items[index - 1] = new_run;
                    inserted_index = index - 1;
                    new_items.insert(new_items.begin() + index, create_gap(new_run.start_position + run_length));
                    if (index == len) inserted_index = new_items.size();
                    return new_items;
                }
                    // the other case is when the new run graph position is the same as the previous non-gap run graph position
                else if (new_run.graph_position.value == new_items[index - 2].graph_position.value) {
//                if (DEEBUG) std::cout << "insert run case 2.2" << std::endl;
                    // in this case we just have to change the starting position of the gap in the index-1
                    new_items[index - 1].start_position += run_length;
                    inserted_index = index - 1;
                    if (index == len) inserted_index = new_items.size();
                    return new_items;
                }

            }


                // third case, when we hit the next run starting point but not the previous one
            else if ((index == 0 || new_run.start_position > new_items[index - 1].start_position) && index != len &&
                     new_run.start_position + run_length == new_items[index].start_position) {
//                std::cout << "insert run case 3" << std::endl;
                // two cases here, first being the new run graph_position not be the same as the next run graph_position
                if (new_run.graph_position.value != new_items[index].graph_position.value) {
//                if (DEEBUG) std::cout << "insert run case 3.1" << std::endl;

                    new_items.insert(new_items.begin() + index, new_run);
                    // in this case we don't need to add a gap after adding the new run
                    inserted_index = index;
                    return new_items;

                }
                    // the other case is when the new run graph_position is the same as the next non-gap run graph_position
                else if (new_run.graph_position.value == new_items[index].graph_position.value) {
//                if (DEEBUG) std::cout << "insert run case 3.2" << std::endl;
                    // in this case we just have to change the starting position of the next non-gap run
                    new_items[index].start_position = new_run.start_position;
                    inserted_index = index;
                    return new_items;
                }
            }




                // forth case, which is when we hit both runs from both sides
            else if (new_run.start_position == new_items[index - 1].start_position &&
                     new_run.start_position + run_length == new_items[index].start_position) {
//                std::cerr << "insert run case 4" << std::endl;

                // if index = 1 then we do not have to check for the index-2 run
                if (index == 1) {
//                if (DEEBUG) std::cout << "insert run case 4.0" << std::endl;
                    // two cases here too! first being the new run graph_position not be the same as the next run graph_position
                    if (new_run.graph_position.value != new_items[index].graph_position.value) {
//                    if (DEEBUG) std::cout << "insert run case 4.0.1" << std::endl;
                        // in this case we don't need to add a gap after adding the new run
                        new_items[index - 1] = new_run;
                        inserted_index = index - 1;
                        return new_items;
                    }
                        // the other case is when the new run graph_position is the same as the next non-gap run graph_position
                    else if (new_run.graph_position.value == new_items[index].graph_position.value) {
//                    if (DEEBUG) std::cout << "insert run case 4.0.2" << std::endl;
                        // in this case we just have to change the starting position of the next non-gap run
                        new_items[index].start_position = new_run.start_position;
                        inserted_index = index - 1;
                        new_items.erase(new_items.begin());
                        return new_items;
                    }
                }

                // two cases here, first being the new run graph_position not be the same as the previous run graph_position
                if ((new_run.graph_position.value != new_items[index - 2].graph_position.value)) {
//                if (DEEBUG) std::cout << "insert run case 4.1" << std::endl;
                    // two cases here too! first being the new run graph_position not be the same as the next run graph_position
                    if (new_run.graph_position.value != new_items[index].graph_position.value) {
//                    if (DEEBUG) std::cout << "insert run case 4.1.1" << std::endl;

                        // in this case we don't merge anything, but we need to remove the index-1 gap run
                        new_items[index - 1] = new_run;
                        inserted_index = index - 1;
                        return new_items;
                    }
                        // the other case is when the new run graph_position is the same as the next non-gap run graph_position
                    else if (new_run.graph_position.value == new_items[index].graph_position.value) {
//                    if (DEEBUG) std::cout << "insert run case 4.1.2" << std::endl;
                        // in this case we just have to change the starting position of the next non-gap run and remove the index-1 gap run
                        new_items[index].start_position = new_run.start_position;
                        new_items.erase(new_items.begin() + index - 1);
                        inserted_index = index - 1;
                        return new_items;

                    }
                }
                    // the other case is when the new run graph_position is the same as the previous non-gap run graph_position
                else if (new_run.graph_position.value == new_items[index - 2].graph_position.value) {
//                if (DEEBUG) std::cout << "insert run case 4.2" << std::endl;
                    // two cases here too! first being the new run graph_position not be the same as the next run graph_position
                    if (new_run.graph_position.value != new_items[index].graph_position.value) {
//                    std::cout << "insert run case 4.2.1" << std::endl;

                        // just remove the index - 1 run
                        new_items.erase(new_items.begin() + index - 1);
                        inserted_index = index - 1;
                        return new_items;
                    }
                        // the other case is when the new run graph_position is the same as the next non-gap run graph_position
                    else if (new_run.graph_position.value == new_items[index].graph_position.value) {
//                    std::cerr << "insert run case 4.2.2" << std::endl;
                        // in this case we have to merge 3 runs together!
                        // in this case we only have to remove the index-1 and index runs
                        new_items.erase(new_items.begin() + index - 1);
                        new_items.erase(new_items.begin() + index - 1);
                        inserted_index = index - 1;
                        return new_items;
                    }
                }
            } else {
//            if (DEEBUG) {
//                std::cerr << "Error: the new run is not inserted correctly" << std::endl;
//                // printing some DEEBUG information
//                std::cerr << "new run: " << new_run << std::endl;
//                std::cerr << "index: " << index << std::endl;
//                std::cerr << "new items: ";
//                for (auto &item: new_items) {
//                    std::cerr << item << " ";
//                }
//                std::cerr << std::endl;
//
//            }
//                cout << "================555555555" << endl;
                inserted_index = std::numeric_limits<size_t>::max();
                return new_items;
            }
//        return arr;
            return new_items;
        }


        bool merge_item_next(vector <T> &new_items) {
//            cout << "merging with the next node" << endl;
            assert(next != nullptr);
            vector <T> next_node_items = next->get_items();

            if (!is_gap(new_items[new_items.size() - 1])) {
//                cerr << "MERGE: GAP";
                return true;
            }
            if (next_node_items[0].start_position < new_items[new_items.size() - 1].start_position) {
//                cerr << "MERGE: CASE2" << endl;
//                next->print();
//            new_items.pop_back();

                return false;
            } else if (next_node_items[0].start_position == new_items[new_items.size() - 1].start_position) {
                if (next_node_items[0].graph_position.value == new_items[new_items.size() - 2].graph_position.value) {
                    new_items.pop_back();
                    next->remove_first_item();
//                cout << "++++++++ " << new_items[new_items.size() - 1] << "----" << next_node_items[0] << endl;
                    return true;
                } else {
                    new_items.pop_back();
                    return true;
                }
            }
            return true;
        }


        bool merge_item_prev(vector <T> &new_items) {
            assert(prev != nullptr);
            assert(!new_items.empty());
//            cout << "merging with the prev node" << endl;
            vector <T> prev_node_items = prev->get_items();

//            std::cerr << new_items[0] << " " <<  prev_node_items[prev_node_items.size() - 1] << endl;
//        if (!is_gap(prev_node_items[prev_node_items.size() - 1])){
//
//            return false;
//        }
            if (prev_node_items[prev_node_items.size() - 1].graph_position.value == new_items[0].graph_position.value) {
                new_items.erase(new_items.begin());
                return true;
            }
            return true;

        }


    };


    template<typename T>
    class BplusTree {
    private:
        bpNode<T> *root;
        std::size_t degree;


        // this function handles the insertion of a new item/child into a parent node
        void parent_insert(bpNode<T> *node, T data, bpNode<T> *child) {
//        if (DEEBUG) cout << "inserting into parent node" << endl;

            // parent overflow case
            if (node->get_size() == degree) {
//            if (DEEBUG) cout << "overflow on the parent node" << endl;
                // if the parent is full first check for the left sibling

                bpNode<T> *new_node = new bpNode<T>(degree, false);
//            T parent_item = node->get_item(degree / 2);
                vector <T> all_items;
                vector < bpNode<T> * > all_children;
                all_items.reserve(degree + 1);
                all_children.reserve(degree + 2);
                // TODO: complete this. In the current version I am adding the new children at the end of the new node!!!
                vector <T> new_node_items;
                vector < bpNode<T> * > new_node_children;

                for (auto &items: node->get_items()) {
                    all_items.push_back(items);
                }

                for (auto &child_: node->get_children()) {
                    all_children.push_back(child_);
                }

                // finding the adding index
                int add_index = node->search_non_leaf(data);
                all_items.insert(all_items.begin() + add_index, data);
                all_children.insert(all_children.begin() + add_index + 1, child);


                vector <T> old_node_items;
                vector < bpNode<T> * > old_node_children;
                for (int i = 0; i < degree / 2 + 1; i++) {
                    old_node_items.push_back(all_items[i]);
                }

                for (int i = 0; i < degree / 2 + 2; i++) {
                    all_children[i]->set_parent(node);
                    old_node_children.push_back(all_children[i]);
                }


                for (int i = degree / 2 + 2; i < degree + 1; i++) {
                    new_node_items.push_back(all_items[i]);
                }


                for (int i = degree / 2 + 2; i < degree + 2; i++) {
                    all_children[i]->set_parent(new_node);
                    new_node_children.push_back(all_children[i]);

                }

                // add the new child to the new node
//            new_node_items.push_back(data);
//            new_node_children.push_back(child);

//            child->set_parent(new_node);


                new_node->replace_items_leaf(new_node_items);
                new_node->replace_children(new_node_children);
//            for (auto &child1: new_node_children) {
//                new_node->add_child(child1);
//            }

                node->replace_items_leaf(old_node_items);
                node->replace_children(old_node_children);
//            for (auto &child1: old_node_children) {
//                node->add_child(child1);
//            }

                T parent_item = all_items[degree / 2 + 1];


//            T parent_item = new_node->get_item(0);

                if (node->is_root()) {
                    bpNode<T> *new_root = new bpNode<T>(degree, false);

                    new_root->add_child(node);
                    new_root->add_child(new_node);
                    new_root->add_item(parent_item);

                    node->set_parent(new_root);
                    new_node->set_parent(new_root);

                    root = new_root;
                } else {
                    parent_insert(node->get_parent(), parent_item, new_node);
                }
            } else {
                // if the parent is not full then just add the item and the child

                // finding the index to add the item and the children
                int add_index = node->search_non_leaf(data);

                child->set_parent(node);
                node->add_item(data, add_index);
                node->add_child(child, add_index + 1);

            }
        }


        // this function handles the removing the item from a non-leaf node, for example when merging two nodes we need to
        // remove the item from the parent node and it might cause underflow in the parent node, ... this function handles
        // all of these cases. The input node is the node that we want to remove from the tree.

        void remove_from_parent(bpNode<T> *node) {
//            cerr << "removing from parent node" << endl;

            bpNode<T> *parent = node->get_parent();
            assert(parent != nullptr);

            int index = node->find_index_in_parent();

//            std::cerr << "Index and The parent node is " << index << endl;
//            parent->print();

//            std::cerr << "The item we are deleting " << parent->get_item(index - 1) << endl;
//            std::cerr << "parent child size " << parent->get_children().size() << endl;
//            std::cerr << "The child that we are removing is " << endl;
//            parent->get_child(index)->print();
            // Remove the item and the corresponding child from the parent

//            std::cerr << "index " << index << endl;

            if (index == 0) {
                parent->remove_item(0);
            } else {
                parent->remove_item(index - 1);
            }
            parent->remove_child(index);



//            std::cerr << "PARENT AFTER REMOVING size " << parent->get_children().size() << endl;
//            parent->print();


//            std::cerr << "PARENT CHILDREN AFTER REMOVING " << endl;
//            for (auto &chil: parent->get_children()) {
//                chil.get_items()[0];
//            }

//        if (DEEBUG) {
//            if (parent->get_size() > 0) {
//                // print all items in the parent node
//                cout << "parent items: ";
//                for (auto &item: parent->get_items()) {
//                    cout << item << " ";
//                }
//                cout << endl;
//
//            }
//
//        }

            // If the parent is the root and it becomes empty, make the current node the new root
            if (parent->is_root() && parent->get_size() == 0) {

//                std::cerr << "FUCKED UP 8" << endl;
                bpNode<T> *new_root = parent->get_child(0);

                root = new_root;
                root->set_parent(nullptr);
                root->set_next(nullptr);
                root->set_prev(nullptr);

                delete parent;
                return;
            }

            // If the parent is underflowing, handle the underflow
            if (parent->is_underflowing()) {
//                std::cerr << "recursively underflowing" << endl;
                handle_parent_underflow(parent);
            }
        }

        // This function handles the underflow of a non-leaf (parent) node.
        void handle_parent_underflow(bpNode<T> *node) {
//            cout << "handling parent underflow" << endl;

            assert(!node->is_leaf());

            // If the node is the root and it is empty, make its only child the new root
            // This will probably never happen!
            if (node->is_root() && node->get_size() == 0) {

//                std::cerr << "FUCKED UP 1" << endl;

                root = node->get_child(0);
                root->set_parent(nullptr);
                root->set_next(nullptr);
                root->set_prev(nullptr);

                delete node;
                return;
            }

            if (node->is_root()) {
//                std::cerr << "FUCKED UP 2" << endl;
                return;
            }

            assert(node->get_parent() != nullptr);

            bpNode<T> *parent = node->get_parent();
            int index = node->find_index_in_parent();

            // Borrow from the left sibling
            if (index > 0 && parent->get_child(index - 1)->get_size() > degree / 2) {
//                std::cerr << "FUCKED UP 3" << endl;
//            if (DEEBUG) cout << "borrowing from the left sibling (parent underflow)" << endl;

                bpNode<T> *left_sibling = parent->get_child(index - 1);


                T borrowed_item = left_sibling->get_item(left_sibling->get_size() - 1);
                bpNode<T> *borrowed_child = left_sibling->get_child(left_sibling->get_size());

                left_sibling->remove_last_item();
                left_sibling->remove_last_child();

                node->add_item(parent->get_item(index - 1), 0);

                borrowed_child->set_parent(node);
                node->add_child(borrowed_child, 0);


//            if (DEEBUG) {
//                // print all items in the node
//                cout << "node items in the node that just got its item/child from its leftsibling: ";
//                for (auto &item: node->get_items()) {
//                    cout << item << " ";
//                }
//            }

                parent->change_item(borrowed_item, index - 1);


                // Borrow from the right sibling
            } else if (index < parent->get_size() && parent->get_child(index + 1)->get_size() > degree / 2) {
//                std::cerr << "FUCKED UP 4" << endl;
//            if (DEEBUG) cout << "borrowing from the right sibling (parent underflow)" << endl;

                bpNode<T> *right_sibling = parent->get_child(index + 1);
                T borrowed_item = right_sibling->get_item(0);
                bpNode<T> *borrowed_child = right_sibling->get_child(0);

                right_sibling->remove_first_item();
                right_sibling->remove_child(0);

                node->add_item(parent->get_item(index));
                node->add_child(borrowed_child);
                parent->change_item(borrowed_item, index);

                borrowed_child->set_parent(node);

                // Merge with the left sibling
            } else if (index > 0) {
//                std::cerr << "FUCKED UP 5" << endl;
//            if (DEEBUG) cout << "merging with the left sibling (parent underflow)" << endl;

                bpNode<T> *left_sibling = parent->get_child(index - 1);
                left_sibling->add_item(parent->get_item(index - 1));
                left_sibling->add_child(node->get_child(0));

                for (int i = 0; i < node->get_size(); i++) {
                    left_sibling->add_item(node->get_item(i));
                    left_sibling->add_child(node->get_child(i + 1));
                }

//            if (DEEBUG) {
//                cout << "left sibling items" << endl;
//                for (auto &item: left_sibling->get_items()) {
//                    cout << item << " ";
//                }
//            }

                // handle the parent of the newly added children
                for (int i = 0; i < left_sibling->get_size() + 1; i++) {
                    left_sibling->get_child(i)->set_parent(left_sibling);
                }

                remove_from_parent(node);
                delete node;

                // Merge with the right sibling
            } else if (index < parent->get_size()) {
//                std::cerr << "FUCKED UP 6" << endl;
//            if (DEEBUG) cout << "merging with the right sibling (parent underflow)" << endl;

                bpNode<T> *right_sibling = parent->get_child(index + 1);


                node->add_item(parent->get_item(index));
                node->add_child(right_sibling->get_child(0));


                for (int i = 0; i < right_sibling->get_size(); i++) {
                    node->add_item(right_sibling->get_item(i));
                    node->add_child(right_sibling->get_child(i + 1));
                }

                // change the parents of the newly added children
                for (int i = 0; i < node->get_size() + 1; i++) {
                    node->get_child(i)->set_parent(node);
                }

                remove_from_parent(right_sibling);
                delete right_sibling;
            }
        }


        // this function handles the underflow of a leaf node. There are multiple that we can handle the underflow, borrowing
        // from the left sibling, borrowing from the right sibling, merging with the left sibling and merging with the right sibling
        void leaf_underflow(bpNode<T> *node) {

            assert(node->is_leaf());

//        if (DEEBUG) cout << "underflow in the leaf node" << endl;
//        if (DEEBUG)
//            cout << "node size: " << node->get_size() << " node first item " << node->get_item(0) << node->is_root()
//                 << endl;

            if (node->is_root()) {
                if (node->get_size() == 0) {
                    // the tree is empty
                    delete node;

                }
                // there is nothing to do, the tree just has one node
                return;
            }


            // if reached here then the node is not the root and has a parent
            assert(node->get_parent() != nullptr);

            bpNode<T> *parent = node->get_parent();
            int index = node->find_index_in_parent();

            if (index == 0) assert(parent->get_children().size() > 0);


//        if (DEEBUG) cout << index << endl;

            // the left siblings exists and has enough items to borrow
            if (index > 0 && parent->get_child(index - 1)->get_size() > (degree + 1) / 2 - 1) {
                // borrow from the left sibling
//            if (DEEBUG) cout << "borrowing from the left sibling" << endl;

                bpNode<T> *left_sibling = parent->get_child(index - 1);


                assert(left_sibling->get_next() == node);
                assert(node->get_prev() == left_sibling);

                T borrowed_item = left_sibling->get_item(left_sibling->get_size() - 1);
                node->add_item(borrowed_item, 0);
                left_sibling->remove_last_item();
                parent->change_item(borrowed_item, index - 1);


                // the right sibling exists and has enough items to borrow
            } else if (index < parent->get_size() && parent->get_child(index + 1)->get_size() > (degree + 1) / 2 - 1) {
                // borrow from the right sibling
//            if (DEEBUG) cout << "borrowing from the right sibling" << endl;

                bpNode<T> *right_sibling = parent->get_child(index + 1);

                assert(right_sibling->get_prev() == node);
                assert(node->get_next() == right_sibling);

                T borrowed_item = right_sibling->get_item(0);

                node->add_item(borrowed_item);
                right_sibling->remove_first_item();

                parent->change_item(right_sibling->get_item(0), index);


                // left sibling exists and we can merge with it
            } else if (index > 0) {
                // merge with the left sibling
//            if (DEEBUG) cout << "merging with the left sibling" << endl;
                bpNode<T> *left_sibling = parent->get_child(index - 1);

//                std::cerr << "LEFT SIBLING " << endl;
//                left_sibling->print();


                assert(left_sibling->get_next() == node);
                assert(node->get_prev() == left_sibling);


                for (int i = 0; i < node->get_size(); i++) {
                    left_sibling->add_item(node->get_item(i));
                }

                left_sibling->set_next(node->get_next());
                if (node->get_next() != nullptr) {
                    node->get_next()->set_prev(left_sibling);
                }

//            if (DEEBUG) cout << "removing the node with the first item: " << node->get_item(0) << endl;
                remove_from_parent(node);
                delete node;

//                std::cerr << "==============================Node Deleted=======================" << std::endl;

                // TODO: call the remove from a non-leaf node function
                // TODO: remove the actual node itself

                if (left_sibling->get_next() != nullptr) assert(left_sibling->get_next()->get_prev() == left_sibling);

                // right sibling exists and we can merge with it
            } else if (index <= parent->get_size()) {
                // merge with the right sibling
//            if (DEEBUG) cout << "merging with the right sibling" << endl;
                bpNode<T> *right_sibling = parent->get_child(index + 1);

//                std::cerr << "RIGHT SIBLING " << endl;
//                right_sibling->print();

                assert(right_sibling->get_prev() == node);
                assert(node->get_next() == right_sibling);

                for (int i = 0; i < right_sibling->get_size(); i++) {
                    node->add_item(right_sibling->get_item(i));
                }

                node->set_next(right_sibling->get_next());
                if (node->get_next() != nullptr) {
                    node->get_next()->set_prev(node);
                }

                remove_from_parent(right_sibling);
                delete right_sibling;

                // TODO: call the remove from a non-leaf node function
                // TODO: remove the actual node itself

//                std::cerr << "++++++++++++++++++++++++++ RIGHT SIBLING DELETED +++++++++++++++++++++++" << std::endl;
                assert(node->get_next()->get_prev() == node);
            }
        }


    public:
        BplusTree(std::size_t _degree) {// Constructor
            this->root = nullptr;
            this->degree = _degree;
        }

        // Destructor
        ~BplusTree() {
            delete root;
        }

        // get root
        bpNode<T> *get_root() {
            return root;
        }

        bpNode<T> *leaf_search(bpNode<T> *node, T data) {
            assert(node != nullptr);
            bpNode<T> *current = node;

            while (!current->is_leaf()) {
                current = current->get_child(current->search_non_leaf(data));
            }
            return current;

        }

        size_t *search_index(T data) {
            bpNode<T> *current = leaf_search(root, data);

            // find the index of the maximum data before the data
            return current->search(data);

        }

        T search(size_t position) {
            T data = {position, 0};
            bpNode<T> *current = leaf_search(root, data);
            size_t index = current->search(data);
            if (index == 0) {
                if (current->get_prev() != nullptr) {
                    return current->get_prev()->get_item(current->get_prev()->get_size() - 1);
                } else {
                    return {0, 0};
                }

            }
            return current->get_item(index - 1);

        }


        void insert(const T &data, size_t run_length) {
            bool temp = insert_success(data, run_length);
        }


        bool insert_success(const T &data, size_t run_length) {
//            std::cerr << "inserting: " << data << " " << run_length << endl;
            bool success = true;
            if (root == nullptr) {
                root = new bpNode<T>(degree, true);
                vector <T> temp = root->insert(data, run_length);
                root->replace_items_leaf(temp);

            } else {
//                Run temp_data = {9211, 0};
//                Run tt_data = data;
//                if (data.start_position == 9211) {
////                    std::cerr << "000000000000000000" << endl;
//                    std::cerr << "inserting: " << data << " " << run_length << endl;
//                    bpNode<T> *test = leaf_search(root, temp_data);
//                    std::cerr << "000000000000000000TEST" << endl;
//                    test->print();
//                    if (test->get_next() != nullptr) {
//                        std::cerr << "000000000000000000TEST NEXT" << endl;
//                        test->get_next()->print();
//                    }
//                }
//                bpNode<T> *test = leaf_search(root, temp_data);
//                std::cerr << "000000000000000000TEST" << endl;
//                test->print();
//                if (test->get_next() != nullptr) {
//                    std::cerr << "000000000000000000TEST NEXT" << endl;
//                    test->get_next()->print();
//                }
//                std::cerr << "000000000000000000TEST END" << endl;

                bpNode<T> *leaf = leaf_search(root, data);
//                leaf->print();
                vector <T> inserted_items = leaf->insert_success(data, run_length, success);
                if (!success) {
//                    std::cerr << "NOT SUCCESS" << endl;
                    return success;
                }


                // in this case we want to completely remove this node from the tree
                if (inserted_items.size() == 0){

//                    std::cerr << "There is a completely removed node" << endl;

//                    auto *next1 = leaf->get_next();
//                    if (next1 != nullptr){
//                        leaf_underflow(next1);
//                    }
//                    leaf->replace_items_leaf(inserted_items);
//                    leaf_underflow(leaf);





                    if (leaf->get_next() != nullptr) {
                        leaf->get_next()->set_prev(leaf->get_prev());
                    }

                    assert(leaf->get_next()->get_prev() == leaf);


                    if (leaf->get_prev() != nullptr) {
                        leaf->get_prev()->set_next(leaf->get_next());
                    }
                    assert(leaf->get_prev()->get_next() == leaf);
                    remove_from_parent(leaf);

                    delete leaf;

                    return success;

                }

//                if (leaf->get_next() != nullptr && leaf->get_next()->is_underflowing()){
//                    std::cerr << "Next NODE IS UNDERFLOWING AFTER INSERT SUCCESS" << std::endl;
//                }
//                if (leaf->get_prev() != nullptr && leaf->get_prev()->is_underflowing()){
//                    std::cerr << "PREVVV NODE IS UNDERFLOWING AFTER INSERT SUCCESS" << std::endl;
//                }
                // print the items in inserted_items
//            if (DEEBUG) for (auto &item: inserted_items) { cout << "item: " << item << endl; }

                if (inserted_items.size() > 0 && inserted_items.size() <= degree) {
                    leaf->replace_items_leaf(inserted_items);
                } else if (inserted_items.size() > 0){
//                    if (leaf->get_next() != nullptr && leaf->get_next()->is_underflowing()){
//                        std::cerr << "((((((((((FURST)))))))))))" << endl;
//                        int index = leaf->get_next()->find_index_in_parent();
//                        T borrowed_item = inserted_items.back();
//                        leaf->get_next()->add_item(borrowed_item, 0);
//                        inserted_items.pop_back();
//                        leaf->get_parent()->change_item(borrowed_item, index - 1);
//                    }
//
//                    if (leaf->get_prev() != nullptr && leaf->get_prev()->is_underflowing()){
//                        std::cerr << "((((((((((Second)))))))))))" << endl;
//
//                        int index = leaf->find_index_in_parent();
//                        T borrowed_item = inserted_items[0];
//                        inserted_items.erase(inserted_items.begin());
//                        leaf->get_prev()->add_item(borrowed_item);
//
//                        leaf->get_parent()->change_item(inserted_items[0], index - 1);
//                    }
//
//                    if (inserted_items.size() <= degree){
//                        std::cerr << "((((((((((Third)))))))))))" << endl;
//                        leaf->replace_items_leaf(inserted_items);
//                    } else {

//                    std::cerr << "EASY CASE" << endl;



                    // split the node
//                if (DEEBUG) cout << "splitting the node" << endl;
                    bpNode<T> *new_node = new bpNode<T>(degree, true);
                    // the new_node items are the second half of the leaf items
                    vector <T> new_node_items;


                    int all_items_len = inserted_items.size();
                    // store the second half of the items in the new node
                    for (int i = all_items_len / 2; i < all_items_len; i++) {
                        new_node_items.push_back(inserted_items[i]);
//                    cout << "---------: " << inserted_items[i] << endl;
                    }

                    // remove the second half of the items from the leaf node
                    for (int i = all_items_len / 2; i < all_items_len; i++) {
                        inserted_items.pop_back();
                    }

                    // print all items in inserted_items
//                if (DEEBUG) for (auto &item: inserted_items) { cout << "ITEM: " << item << endl; }

                    // change the items of the nodes
                    leaf->replace_items_leaf(inserted_items);
                    new_node->replace_items_leaf(new_node_items);

                    // set the next and prev pointers
                    new_node->set_next(leaf->get_next());
                    new_node->set_prev(leaf);
                    leaf->set_next(new_node);
                    if (new_node->get_next() != nullptr) {
                        new_node->get_next()->set_prev(new_node);
                    }
                    assert(leaf->get_next()->get_prev() == leaf);


                    T parent_item = new_node->get_item(0);


                    // handle the parents
                    if (leaf->is_root()) {
                        bpNode<T> *new_root = new bpNode<T>(degree, false);

                        // handle the children
                        new_root->add_child(leaf);
                        new_root->add_child(new_node);


                        // handle the items
                        new_root->add_item(parent_item);

                        // handle the parents
                        leaf->set_parent(new_root);
                        new_node->set_parent(new_root);


                        root = new_root;

                    } else {
                        // a parent exists
                        parent_insert(leaf->get_parent(), parent_item, new_node);
                    }
//                    }
                }

//                std::cerr << "HERE" << endl;

                if (leaf->get_next() != nullptr) assert(leaf->get_next()->get_prev() == leaf);

//                if (data.start_position == 9211) {
//                    std::cerr << "XXXXXXPREVXXXXXXPREVXXXXXXPREVXXXXXXPREVXXXXXXPREVXXXXXXPREVXXXXXXPREV" << endl;
//                    leaf->get_prev()->print();
//                    std::cerr << "LEAF" << endl;
//                    leaf->print();
//                    std::cerr << "NEXT" << endl;
//                    leaf->get_next()->print();
//
//                }

                // handle the underflow of the next node
                if (leaf->get_next() != nullptr && leaf->get_next()->is_underflowing()) {

//                    cerr << "underflow in the next node " << leaf->get_size() << " " << leaf->get_next()->get_size() << endl;
                    leaf_underflow(leaf->get_next());
//                    cerr << leaf->get_size() << " " << leaf->get_next()->get_size() << endl;
                }


                // handle the underflow
                if (leaf->is_underflowing() && leaf != nullptr && leaf->get_size() != 0) {
//                    cerr << "working node " << leaf->get_size() << endl;
                    leaf_underflow(leaf);
//                    cerr << "WORKING NODE AFTER " << leaf->get_size() << endl;
                }

//                if (leaf->get_prev() != nullptr && leaf->get_prev()->is_underflowing()) {
//                    leaf_underflow(leaf->get_prev());
//                }


            }
            return success;

        }


        void print(bpNode<T> *node) {
            bpNode<T> *cursor = node;
            if (cursor != NULL) {
                for (int i = 0; i < cursor->get_size(); ++i) {
                    std::cout << cursor->get_item(i) << " ";
                }
                if (cursor->get_next() != nullptr)
                    std::cout << " node next first item " << cursor->get_next()->get_item(0) << endl;
                std::cout << cursor->is_root() << "\n";

                if (!cursor->is_leaf()) {
                    for (int i = 0; i < cursor->get_size() + 1; ++i) {
                        print(cursor->get_child(i));
                    }
                }
            }
        }

        void print_tree(bpNode<Run> *node, int level = 0) {
            if (node == nullptr) return;

            // Indent to represent tree level
            std::string indent(level * 4, ' ');

            // Label the node
            std::cout << indent << (node->is_leaf() ? "Leaf: " : "Internal: ");
            for (int i = 0; i < node->get_size(); ++i) {
                std::cout << "[" << node->get_item(i).start_position << ":"
                          << node->get_item(i).graph_position.value << "] ";
            }
            std::cout << (node->is_root() ? "(root)" : "") << "\n";

            if (!node->is_leaf()) {
                for (int i = 0; i <= node->get_size(); ++i) {
                    print_tree(node->get_child(i), level + 1);
                }
            }
        }

        void print_whole_tree() {
            print_tree(root);
        }

        void print_whole() {
            print(root);
        }

        size_t get_bpt_size() {
            return get_size(this->root);
        }

        size_t get_size(bpNode<T> *cursor) {
            size_t size = 0;
            if (cursor != NULL) {
                if (cursor->is_leaf()) {
                    size += cursor->get_size();
                } else {
                    for (int i = 0; i < cursor->get_size() + 1; ++i) {
                        size += get_size(cursor->get_child(i));
                    }
                }
            }
            return size;
        }

        size_t get_bpt_run_count() {
            return get_run_count(this->root);
        }

        size_t get_run_count(bpNode<T> *cursor) {
            size_t count = 0;
            if (cursor != NULL) {
                if (cursor->is_leaf()) {
                    for (int i = 0; i < cursor->get_size(); ++i) {
                        if (cursor->get_item(i).graph_position.value != 0) {
                            count++;
                        }
                    }
                } else {
                    for (int i = 0; i < cursor->get_size() + 1; ++i) {
                        count += get_run_count(cursor->get_child(i));
                    }
                }
            }
            return count;
        }


        class Iterator {
            bpNode<T> *node;
            int index;

        public:
            Iterator(bpNode<T> *node, int index) : node(node), index(index) {}

            T operator*() const {
                return node->get_item(index);
            }

            Iterator &operator++() {
                if (node != nullptr) {
                    index++;
                    if (index >= node->get_size()) {
                        node = node->get_next();
                        index = 0;
                    }
                }
                return *this;
            }

            Iterator &operator--() {
                if (node != nullptr) {
                    index--;
                    if (index < 0) {
                        node = node->get_prev();
                        index = node->get_size() - 1;
                    }
                }
                return *this;
            }

            bool operator!=(const Iterator &other) const {
                return node != other.node || index != other.index;
            }

            bool operator==(const Iterator &other) const {
                return node == other.node && index == other.index;
            }
        };

// Ensure the root is not null and find the leftmost leaf node
        Iterator begin() {
            bpNode<T> *current = root;
            if (!current) return Iterator(nullptr, 0);
            while (current && !current->is_leaf()) {
                current = current->get_child(0);
            }
            return Iterator(current, 0);
        }

// The end iterator represents an iterator past the last element
        Iterator end() {
            return Iterator(nullptr, 0);
        }


    };


}


#endif // PANGENOME_INDEX_BPLUS_TREE_HPP
