//
// Created by seeskand on 5/22/24.
//
/*
   A B+ tree that stores the run-length encoding of a string.

   Each item is a Run. The key is the position where the run starts and the value is its graph position. Inserting
   a run of some length adds two items: the run itself and a "gap" at the position where the run ends. For example,
   inserting Run(1, G) with length 5 adds Run(1, G) and the gap Run(6, no position). A gap has no graph position.
   This way the lengths of the runs are never stored. Runs next to each other with the same graph position are
   merged into one, and the tree splits and merges its nodes as items come and go.

   The implementation is in src/bplus_tree.cpp. The tree is only used with Run, so it is instantiated there.
*/
#ifndef PANGENOME_INDEX_BPLUS_TREE_HPP
#define PANGENOME_INDEX_BPLUS_TREE_HPP

#include <cstddef>
#include <ostream>
#include <vector>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/index.h>

namespace panindexer {

    // An item of the B+ tree: where a run starts and the graph position it has.
    struct Run {
        size_t start_position;
        gbwtgraph::Position graph_position;

        // A gap that starts at the given position.
        static Run gap(size_t start_position) {
            return {start_position, gbwtgraph::Position::no_pos()};
        }

        bool is_gap() const {
            return graph_position.value == 0;
        }

        // Runs are ordered by their start position only.
        bool operator<(const Run &other) const {
            return start_position < other.start_position;
        }

        bool operator==(const Run &other) const {
            return start_position == other.start_position && graph_position.value == other.graph_position.value;
        }

        friend std::ostream &operator<<(std::ostream &os, const Run &run) {
            os << "start_position: " << run.start_position
               << ", graph_position: " << run.graph_position.value;
            return os;
        }
    };

    // A node of the B+ tree. A leaf holds the items and is linked to the leaves before and after it. An internal
    // node holds keys and has one more child than it has keys. A node does not own its children, the tree frees
    // all of the nodes.
    template<typename T>
    class bpNode {
    public:
        bpNode(size_t degree, bool leaf = false) : leaf(leaf), degree(degree) {
            items.reserve(degree);
            if (!leaf) {
                children.reserve(degree + 1);
            }
        }

        void set_next(bpNode<T> *next_node) {
            next = next_node;
        }

        void set_prev(bpNode<T> *prev_node) {
            prev = prev_node;
        }

        void set_parent(bpNode<T> *parent_node) {
            parent = parent_node;
        }

        bpNode<T> *get_next() const {
            return next;
        }

        bpNode<T> *get_prev() const {
            return prev;
        }

        bpNode<T> *get_parent() const {
            return parent;
        }

        bool is_leaf() const {
            return leaf;
        }

        bool is_root() const {
            return parent == nullptr;
        }

        // The number of items in the node.
        size_t get_size() const {
            return items.size();
        }

        const T &get_item(size_t i) const {
            return items[i];
        }

        bpNode<T> *get_child(size_t i) const {
            return children[i];
        }

        const std::vector<T> &get_items() const {
            return items;
        }

        const std::vector<bpNode<T> *> &get_children() const {
            return children;
        }

        void add_item(T item) {
            items.push_back(item);
        }

        // add the item at the ith position
        void add_item(T item, size_t i) {
            items.insert(items.begin() + i, item);
        }

        void add_child(bpNode<T> *child) {
            children.push_back(child);
        }

        // add the child at the ith position
        void add_child(bpNode<T> *child, size_t i) {
            children.insert(children.begin() + i, child);
        }

        void change_item(T new_item, size_t i) {
            items[i] = new_item;
        }

        void remove_item(size_t i) {
            items.erase(items.begin() + i);
        }

        void remove_first_item() {
            items.erase(items.begin());
        }

        void remove_last_item() {
            items.pop_back();
        }

        void remove_child(size_t i) {
            children.erase(children.begin() + i);
        }

        void remove_last_child() {
            children.pop_back();
        }

        void replace_items(const std::vector<T> &new_items) {
            items.assign(new_items.begin(), new_items.end());
        }

        void replace_children(const std::vector<bpNode<T> *> &new_children) {
            children.assign(new_children.begin(), new_children.end());
        }

        // The index of the first item that is greater than data, or the number of items if there is none.
        // In an internal node this is the index of the child to go down to.
        size_t search(const T &data) const;

        // find the index of this node in its parent
        size_t find_index_in_parent() const;

        bool is_underflowing() const;

        /*
           Insert a new run into this leaf. The items that the leaf should have after the insertion are stored in
           new_items, and it is on the tree to put them in the leaf, or to split or remove the leaf if they are
           too many or none. If the run reaches the previous or the next leaf it is merged with the items there.
           Returns false if the run overlaps another run. In that case nothing is changed.
        */
        bool insert(const T &data, size_t run_length, std::vector<T> &new_items);

    private:
        bool leaf;
        size_t degree;
        std::vector<T> items;
        std::vector<bpNode<T> *> children;
        bpNode<T> *parent = nullptr;
        bpNode<T> *next = nullptr;
        bpNode<T> *prev = nullptr;

        /*
           Insert a new run into the items of this leaf and merge it with the runs it touches in there. The items
           of the other leaves are not changed. Returns false if the new run overlaps another run. inserted_index
           is set to the index the run was inserted at, or to the number of items if it was added at the end.
        */
        bool run_insert(std::vector<T> &new_items, const T &new_run, size_t run_length, size_t &inserted_index);

        // Called when new_items may end right where the next leaf starts. If new_items ends with a gap that
        // starts at the first item of the next leaf, the gap is not needed anymore. If that first item also has
        // the same graph position as our last run, it is merged into our run and removed from the next leaf.
        // Returns false if the gap goes past the start of the next leaf.
        bool merge_item_next(std::vector<T> &new_items);

        // If the previous leaf ends with the same graph position that new_items starts with, the first item of
        // new_items just continues that one (a run after the same run, or a gap after a gap) and is removed.
        void merge_item_prev(std::vector<T> &new_items);
    };


    // The tree itself. All of the items are in the leaves, and the internal nodes only have keys that guide the
    // search. The tree finds the leaf for a run, and keeps itself balanced by splitting the nodes that get full
    // and by borrowing or merging when they get too small.
    template<typename T>
    class BplusTree {
    public:
        // The degree is the maximum number of items in a node. It has to be at least 5, so that a leaf always
        // has two items or more.
        explicit BplusTree(size_t degree);

        ~BplusTree();

        // The tree owns its nodes through raw pointers, so it cannot be copied.
        BplusTree(const BplusTree &) = delete;
        BplusTree &operator=(const BplusTree &) = delete;

        // Insert a run with the given length. Returns false if the run overlaps a run that is already in the
        // tree, and the tree is not changed in that case.
        bool insert(const T &data, size_t run_length);

        // Returns the last item that starts at or before the position. This is the run that covers the position,
        // or a gap if no run covers it.
        T search(size_t position) const;

        // The number of items in the tree, with the gaps.
        size_t get_bpt_size() const;

        // The number of runs in the tree, without the gaps.
        size_t get_bpt_run_count() const;

        // Print the tree, one node per line. Only for debugging.
        void print(std::ostream &out) const;

        // Goes over the items of the leaves in order.
        class Iterator {
            bpNode<T> *node;
            size_t index;

        public:
            Iterator(bpNode<T> *node, size_t index) : node(node), index(index) {}

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
                    if (index > 0) {
                        index--;
                    } else {
                        node = node->get_prev();
                        index = (node != nullptr) ? node->get_size() - 1 : 0;
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

        Iterator begin() const {
            return Iterator(first_leaf(), 0);
        }

        // the end iterator is the one past the last item
        Iterator end() const {
            return Iterator(nullptr, 0);
        }

    private:
        bpNode<T> *root;
        size_t degree;

        void destroy(bpNode<T> *node);

        // The leaf that the data belongs to.
        bpNode<T> *leaf_search(const T &data) const;

        bpNode<T> *first_leaf() const;

        // Make a new root over two nodes, where data is the key between them.
        void new_root(bpNode<T> *left, const T &data, bpNode<T> *right);

        // Split a leaf that got too many items. The leaf keeps the first half of the items and a new leaf right
        // after it gets the second half.
        void split_leaf(bpNode<T> *leaf, const std::vector<T> &items);

        // Add a new key and the child that comes after it to an internal node. The node is split if it is full.
        void parent_insert(bpNode<T> *node, const T &data, bpNode<T> *child);

        // Remove a node from its parent, with the key that goes with it. The node itself is not deleted here.
        // This can leave the parent with too few children, which is handled here too.
        void remove_from_parent(bpNode<T> *node);

        // Handles the underflow of an internal node. It borrows a child from the left or the right sibling if
        // they have enough, and merges with one of them otherwise.
        void handle_parent_underflow(bpNode<T> *node);

        // Handles the underflow of a leaf. It borrows an item from the left or the right sibling if they have
        // enough, and merges with one of them otherwise.
        void leaf_underflow(bpNode<T> *node);

        void print(std::ostream &out, const bpNode<T> *node, size_t level) const;
    };

    // Defined in src/bplus_tree.cpp.
    extern template class bpNode<Run>;
    extern template class BplusTree<Run>;

}

#endif // PANGENOME_INDEX_BPLUS_TREE_HPP
