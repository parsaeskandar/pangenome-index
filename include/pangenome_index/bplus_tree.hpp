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
*/
#ifndef PANGENOME_INDEX_BPLUS_TREE_HPP
#define PANGENOME_INDEX_BPLUS_TREE_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <string>
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

    template<typename T>
    class BplusTree;

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
        size_t search(const T &data) const {
            size_t i = 0;
            while (i < items.size() && !(data < items[i])) {
                i++;
            }
            return i;
        }

        // find the index of this node in its parent
        size_t find_index_in_parent() const {
            assert(parent != nullptr);
            for (size_t i = 0; i < parent->children.size(); i++) {
                if (parent->children[i] == this) {
                    return i;
                }
            }
            throw std::runtime_error("BplusTree: node not found in its parent");
        }

        bool is_underflowing() const {
            if (leaf) {
                return items.size() < ((degree + 1) / 2 - 1);
            } else {
                return children.size() < ((degree + 1) / 2 - 1);
            }
        }

        /*
           Insert a new run into this leaf. The items that the leaf should have after the insertion are stored in
           new_items, and it is on the tree to put them in the leaf, or to split or remove the leaf if they are
           too many or none. If the run reaches the previous or the next leaf it is merged with the items there.
           Returns false if the run overlaps another run. In that case nothing is changed.
        */
        bool insert(const T &data, size_t run_length, std::vector<T> &new_items) {
            assert(leaf);
            new_items.reserve(degree + 2);
            new_items.assign(items.begin(), items.end());

            size_t inserted_index;
            if (!run_insert(new_items, data, run_length, inserted_index)) {
                return false;
            }

            // the new run is at one end of the leaf, so it may touch the items of the next leaf
            if ((inserted_index >= new_items.size() || inserted_index == 0) && !new_items.empty() &&
                next != nullptr) {
                if (!merge_item_next(new_items)) {
                    return false;
                }
            }

            // the new run is at the start of the leaf, so it may touch the last item of the previous leaf
            if (inserted_index <= 1 && !new_items.empty() && prev != nullptr) {
                merge_item_prev(new_items);
            }

            return true;
        }

    private:
        bool leaf;
        size_t degree;
        std::vector<T> items;
        std::vector<bpNode<T> *> children;
        bpNode<T> *parent = nullptr;
        bpNode<T> *next = nullptr;
        bpNode<T> *prev = nullptr;

        static bool same_position(const T &a, const T &b) {
            return a.graph_position.value == b.graph_position.value;
        }

        /*
           Insert a new run into the items of this leaf and merge it with the runs it touches in there. The items
           of the other leaves are not changed. Returns false if the new run overlaps another run. inserted_index
           is set to the index the run was inserted at, or to the number of items if it was added at the end.
        */
        bool run_insert(std::vector<T> &new_items, const T &new_run, size_t run_length, size_t &inserted_index) {
            size_t len = new_items.size();
            size_t run_end = new_run.start_position + run_length;
            T gap = T::gap(run_end);

            // index is the first item that starts after the start of the new run and gap_index is the first item
            // that starts at or after its end. If they are not the same, an item starts inside the new run.
            size_t index = std::upper_bound(new_items.begin(), new_items.end(), new_run) - new_items.begin();
            size_t gap_index = std::lower_bound(new_items.begin(), new_items.end(), gap) - new_items.begin();
            if (gap_index > index) {
                return false;
            }

            if (new_items.empty()) {
                new_items.push_back(new_run);
                new_items.push_back(gap);
                inserted_index = 0;
                return true;
            }

            if (index == 0) {
                // the new run starts before all of the items here, so the previous leaf has to end with a gap
                if (prev != nullptr && !prev->get_item(prev->get_size() - 1).is_gap()) {
                    return false;
                }
            } else {
                // a new run can only start in a gap
                if (!new_items[index - 1].is_gap()) {
                    return false;
                }
                // and it should not go into the next leaf
                if (index == len && next != nullptr && run_end > next->get_item(0).start_position) {
                    return false;
                }
            }

            // Does the new run start right where the previous run ends (at its gap), and does it end right
            // where the next item starts?
            bool hits_prev = index != 0 && new_run.start_position == new_items[index - 1].start_position;
            bool hits_next = index != len && run_end == new_items[index].start_position;

            if (!hits_prev && !hits_next) {
                // first case, the new run doesn't touch anything. Add it with its own gap
                new_items.insert(new_items.begin() + index, new_run);
                new_items.insert(new_items.begin() + index + 1, gap);
                inserted_index = (index == len) ? new_items.size() : index;

            } else if (hits_prev && !hits_next) {
                // second case, the new run comes right after the previous run
                if (index > 1 && same_position(new_run, new_items[index - 2])) {
                    // same graph position, so the previous run just gets longer and we only move its gap
                    new_items[index - 1].start_position += run_length;
                } else {
                    // put the new run in place of the gap and add a new gap after it
                    new_items[index - 1] = new_run;
                    new_items.insert(new_items.begin() + index, gap);
                }
                // If the gap was the first item of the leaf, the previous run is in the previous leaf. The index
                // has to stay 0 in this case so that the new run gets checked against that leaf.
                inserted_index = (index == len && index > 1) ? new_items.size() : index - 1;

            } else if (!hits_prev && hits_next) {
                // third case, the next run starts right where the new run ends, so no gap is needed
                if (same_position(new_run, new_items[index])) {
                    // same graph position, the next run now starts earlier
                    new_items[index].start_position = new_run.start_position;
                } else {
                    new_items.insert(new_items.begin() + index, new_run);
                }
                inserted_index = index;

            } else {
                // fourth case, the new run fills the whole gap between the previous run and the next run
                bool same_as_prev = index > 1 && same_position(new_run, new_items[index - 2]);
                bool same_as_next = same_position(new_run, new_items[index]);

                if (same_as_prev && same_as_next) {
                    // the three runs become one, remove the gap and the next run
                    new_items.erase(new_items.begin() + index - 1, new_items.begin() + index + 1);
                } else if (same_as_prev) {
                    // the previous run gets longer, just remove the gap
                    new_items.erase(new_items.begin() + index - 1);
                } else if (same_as_next) {
                    // the next run now starts earlier, and the gap is removed
                    new_items[index].start_position = new_run.start_position;
                    new_items.erase(new_items.begin() + index - 1);
                } else {
                    // nothing to merge, the new run takes the place of the gap
                    new_items[index - 1] = new_run;
                }
                inserted_index = index - 1;
            }

            return true;
        }

        // Called when new_items may end right where the next leaf starts. If new_items ends with a gap that
        // starts at the first item of the next leaf, the gap is not needed anymore. If that first item also has
        // the same graph position as our last run, it is merged into our run and removed from the next leaf.
        // Returns false if the gap goes past the start of the next leaf.
        bool merge_item_next(std::vector<T> &new_items) {
            assert(next != nullptr);
            const T &last = new_items.back();
            if (!last.is_gap()) {
                return true;
            }

            T next_first = next->get_item(0);
            if (next_first.start_position < last.start_position) {
                return false;
            }
            if (next_first.start_position == last.start_position) {
                if (new_items.size() > 1 && same_position(next_first, new_items[new_items.size() - 2])) {
                    next->remove_first_item();
                }
                new_items.pop_back();
            }
            return true;
        }

        // If the previous leaf ends with the same graph position that new_items starts with, the first item of
        // new_items just continues that one (a run after the same run, or a gap after a gap) and is removed.
        void merge_item_prev(std::vector<T> &new_items) {
            assert(prev != nullptr);
            assert(!new_items.empty());
            if (same_position(prev->get_item(prev->get_size() - 1), new_items[0])) {
                new_items.erase(new_items.begin());
            }
        }
    };


    template<typename T>
    class BplusTree {
    public:
        // The degree is the maximum number of items in a node. It has to be at least 5, so that a leaf always
        // has two items or more.
        explicit BplusTree(size_t degree) : root(nullptr), degree(degree) {
            if (degree < 5) {
                throw std::invalid_argument("BplusTree: the degree has to be at least 5");
            }
        }

        ~BplusTree() {
            destroy(root);
        }

        // The tree owns its nodes through raw pointers, so it cannot be copied.
        BplusTree(const BplusTree &) = delete;
        BplusTree &operator=(const BplusTree &) = delete;

        // Insert a run with the given length. Returns false if the run overlaps a run that is already in the
        // tree, and the tree is not changed in that case.
        bool insert(const T &data, size_t run_length) {
            if (root == nullptr) {
                root = new bpNode<T>(degree, true);
            }

            bpNode<T> *leaf = leaf_search(data);
            std::vector<T> new_items;
            if (!leaf->insert(data, run_length, new_items)) {
                return false;
            }

            if (new_items.empty()) {
                // all of the items of this leaf got merged into the leaves around it, so remove the leaf
                if (leaf->get_next() != nullptr) {
                    leaf->get_next()->set_prev(leaf->get_prev());
                }
                if (leaf->get_prev() != nullptr) {
                    leaf->get_prev()->set_next(leaf->get_next());
                }
                remove_from_parent(leaf);
                delete leaf;
                return true;
            }

            if (new_items.size() <= degree) {
                leaf->replace_items(new_items);
            } else {
                split_leaf(leaf, new_items);
            }

            // merging with the next leaf can take an item away from it
            if (leaf->get_next() != nullptr && leaf->get_next()->is_underflowing()) {
                leaf_underflow(leaf->get_next());
            }
            if (leaf->is_underflowing() && leaf->get_size() != 0) {
                leaf_underflow(leaf);
            }

            return true;
        }

        // Returns the last item that starts at or before the position. This is the run that covers the position,
        // or a gap if no run covers it.
        T search(size_t position) const {
            if (root == nullptr) {
                return T::gap(0);
            }

            T data = T::gap(position);
            bpNode<T> *leaf = leaf_search(data);
            size_t index = leaf->search(data);
            if (index > 0) {
                return leaf->get_item(index - 1);
            }
            if (leaf->get_prev() != nullptr) {
                return leaf->get_prev()->get_item(leaf->get_prev()->get_size() - 1);
            }
            return T::gap(0);
        }

        // The number of items in the tree, with the gaps.
        size_t get_bpt_size() const {
            size_t size = 0;
            for (bpNode<T> *leaf = first_leaf(); leaf != nullptr; leaf = leaf->get_next()) {
                size += leaf->get_size();
            }
            return size;
        }

        // The number of runs in the tree, without the gaps.
        size_t get_bpt_run_count() const {
            size_t count = 0;
            for (bpNode<T> *leaf = first_leaf(); leaf != nullptr; leaf = leaf->get_next()) {
                for (const T &item: leaf->get_items()) {
                    if (!item.is_gap()) {
                        count++;
                    }
                }
            }
            return count;
        }

        // Print the tree, one node per line. Only for debugging.
        void print(std::ostream &out) const {
            print(out, root, 0);
        }

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

        void destroy(bpNode<T> *node) {
            if (node == nullptr) {
                return;
            }
            for (bpNode<T> *child: node->get_children()) {
                destroy(child);
            }
            delete node;
        }

        // The leaf that the data belongs to.
        bpNode<T> *leaf_search(const T &data) const {
            assert(root != nullptr);
            bpNode<T> *current = root;
            while (!current->is_leaf()) {
                current = current->get_child(current->search(data));
            }
            return current;
        }

        bpNode<T> *first_leaf() const {
            bpNode<T> *current = root;
            while (current != nullptr && !current->is_leaf()) {
                current = current->get_child(0);
            }
            return current;
        }

        // Make a new root over two nodes, where data is the key between them.
        void new_root(bpNode<T> *left, const T &data, bpNode<T> *right) {
            bpNode<T> *node = new bpNode<T>(degree, false);
            node->add_child(left);
            node->add_child(right);
            node->add_item(data);

            left->set_parent(node);
            right->set_parent(node);
            root = node;
        }

        // Split a leaf that got too many items. The leaf keeps the first half of the items and a new leaf right
        // after it gets the second half.
        void split_leaf(bpNode<T> *leaf, const std::vector<T> &items) {
            size_t mid = items.size() / 2;
            bpNode<T> *new_node = new bpNode<T>(degree, true);
            leaf->replace_items(std::vector<T>(items.begin(), items.begin() + mid));
            new_node->replace_items(std::vector<T>(items.begin() + mid, items.end()));

            new_node->set_next(leaf->get_next());
            new_node->set_prev(leaf);
            leaf->set_next(new_node);
            if (new_node->get_next() != nullptr) {
                new_node->get_next()->set_prev(new_node);
            }

            T parent_item = new_node->get_item(0);
            if (leaf->is_root()) {
                new_root(leaf, parent_item, new_node);
            } else {
                parent_insert(leaf->get_parent(), parent_item, new_node);
            }
        }

        // Add a new key and the child that comes after it to an internal node. The node is split if it is full.
        void parent_insert(bpNode<T> *node, const T &data, bpNode<T> *child) {
            size_t add_index = node->search(data);

            if (node->get_size() < degree) {
                child->set_parent(node);
                node->add_item(data, add_index);
                node->add_child(child, add_index + 1);
                return;
            }

            // The node is full. Put all of the keys and children together with the new ones, the node keeps the
            // first half, the key in the middle goes up to the parent and a new node gets the rest.
            std::vector<T> all_items = node->get_items();
            std::vector<bpNode<T> *> all_children = node->get_children();
            all_items.insert(all_items.begin() + add_index, data);
            all_children.insert(all_children.begin() + add_index + 1, child);

            size_t mid = degree / 2 + 1;
            T parent_item = all_items[mid];

            bpNode<T> *new_node = new bpNode<T>(degree, false);
            node->replace_items(std::vector<T>(all_items.begin(), all_items.begin() + mid));
            node->replace_children(
                    std::vector<bpNode<T> *>(all_children.begin(), all_children.begin() + mid + 1));
            new_node->replace_items(std::vector<T>(all_items.begin() + mid + 1, all_items.end()));
            new_node->replace_children(
                    std::vector<bpNode<T> *>(all_children.begin() + mid + 1, all_children.end()));

            for (bpNode<T> *moved: node->get_children()) {
                moved->set_parent(node);
            }
            for (bpNode<T> *moved: new_node->get_children()) {
                moved->set_parent(new_node);
            }

            if (node->is_root()) {
                new_root(node, parent_item, new_node);
            } else {
                parent_insert(node->get_parent(), parent_item, new_node);
            }
        }

        // Remove a node from its parent, with the key that goes with it. The node itself is not deleted here.
        // This can leave the parent with too few children, which is handled here too.
        void remove_from_parent(bpNode<T> *node) {
            bpNode<T> *parent = node->get_parent();
            assert(parent != nullptr);

            size_t index = node->find_index_in_parent();
            parent->remove_item(index == 0 ? 0 : index - 1);
            parent->remove_child(index);

            // if the root is left with only one child, that child is the new root
            if (parent->is_root() && parent->get_size() == 0) {
                root = parent->get_child(0);
                root->set_parent(nullptr);
                root->set_next(nullptr);
                root->set_prev(nullptr);

                delete parent;
                return;
            }

            if (parent->is_underflowing()) {
                handle_parent_underflow(parent);
            }
        }

        // Handles the underflow of an internal node. It borrows a child from the left or the right sibling if
        // they have enough, and merges with one of them otherwise.
        void handle_parent_underflow(bpNode<T> *node) {
            assert(!node->is_leaf());

            // the root can have any number of children
            if (node->is_root()) {
                return;
            }

            bpNode<T> *parent = node->get_parent();
            size_t index = node->find_index_in_parent();

            if (index > 0 && parent->get_child(index - 1)->get_size() > degree / 2) {
                // Borrow from the left sibling
                bpNode<T> *left_sibling = parent->get_child(index - 1);
                T borrowed_item = left_sibling->get_item(left_sibling->get_size() - 1);
                bpNode<T> *borrowed_child = left_sibling->get_child(left_sibling->get_size());

                left_sibling->remove_last_item();
                left_sibling->remove_last_child();

                node->add_item(parent->get_item(index - 1), 0);
                node->add_child(borrowed_child, 0);
                borrowed_child->set_parent(node);

                parent->change_item(borrowed_item, index - 1);

            } else if (index < parent->get_size() && parent->get_child(index + 1)->get_size() > degree / 2) {
                // Borrow from the right sibling
                bpNode<T> *right_sibling = parent->get_child(index + 1);
                T borrowed_item = right_sibling->get_item(0);
                bpNode<T> *borrowed_child = right_sibling->get_child(0);

                right_sibling->remove_first_item();
                right_sibling->remove_child(0);

                node->add_item(parent->get_item(index));
                node->add_child(borrowed_child);
                borrowed_child->set_parent(node);

                parent->change_item(borrowed_item, index);

            } else if (index > 0) {
                // Merge with the left sibling
                bpNode<T> *left_sibling = parent->get_child(index - 1);
                left_sibling->add_item(parent->get_item(index - 1));
                for (const T &item: node->get_items()) {
                    left_sibling->add_item(item);
                }
                for (bpNode<T> *child: node->get_children()) {
                    left_sibling->add_child(child);
                    child->set_parent(left_sibling);
                }

                remove_from_parent(node);
                delete node;

            } else if (index < parent->get_size()) {
                // Merge with the right sibling
                bpNode<T> *right_sibling = parent->get_child(index + 1);
                node->add_item(parent->get_item(index));
                for (const T &item: right_sibling->get_items()) {
                    node->add_item(item);
                }
                for (bpNode<T> *child: right_sibling->get_children()) {
                    node->add_child(child);
                    child->set_parent(node);
                }

                remove_from_parent(right_sibling);
                delete right_sibling;
            }
        }

        // Handles the underflow of a leaf. It borrows an item from the left or the right sibling if they have
        // enough, and merges with one of them otherwise.
        void leaf_underflow(bpNode<T> *node) {
            assert(node->is_leaf());

            // the tree has only one node, nothing to do
            if (node->is_root()) {
                return;
            }

            bpNode<T> *parent = node->get_parent();
            size_t index = node->find_index_in_parent();
            size_t min_items = (degree + 1) / 2 - 1;

            if (index > 0 && parent->get_child(index - 1)->get_size() > min_items) {
                // borrow from the left sibling
                bpNode<T> *left_sibling = parent->get_child(index - 1);
                assert(left_sibling->get_next() == node);
                assert(node->get_prev() == left_sibling);

                T borrowed_item = left_sibling->get_item(left_sibling->get_size() - 1);
                node->add_item(borrowed_item, 0);
                left_sibling->remove_last_item();
                parent->change_item(borrowed_item, index - 1);

            } else if (index < parent->get_size() && parent->get_child(index + 1)->get_size() > min_items) {
                // borrow from the right sibling
                bpNode<T> *right_sibling = parent->get_child(index + 1);
                assert(right_sibling->get_prev() == node);
                assert(node->get_next() == right_sibling);

                node->add_item(right_sibling->get_item(0));
                right_sibling->remove_first_item();
                parent->change_item(right_sibling->get_item(0), index);

            } else if (index > 0) {
                // merge with the left sibling
                bpNode<T> *left_sibling = parent->get_child(index - 1);
                assert(left_sibling->get_next() == node);
                assert(node->get_prev() == left_sibling);

                for (const T &item: node->get_items()) {
                    left_sibling->add_item(item);
                }
                left_sibling->set_next(node->get_next());
                if (node->get_next() != nullptr) {
                    node->get_next()->set_prev(left_sibling);
                }

                remove_from_parent(node);
                delete node;

            } else if (index < parent->get_size()) {
                // merge with the right sibling
                bpNode<T> *right_sibling = parent->get_child(index + 1);
                assert(right_sibling->get_prev() == node);
                assert(node->get_next() == right_sibling);

                for (const T &item: right_sibling->get_items()) {
                    node->add_item(item);
                }
                node->set_next(right_sibling->get_next());
                if (node->get_next() != nullptr) {
                    node->get_next()->set_prev(node);
                }

                remove_from_parent(right_sibling);
                delete right_sibling;
            }
        }

        void print(std::ostream &out, const bpNode<T> *node, size_t level) const {
            if (node == nullptr) {
                return;
            }

            out << std::string(level * 4, ' ') << (node->is_leaf() ? "Leaf: " : "Internal: ");
            for (const T &item: node->get_items()) {
                out << "[" << item.start_position << ":" << item.graph_position.value << "] ";
            }
            out << (node->is_root() ? "(root)" : "") << "\n";

            for (const bpNode<T> *child: node->get_children()) {
                print(out, child, level + 1);
            }
        }
    };

}

#endif // PANGENOME_INDEX_BPLUS_TREE_HPP
