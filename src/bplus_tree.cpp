//
// Created by seeskand on 5/22/24.
//

#include "pangenome_index/bplus_tree.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>

namespace panindexer {

    namespace {

        // Do the two items have the same graph position? All gaps are the same in this sense.
        template<typename T>
        bool same_position(const T &a, const T &b) {
            return a.graph_position.value == b.graph_position.value;
        }

    }

    // The node part: what happens inside one node.

    // A simple scan over the items, a node only has a few of them.
    template<typename T>
    size_t bpNode<T>::search(const T &data) const {
        size_t i = 0;
        while (i < items.size() && !(data < items[i])) {
            i++;
        }
        return i;
    }

    // Look for this node in the children of its parent.
    template<typename T>
    size_t bpNode<T>::find_index_in_parent() const {
        assert(parent != nullptr);
        for (size_t i = 0; i < parent->children.size(); i++) {
            if (parent->children[i] == this) {
                return i;
            }
        }
        throw std::runtime_error("BplusTree: node not found in its parent");
    }

    // A node underflows when it has less than about half of the degree. A leaf counts its items and an internal
    // node counts its children.
    template<typename T>
    bool bpNode<T>::is_underflowing() const {
        if (leaf) {
            return items.size() < ((degree + 1) / 2 - 1);
        } else {
            return children.size() < ((degree + 1) / 2 - 1);
        }
    }

    // The insert on a leaf. The run is first added to a copy of the items, and then the two ends of the leaf are
    // checked against the leaves around it.
    template<typename T>
    bool bpNode<T>::insert(const T &data, size_t run_length, std::vector<T> &new_items) {
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

    // The main part of the insertion. It finds where the run goes, rejects it if it overlaps something, and then
    // there are four cases depending on whether it touches the run before it and the run after it.
    template<typename T>
    bool bpNode<T>::run_insert(std::vector<T> &new_items, const T &new_run, size_t run_length,
                               size_t &inserted_index) {
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

    // Fix the border between this leaf and the next one.
    template<typename T>
    bool bpNode<T>::merge_item_next(std::vector<T> &new_items) {
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

    // Fix the border between the previous leaf and this one.
    template<typename T>
    void bpNode<T>::merge_item_prev(std::vector<T> &new_items) {
        assert(prev != nullptr);
        assert(!new_items.empty());
        if (same_position(prev->get_item(prev->get_size() - 1), new_items[0])) {
            new_items.erase(new_items.begin());
        }
    }


    // The tree part: finding the leaf and keeping the tree balanced.

    // The tree starts empty, the root is made at the first insert.
    template<typename T>
    BplusTree<T>::BplusTree(size_t degree) : root(nullptr), degree(degree) {
        if (degree < 5) {
            throw std::invalid_argument("BplusTree: the degree has to be at least 5");
        }
    }

    template<typename T>
    BplusTree<T>::~BplusTree() {
        destroy(root);
    }

    // Delete a node and everything under it.
    template<typename T>
    void BplusTree<T>::destroy(bpNode<T> *node) {
        if (node == nullptr) {
            return;
        }
        for (bpNode<T> *child: node->get_children()) {
            destroy(child);
        }
        delete node;
    }

    // Find the leaf for the run and let it work out its new items. Then fix the tree: the leaf is removed if it
    // has nothing left, split if it has too many items, and the underflows are handled.
    template<typename T>
    bool BplusTree<T>::insert(const T &data, size_t run_length) {
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

    // Go down to the leaf of the position and return the item just before the first one that starts after it.
    // That item can be the last one of the previous leaf.
    template<typename T>
    T BplusTree<T>::search(size_t position) const {
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

    // Walk over the leaves and count the items.
    template<typename T>
    size_t BplusTree<T>::get_bpt_size() const {
        size_t size = 0;
        for (bpNode<T> *leaf = first_leaf(); leaf != nullptr; leaf = leaf->get_next()) {
            size += leaf->get_size();
        }
        return size;
    }

    // The same walk, but the gaps are not counted.
    template<typename T>
    size_t BplusTree<T>::get_bpt_run_count() const {
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

    // Go down from the root. At each internal node the keys tell which child to take.
    template<typename T>
    bpNode<T> *BplusTree<T>::leaf_search(const T &data) const {
        assert(root != nullptr);
        bpNode<T> *current = root;
        while (!current->is_leaf()) {
            current = current->get_child(current->search(data));
        }
        return current;
    }

    // The leftmost leaf, this is where the iterator starts.
    template<typename T>
    bpNode<T> *BplusTree<T>::first_leaf() const {
        bpNode<T> *current = root;
        while (current != nullptr && !current->is_leaf()) {
            current = current->get_child(0);
        }
        return current;
    }

    // The old root was split, so the tree gets one level taller.
    template<typename T>
    void BplusTree<T>::new_root(bpNode<T> *left, const T &data, bpNode<T> *right) {
        bpNode<T> *node = new bpNode<T>(degree, false);
        node->add_child(left);
        node->add_child(right);
        node->add_item(data);

        left->set_parent(node);
        right->set_parent(node);
        root = node;
    }

    // Cut the items in half, link the new leaf after the old one, and give the first item of the new leaf to the
    // parent as the key between the two.
    template<typename T>
    void BplusTree<T>::split_leaf(bpNode<T> *leaf, const std::vector<T> &items) {
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

    // Called after a split, to add the new node to the parent. If the parent is full it is split too, and this
    // can go all the way up to the root.
    template<typename T>
    void BplusTree<T>::parent_insert(bpNode<T> *node, const T &data, bpNode<T> *child) {
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
        node->replace_children(std::vector<bpNode<T> *>(all_children.begin(), all_children.begin() + mid + 1));
        new_node->replace_items(std::vector<T>(all_items.begin() + mid + 1, all_items.end()));
        new_node->replace_children(std::vector<bpNode<T> *>(all_children.begin() + mid + 1, all_children.end()));

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

    // Take a node out of its parent. If the root is left with one child, the tree gets one level shorter.
    template<typename T>
    void BplusTree<T>::remove_from_parent(bpNode<T> *node) {
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

    // An internal node has too few children. First try to borrow one from the left sibling, then from the right
    // one, and if they have nothing to give merge with one of them. The keys move through the parent.
    template<typename T>
    void BplusTree<T>::handle_parent_underflow(bpNode<T> *node) {
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

    // The same idea for a leaf. Here the items move directly between the leaves and only the key in the parent is
    // updated.
    template<typename T>
    void BplusTree<T>::leaf_underflow(bpNode<T> *node) {
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

    template<typename T>
    void BplusTree<T>::print(std::ostream &out) const {
        print(out, root, 0);
    }

    // One line per node, indented by its level in the tree.
    template<typename T>
    void BplusTree<T>::print(std::ostream &out, const bpNode<T> *node, size_t level) const {
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


    // The tree is only used with Run.
    template class bpNode<Run>;
    template class BplusTree<Run>;

}
