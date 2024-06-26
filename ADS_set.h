#ifndef ADS_SET_ADS_SET_H
#define ADS_SET_ADS_SET_H

#include <initializer_list>
#include <iostream>
#include <utility>

#define NDEBUG

#include <cassert>

template<class Type, uint32_t BlockSize>
class NodeDataBlock {
public:
    using key_type = Type;
    using key_compare_values = std::less<key_type>;
    using key_equal_values = std::equal_to<key_type>;
    using key_compare_pointers = std::less<std::remove_pointer_t<key_type> >;
    using key_equal_pointers = std::equal_to<std::remove_pointer_t<key_type> >;

    using size_type = std::uint32_t;

    class View {
    public:
        View(const key_type *const parr_, const size_type sz_) : parr{parr_}, sz{sz_} {}

        [[nodiscard]] inline const key_type *begin() const noexcept {
            return parr;
        }

        [[nodiscard]] inline const key_type *end() const noexcept {
            return parr + sz;
        }

        [[nodiscard]] inline size_type size() const noexcept {
            return sz;
        }

        [[nodiscard]] inline const key_type &operator[](const size_type idx) const noexcept {
            assert(idx < sz);
            return parr[idx];
        }

    private:
        const key_type *parr = nullptr;
        size_type sz = 0;
    };

    [[nodiscard]] inline operator View() const noexcept {  // NOLINT(google-explicit-constructor)
        return {block, sz};
    }

    [[nodiscard]] inline const key_type *begin() const noexcept {
        return block;
    }

    [[nodiscard]] inline const key_type *end() const noexcept {
        return block + sz;
    }

    [[nodiscard]] inline const key_type &operator[](const size_type idx) const noexcept {
        assert(idx < sz && !empty());
        return block[idx];
    }

    [[nodiscard]] inline key_type &operator[](const size_type idx) noexcept {
        assert(idx < sz && !empty());
        return block[idx];
    }

    [[nodiscard]] inline const key_type &front() const noexcept {
        assert(!empty());
        return block[0];
    }

    [[nodiscard]] inline const key_type &back() const noexcept {
        assert(!empty());
        return block[sz - 1];
    }

    [[nodiscard]] inline size_type size() const noexcept {
        return sz;
    }

    [[nodiscard]] inline bool empty() const noexcept {
        return sz == 0;
    }

    [[nodiscard]] inline bool full() const noexcept {
        return sz >= BlockSize;
    }

    void insertAndShiftRight(const size_type idx, const key_type &value) noexcept {
        assert(idx < BlockSize && !full());
        for (size_type i = sz; i > idx; --i)
#pragma GCC diagnostic warning "-Warray-bounds"
                block[i] = std::move(block[i - 1]); // no overflow because now i >= 1
#pragma GCC diagnostic pop
        block[idx] = value;
        ++sz;
    }

    inline void pushFront(const key_type &value) noexcept {
        insertAndShiftRight(0, value);
    }

    inline void pushBack(const key_type &value) noexcept {
        insertAndShiftRight(sz, value);
    }

    void eraseAndShiftLeft(const size_type idx) noexcept {
        assert(idx < BlockSize && !empty());
        for (size_type i = idx; i < sz - 1; ++i)   // no overflow because now sz >= 1
            block[i] = std::move(block[i + 1]);
        --sz;
    }

    inline void popFront() noexcept {
        eraseAndShiftLeft(0);
    }

    inline void popBack() noexcept {
        eraseAndShiftLeft(sz - 1);
    }

    NodeDataBlock splitAddNewValue(const size_type newValueIdx, const key_type &newValue) noexcept {
        assert(newValueIdx <= BlockSize && full());
        // oldBlockSize >= newBlock_size
        size_type oldBlockSize = (BlockSize + 1) / 2;
        NodeDataBlock newBlock;
        if (newValueIdx < oldBlockSize) {  // insert value in old block
            --oldBlockSize;
            newBlock.sz = sz - oldBlockSize;
            sz = oldBlockSize;
            for (size_type i = 0; i < newBlock.sz; ++i)
                newBlock.block[i] = std::move(block[i + oldBlockSize]);
            insertAndShiftRight(newValueIdx, newValue);
        } else {  // insert value in new block
            newBlock.sz = sz - oldBlockSize + 1;
            sz = oldBlockSize;
            for (size_type i = 0; i < newValueIdx - oldBlockSize; ++i)
                newBlock.block[i] = std::move(block[i + oldBlockSize]);
            newBlock.block[newValueIdx - oldBlockSize] = newValue;
            for (size_type i = newValueIdx - oldBlockSize + 1; i < newBlock.sz; ++i)
                newBlock.block[i] = std::move(block[i + oldBlockSize - 1]);
        }
        return newBlock;
    }

    void merge(View rightBlock) noexcept {
        assert(size() + rightBlock.size() <= BlockSize);
        for (size_type i = 0; i < rightBlock.size(); ++i)
            block[i + sz] = std::move(rightBlock[i]);
        sz += rightBlock.size();
    }

    // only for sorted values, SFINAE
    // returns (true, idx) if contains value, else (false, idx) of first element, so that value <= element
    [[nodiscard]] std::pair<bool, size_type> contains(const std::remove_pointer_t<key_type> &value) const noexcept {
        if constexpr (std::is_pointer_v<key_type>) {
            if constexpr (BlockSize >= 256) {   // binary search
                const key_type *it = binarySearch(begin(), end(), value);  // lower bound, value <= *it
                const size_type idx = it - begin();
                // contains value
                return {it != end() && key_equal_pointers{}(value, *(*it)), idx};
            } else {   // linear search
                for (size_type idx = 0; idx < sz; ++idx) {
                    if (key_compare_pointers{}(*block[idx], value))   // value >= block[idx]
                        continue;
                    if (key_equal_pointers{}(value, *block[idx]))
                        return {true, idx};
                    if (key_compare_pointers{}(value, *block[idx]))
                        return {false, idx};
                }
                return {false, sz};
            }
        } else {
            if constexpr (BlockSize >= 256) {   // binary search
                const key_type *it = binarySearch(begin(), end(), value);  // lower bound, value <= *it
                const size_type idx = it - begin();
                // contains value
                return {it != end() && key_equal_values{}(value, *it), idx};
            } else {   // linear search
                for (size_type idx = 0; idx < sz; ++idx) {
                    if (key_compare_values{}(block[idx], value))   // value >= block[idx]
                        continue;
                    if (key_equal_values{}(value, block[idx]))
                        return {true, idx};
                    if (key_compare_values{}(value, block[idx]))
                        return {false, idx};
                }
                return {false, sz};
            }
        }
    }

    inline void clear() noexcept {
        sz = 0;
    }

private:
    key_type block[BlockSize]{};
    size_type sz = 0;

    // only for sorted values, SFINAE, lower bound
    static const key_type *
    binarySearch(const key_type *first, const key_type *last, const std::remove_pointer_t<key_type> &value) {
        const key_type *it;
        for (ptrdiff_t count = last - first, step; count > 0;) {
            it = first;
            step = count / 2;
            it += step;
            if constexpr (std::is_pointer_v<key_type>) {
                if (key_compare_pointers{}(*(*it), value)) {
                    first = ++it;
                    count -= step + 1;
                } else
                    count = step;
            } else {
                if (key_compare_values{}(*it, value)) {
                    first = ++it;
                    count -= step + 1;
                } else
                    count = step;
            }
        }
        return first;
    }
};

template<class Type>
class Stack {
public:
    using size_type = std::size_t;

    explicit Stack(const size_type max_sz_) noexcept: data{new Type[max_sz_]}, max_sz{max_sz_}, sz{0} {}

    Stack(const Stack &other) noexcept: data{new Type[other.max_sz]}, max_sz{other.max_sz}, sz{other.sz} {
        std::copy(other.data, other.data + other.sz, data);
    }

    Stack(Stack &&other) noexcept:
            data{std::exchange(other.data, nullptr)}, max_sz{std::exchange(other.max_sz, 0)},
            sz{std::exchange(other.sz, 0)} {}

    // copy and move assignment
    Stack &operator=(Stack other) noexcept {
        std::swap(data, other.data);
        std::swap(max_sz, other.max_sz);
        std::swap(sz, other.sz);
        return *this;
    }

    ~Stack() noexcept {
        delete[] data;
    }

    [[nodiscard]] inline size_type size() const noexcept {
        return sz;
    }

    [[nodiscard]] inline bool empty() const noexcept {
        return sz == 0;
    }

    [[nodiscard]] inline size_type capacity() const noexcept {
        return max_sz;
    }

    template<class... Args>
    inline void emplace(Args &&... args) {
        assert(sz < max_sz);
        data[sz++] = Type(std::forward<Args>(args)...);
    }

    inline void push(Type &&value) noexcept {
        assert(sz < max_sz);
        emplace(std::move(value));
    }

    inline void push(const Type &value) noexcept {
        assert(sz < max_sz);
        emplace(value);
    }

    inline const Type &top() const noexcept {
        assert(!empty());
        return data[sz - 1];
    }

    inline Type &top() noexcept {
        assert(!empty());
        return data[sz - 1];
    }

    inline void pop() noexcept {
        assert(!empty());
        --sz;
    }

private:
    Type *data = nullptr;
    size_type max_sz = 0;
    size_type sz = 0;
};

template<class KeyType, size_t Order>
class BaseNode;

template<class KeyType, size_t Order>
class InternalNode;

template<class KeyType, size_t Order>
class LeafNode;

template<typename KeyType, size_t Order>
struct keys_block_helper {
    using type = std::conditional_t<(sizeof(KeyType) > 16),
            NodeDataBlock<KeyType *, 2 * Order>, NodeDataBlock<KeyType, 2 * Order> >;
};

template<class KeyType, size_t Order>
class BaseNode {
public:
    using key_type = KeyType;
    using key_compare = std::less<key_type>;
    using key_equal = std::equal_to<key_type>;

    using size_type = std::size_t;

    using base_node = BaseNode<key_type, Order>;
    using internal_node = InternalNode<key_type, Order>;
    using leaf_node = LeafNode<key_type, Order>;

    using keys_block = typename keys_block_helper<key_type, Order>::type;

    static_assert(Order >= 1);

    friend class InternalNode<key_type, Order>;

    friend class LeafNode<key_type, Order>;

    const bool isLeaf = true;

    explicit BaseNode(const bool isLeaf_) : isLeaf{isLeaf_} {}

    virtual ~BaseNode() noexcept = default;

    [[nodiscard]] base_node *createNewRootNodeFromOld(base_node &newChild) noexcept {
        auto newRootNode = new internal_node;
        newRootNode->children.pushFront(this);
        newRootNode->insertChild(0, newChild);
        return newRootNode;
    }
};

template<class KeyType, size_t Order>
class InternalNode : public BaseNode<KeyType, Order> {
public:
    using key_type = KeyType;
    using key_compare = std::less<key_type>;
    using key_equal = std::equal_to<key_type>;

    using size_type = std::size_t;

    using base_node = BaseNode<KeyType, Order>;
    using internal_node = InternalNode<key_type, Order>;
    using leaf_node = LeafNode<KeyType, Order>;

    using keys_block = typename keys_block_helper<key_type, Order>::type;

    friend class BaseNode<KeyType, Order>;

    friend class LeafNode<KeyType, Order>;

    InternalNode() noexcept: base_node{false} {}

    ~InternalNode() noexcept override {
        for (const auto child: children)
            delete child;
    }

    friend std::ostream &operator<<(std::ostream &os, const internal_node &node) noexcept {
        node.print(os);
        return os;
    }

    void print(std::ostream &os) const noexcept {
        os << "InternalNode: this = " << this << " keys (size " << keys.size() << "): ";
        if constexpr (sizeof(KeyType) > 16) {
            for (const auto &key: keys)
                os << *key << ' ';
        } else {
            for (const auto &key: keys)
                os << key << ' ';
        }
        os << " children (size " << children.size() << "): ";
        for (const auto child: children)
            os << child << ' ';
    }

    [[nodiscard]] inline bool needsMerging() const noexcept {
        return keys.size() < Order;
    }

    [[nodiscard]] inline bool hasEnoughKeysToBorrow() const noexcept {
        assert(keys.size() >= Order - 1);
        return keys.size() > Order;
    }

    [[nodiscard]] inline const keys_block &getKeys() const noexcept {
        return keys;
    }

    // returns nullptr if was split
    [[nodiscard]] internal_node *addChildToInternalNode(const size_type idx, base_node &newChild) noexcept {
        if (!keys.full()) {
            insertChild(idx, newChild);
            return nullptr;
        }
        internal_node *newRightINode = splitInternalNode(idx, newChild);
        return newRightINode;
    }

    [[nodiscard]] inline const NodeDataBlock<base_node *, 2 * Order + 1> &getChildren() const noexcept {
        return children;
    }

    void fixProblemChild(const size_type problemChildIdx) noexcept {
        assert(!keys.empty() && !children.empty());
        if (children[problemChildIdx]->isLeaf) {
            if (problemChildIdx > 0 && static_cast<leaf_node &>(*children[problemChildIdx - 1]).hasEnoughKeysToBorrow())
                borrowKeyFromLeftSibling(problemChildIdx);
            else if (problemChildIdx < children.size() - 1 &&
                     static_cast<leaf_node &>(*children[problemChildIdx + 1]).hasEnoughKeysToBorrow())
                borrowKeyFromRightSibling(problemChildIdx);
            else {
                size_type leftChildIdx = std::min(problemChildIdx, static_cast<size_type>(keys.size()) - 1);
                mergeWithRightSibling(leftChildIdx);
            }
        } else {
            if (problemChildIdx > 0 &&
                static_cast<internal_node &>(*children[problemChildIdx - 1]).hasEnoughKeysToBorrow())
                borrowKeyFromLeftSibling(problemChildIdx);
            else if (problemChildIdx < children.size() - 1 &&
                     static_cast<internal_node &>(*children[problemChildIdx + 1]).hasEnoughKeysToBorrow())
                borrowKeyFromRightSibling(problemChildIdx);
            else {
                size_type leftChildIdx = std::min(problemChildIdx, static_cast<size_type>(keys.size()) - 1);
                mergeWithRightSibling(leftChildIdx);
            }
        }
    }

    base_node *pullUpChildToRootNode() noexcept {
        assert(keys.empty() && children.size() == 1);
        base_node *newRootNode = children.front();
        children.clear();   // nullify children
        return newRootNode;
    }

private:
    keys_block keys;
    NodeDataBlock<base_node *, 2 * Order + 1> children;

    void insertChild(const size_type idx, base_node &newChild) noexcept {
        assert(!keys.full() && !children.full() && !children.empty() &&
               ((newChild.isLeaf && !static_cast<leaf_node &>(newChild).keys.empty()) ||
                (!newChild.isLeaf && !static_cast<internal_node &>(newChild).keys.empty())));
        assert(keys.contains(newChild.isLeaf ? static_cast<leaf_node &>(newChild).keys.front()
                                             : static_cast<internal_node &>(newChild).keys.front()).second == idx);
        if (newChild.isLeaf) {
            if constexpr (sizeof(KeyType) > 16)
                keys.insertAndShiftRight(idx, const_cast<key_type *>(static_cast<leaf_node &>(newChild).keys.begin()));
            else
                keys.insertAndShiftRight(idx, static_cast<leaf_node &>(newChild).keys.front());
            children.insertAndShiftRight(idx + 1, &newChild);
        } else {
            keys.insertAndShiftRight(idx, static_cast<internal_node &>(newChild).keys.front());
            children.insertAndShiftRight(idx + 1, &newChild);
            // deletes the temporary first key of new right internal node child
            static_cast<internal_node &>(newChild).keys.popFront();
        }
    }

    [[nodiscard]] internal_node *splitInternalNode(const size_type idx, base_node &newChild) noexcept {
        assert(keys.full() && children.full());
        auto newRightInternalNode = new internal_node;
        if (newChild.isLeaf) {
            // the first key of new right internal node is saved in order to delete it later in insertChild or splitInternalNode parent method
            if constexpr (sizeof(KeyType) > 16)
                newRightInternalNode->keys = keys.splitAddNewValue(idx,
                                                                   const_cast<key_type *>(static_cast<leaf_node &>(newChild).keys.begin()));
            else
                newRightInternalNode->keys = keys.splitAddNewValue(idx,
                                                                   static_cast<leaf_node &>(newChild).keys.front());
            newRightInternalNode->children = children.splitAddNewValue(idx + 1, &newChild);
        } else {
            // the first key of new right internal node is saved in order to delete it later in insertChild or splitInternalNode parent method
            newRightInternalNode->keys = keys.splitAddNewValue(idx,
                                                               static_cast<internal_node &>(newChild).keys.front());
            newRightInternalNode->children = children.splitAddNewValue(idx + 1, &newChild);
            // deletes the temporary first key of new right internal node child
            static_cast<internal_node &>(newChild).keys.popFront();
        }
        return newRightInternalNode;
    }

    void borrowKeyFromLeftSibling(const size_type problemChildIdx) noexcept {
        assert(problemChildIdx > 0 &&
               ((children[problemChildIdx]->isLeaf &&
                 static_cast<leaf_node &>(*children[problemChildIdx - 1]).hasEnoughKeysToBorrow()) ||
                (!children[problemChildIdx]->isLeaf &&
                 static_cast<internal_node &>(*children[problemChildIdx - 1]).hasEnoughKeysToBorrow())));
        if (children[problemChildIdx]->isLeaf) {
            auto &leftSibling = static_cast<leaf_node &>(*children[problemChildIdx - 1]);
            auto &rightProblemChild = static_cast<leaf_node &>(*children[problemChildIdx]);
            rightProblemChild.keys.pushFront(leftSibling.keys.back());
            leftSibling.keys.popBack();
            // update index key
            if constexpr (sizeof(KeyType) > 16)
                keys[problemChildIdx - 1] = const_cast<key_type *>(rightProblemChild.keys.begin());
            else
                keys[problemChildIdx - 1] = rightProblemChild.keys.front();
        } else {
            auto &leftSibling = static_cast<internal_node &>(*children[problemChildIdx - 1]);
            auto &rightProblemChild = static_cast<internal_node &>(*children[problemChildIdx]);
            // left sibling last key ->/ left to right problem child parent key \-> right problem child first key
            rightProblemChild.keys.pushFront(keys[problemChildIdx - 1]);  // move down index key
            keys[problemChildIdx - 1] = leftSibling.keys.back();  // update index key
            rightProblemChild.children.pushFront(leftSibling.children.back());
            leftSibling.keys.popBack();
            leftSibling.children.popBack();
        }
    }

    void borrowKeyFromRightSibling(const size_type problemChildIdx) noexcept {
        assert(problemChildIdx < children.size() - 1 &&
               ((children[problemChildIdx]->isLeaf &&
                 static_cast<leaf_node &>(*children[problemChildIdx + 1]).hasEnoughKeysToBorrow()) ||
                (!children[problemChildIdx]->isLeaf &&
                 static_cast<internal_node &>(*children[problemChildIdx + 1]).hasEnoughKeysToBorrow())));
        if (children[problemChildIdx]->isLeaf) {
            auto &leftProblemChild = static_cast<leaf_node &>(*children[problemChildIdx]);
            auto &rightSibling = static_cast<leaf_node &>(*children[problemChildIdx + 1]);
            leftProblemChild.keys.pushBack(rightSibling.keys.front());
            rightSibling.keys.popFront();
            // update index key
            if constexpr (sizeof(KeyType) > 16)
                keys[problemChildIdx] = const_cast<key_type *>(rightSibling.keys.begin());
            else
                keys[problemChildIdx] = rightSibling.keys.front();
        } else {
            auto &leftProblemChild = static_cast<internal_node &>(*children[problemChildIdx]);
            auto &rightSibling = static_cast<internal_node &>(*children[problemChildIdx + 1]);
            // left problem child last key <-/ left to right sibling parent key \<- right sibling first key
            leftProblemChild.keys.pushBack(keys[problemChildIdx]);  // move down index key
            keys[problemChildIdx] = rightSibling.keys.front();  // update index key
            leftProblemChild.children.pushBack(rightSibling.children.front());
            rightSibling.keys.popFront();
            rightSibling.children.popFront();
        }
    }

    void mergeWithRightSibling(const size_type leftChildIdx) noexcept {
        assert(leftChildIdx < keys.size());
        if (children[leftChildIdx]->isLeaf) {
            auto &leftChild = static_cast<leaf_node &>(*children[leftChildIdx]);
            auto &rightChild = static_cast<leaf_node &>(*children[leftChildIdx + 1]);
            assert(!leftChild.hasEnoughKeysToBorrow() && !rightChild.hasEnoughKeysToBorrow());
            // left child keys + left to right sibling parent key (than deleted in parent) + right sibling keys
            leftChild.keys.merge(rightChild.keys);
            leftChild.next = rightChild.next;
            delete &rightChild;
        } else {
            auto &leftChild = static_cast<internal_node &>(*children[leftChildIdx]);
            auto &rightChild = static_cast<internal_node &>(*children[leftChildIdx + 1]);
            assert(!leftChild.hasEnoughKeysToBorrow() && !rightChild.hasEnoughKeysToBorrow());
            // left child keys + left to right sibling parent key (than deleted in parent) + right sibling keys
            leftChild.keys.pushBack(keys[leftChildIdx]);
            leftChild.keys.merge(rightChild.keys);
            leftChild.children.merge(rightChild.children);
            rightChild.children.clear();   // nullify children
            delete &rightChild;
        }
        keys.eraseAndShiftLeft(leftChildIdx);
        children.eraseAndShiftLeft(leftChildIdx + 1);
    }
};

template<class KeyType, size_t Order>
class LeafNode : public BaseNode<KeyType, Order> {
public:
    using key_type = KeyType;
    using key_compare = std::less<key_type>;
    using key_equal = std::equal_to<key_type>;

    using size_type = std::size_t;

    using base_node = BaseNode<KeyType, Order>;
    using internal_node = InternalNode<key_type, Order>;
    using leaf_node = LeafNode<KeyType, Order>;

    using keys_block = typename keys_block_helper<key_type, Order>::type;

    friend class BaseNode<KeyType, Order>;

    friend class InternalNode<KeyType, Order>;

    LeafNode() noexcept: base_node{true} {}

    friend std::ostream &operator<<(std::ostream &os, const leaf_node &node) noexcept {
        node.print(os);
        return os;
    }

    void print(std::ostream &os) const noexcept {
//        os << "LeafNode: prev = " << prev << " this = " << this << " next = " << next << " keys (size "
//           << keys.size() << "): ";
        os << "LeafNode: this = " << this << " next = " << next << " keys (size " << keys.size() << "): ";
        for (const auto &key: keys)
            os << key << ' ';
    }

    [[nodiscard]] inline bool needsMerging() const noexcept {
        return keys.size() < Order;
    }

    [[nodiscard]] inline bool hasEnoughKeysToBorrow() const noexcept {
        assert(keys.size() >= Order - 1);
        return keys.size() > Order;
    }

    [[nodiscard]] inline const NodeDataBlock<key_type, 2 * Order> &getKeys() const noexcept {
        return keys;
    }

    // returns nullptr if was split
    [[nodiscard]] leaf_node *addKeyToLeaf(const size_type idx, const key_type &key) noexcept {
        if (!keys.full()) {
            keys.insertAndShiftRight(idx, key);
            return nullptr;
        }
        leaf_node *newRightLeaf = splitLeaf(idx, key);
        return newRightLeaf;
    }

    [[nodiscard]] inline const leaf_node *getNext() const noexcept {
        return next;
    }

//    [[nodiscard]] inline const leaf_node *getPrev() const noexcept {
//        return prev;
//    }

    inline void removeKeyFromLeaf(const size_type idx) noexcept {
        keys.eraseAndShiftLeft(idx);
    }

private:
    NodeDataBlock<key_type, 2 * Order> keys;
    leaf_node *next = nullptr;

    [[nodiscard]] leaf_node *splitLeaf(const size_type idx, const key_type &newKey) noexcept {
        assert(keys.full());
        auto newRightLeaf = new leaf_node;
        newRightLeaf->keys = keys.splitAddNewValue(idx, newKey);
//        if (next != nullptr)
//            next->prev = newRightLeaf;
        newRightLeaf->next = next;
        next = newRightLeaf;
//        newRightLeaf->prev = this;
        return newRightLeaf;
    }
};

template<class Key, uint32_t N = 2>
class Iterator {
public:
    using value_type = Key;
    using reference = const value_type &;
    using pointer = const value_type *;
    using difference_type = std::ptrdiff_t;
//    using iterator_category = std::bidirectional_iterator_tag;
    using iterator_category = std::forward_iterator_tag;

    using size_type = std::size_t;

    using base_node = BaseNode<value_type, N>;
    using internal_node = InternalNode<value_type, N>;
    using leaf_node = LeafNode<value_type, N>;

    using keys_block = typename keys_block_helper<value_type, N>::type;

    Iterator() = default;

    Iterator(const leaf_node *const currentLeaf_, const size_type currentKeyIdx_) noexcept:
            currentLeaf{currentLeaf_}, currentKeyIdx{currentKeyIdx_} {}

    [[nodiscard]] inline reference operator*() const noexcept {
        return currentLeaf->getKeys()[currentKeyIdx];
    }

    [[nodiscard]] inline pointer operator->() const noexcept {
        return &currentLeaf->getKeys()[currentKeyIdx];
    }

    Iterator &operator++() noexcept {
        if (currentKeyIdx == currentLeaf->getKeys().size() - 1) {
            currentLeaf = currentLeaf->getNext();
            currentKeyIdx = 0;
        } else
            ++currentKeyIdx;
        return *this;
    }

    Iterator operator++(int) noexcept { // NOLINT(*-dcl21-cpp)
        Iterator old = *this;
        operator++();
        return old;
    }

//    Iterator &operator--() noexcept {
//        if (currentKeyIdx == 0) {
//            currentLeaf = currentLeaf->getPrev();
//            currentKeyIdx = currentLeaf->getKeys().size() - 1;
//        } else
//            --currentKeyIdx;
//        return *this;
//    }
//
//    Iterator operator--(int) noexcept { // NOLINT(*-dcl21-cpp)
//        Iterator old = *this;
//        operator--();
//        return old;
//    }

    friend inline bool operator==(const Iterator &lhs, const Iterator &rhs) noexcept {
        return lhs.currentLeaf == rhs.currentLeaf && lhs.currentKeyIdx == rhs.currentKeyIdx;
    }

    friend inline bool operator!=(const Iterator &lhs, const Iterator &rhs) noexcept {
        return lhs.currentLeaf != rhs.currentLeaf || lhs.currentKeyIdx != rhs.currentKeyIdx;
    }

private:
    const leaf_node *currentLeaf = nullptr;
    size_type currentKeyIdx = 0;
};

// B+-Tree
// C++17, no STL
template<class Key, size_t N = 2>
class ADS_set {
public:
    using value_type = Key;
    using key_type = Key;
    using reference = value_type &;
    using const_reference = const value_type &;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using const_iterator = Iterator<key_type, N>;
    using iterator = const_iterator;
    using key_compare = std::less<key_type>;   // B+-Tree
    using key_equal = std::equal_to<key_type>; // Hashing
    using hasher = std::hash<key_type>;        // Hashing

    using base_node = BaseNode<Key, N>;
    using internal_node = InternalNode<Key, N>;
    using leaf_node = LeafNode<Key, N>;

    using keys_block = typename keys_block_helper<value_type, N>::type;

    // O(1)
    ADS_set() = default;

    // O(init_size * log init_size)
    ADS_set(const std::initializer_list<key_type> init) noexcept {
        for (const auto &key: init)
            tryAddKey(key);
    }

    // O(size + init_size * log init_size)
    ADS_set &operator=(const std::initializer_list<key_type> init) noexcept {
        clear();
        for (const auto &key: init)
            tryAddKey(key);
        return *this;
    }

    // O(other_size * log other_size) (or faster)
    ADS_set(const ADS_set &other) noexcept {
        for (const auto &key: other)
            tryAddKey(key);
    }

    // O(size + other_size * log other_size) (or faster)
    ADS_set &operator=(ADS_set other) noexcept {
        swap(other);
        return *this;
    }

    // O(range_size * log range_size)
    template<typename InputIt>
    ADS_set(InputIt first, InputIt last) noexcept {
        for (auto it = first; it != last; ++it)
            tryAddKey(*it);
    }

    // O(size)
    ~ADS_set() noexcept {
        delete roodNode;
    }

    // O(1)
    [[nodiscard]] inline size_type size() const noexcept {
        return sz;
    }

    // O(1)
    [[nodiscard]] inline bool empty() const noexcept {
        return sz == 0;
    }

    // O(init_size * log (init_size + size))
    void insert(const std::initializer_list<key_type> init) noexcept {
        for (const auto &key: init)
            tryAddKey(key);
    }

    // O(log size)
    std::pair<const_iterator, bool> insert(const key_type &key) noexcept {
        Stack<std::pair<base_node *, size_type> > stack{height + 1};
        if (leafSearchWithPath(key, stack))
            return {const_iterator{static_cast<leaf_node *>(stack.top().first), stack.top().second}, false};
        addKeyWithPath(key, stack);
        return {find(key), true};
    }

    // O(range_size * log (range_size + size))
    template<class InputIt>
    void insert(InputIt first, InputIt last) noexcept {
        for (auto it = first; it != last; ++it)
            tryAddKey(*it);
    }

    //  O(size)
    void clear() noexcept {
        delete roodNode;
        roodNode = new leaf_node;
        sz = 0;
        height = 0;
    }

    // O(log size)
    inline size_type erase(const key_type &key) noexcept {
        return static_cast<size_type>(tryRemoveKey(key));  // removed
    }

    // O(1)
    inline void swap(ADS_set &other) noexcept {
        std::swap(roodNode, other.roodNode);
        std::swap(sz, other.sz);
        std::swap(height, other.height);
    }

    // O(log size))
    [[nodiscard]] inline size_type count(const key_type &key) const noexcept {
        return static_cast<size_type>(find(key) != end());  // erased
    }

    // B+-Baum O(log size)
    [[nodiscard]] inline const_iterator find(const key_type &key) const noexcept {
        leaf_node *leaf;
        size_type leafIdx;
        if (leafSearch(key, leaf, leafIdx)) // contains
            return {leaf, leafIdx};
        return end();
    }

    [[nodiscard]] inline const_iterator begin() const noexcept {
        if (empty())
            return end();
        return {getFirstLeaf(), 0};
    }

    [[nodiscard]] inline const_iterator end() const noexcept {
        return {nullptr, 0};
    }

    void dump(std::ostream &o = std::cerr) const noexcept {
        Stack<std::pair<base_node *, size_type> > stack{2 * N * height + 1};
        stack.emplace(roodNode, 0);
        o << "B+ Tree: size = " << sz << " height = " << height << '\n';
        while (!stack.empty()) {
            const auto [node, level] = stack.top();
            stack.pop();
            for (size_type i = 0; i < level; ++i)
                o << '\t';
            if (node->isLeaf)
                o << static_cast<leaf_node &>(*node) << '\n';
            else {
                o << static_cast<internal_node &>(*node) << '\n';
                for (auto child: static_cast<internal_node &>(*node).getChildren())
                    stack.emplace(child, level + 1);  // output right to left, "feature"
            }
        }
    }

    // O(size)
    [[nodiscard]] friend bool operator==(const ADS_set<key_type, N> &lhs, const ADS_set<key_type, N> &rhs) noexcept {
        if (lhs.size() != rhs.size())
            return false;
        auto endIt = lhs.end();
        for (auto lhsIt = lhs.begin(), rhsIt = rhs.begin(); lhsIt != endIt; ++lhsIt, ++rhsIt)
            if (!key_equal{}(*lhsIt, *rhsIt))
                return false;
        return true;
    }

    // O(size)
    [[nodiscard]] friend inline bool
    operator!=(const ADS_set<key_type, N> &lhs, const ADS_set<key_type, N> &rhs) noexcept {
        return !operator==(lhs, rhs);
    }

private:
    base_node *roodNode = new leaf_node;
    size_type sz = 0;
    size_type height = 0;

    inline bool tryAddKey(const key_type &key) noexcept {
        Stack<std::pair<base_node *, size_type> > stack{height + 1};
        if (leafSearchWithPath(key, stack))
            return false;
        addKeyWithPath(key, stack);
        return true;
    }

    void addKeyWithPath(const key_type &key, Stack<std::pair<base_node *, size_type> > &stack) noexcept {
        assert(!stack.empty());
        const auto [leaf, leafIdx] = stack.top();
        stack.pop();
        ++sz;
        base_node *newRightChild = static_cast<leaf_node &>(*leaf).addKeyToLeaf(leafIdx, key);
        if (newRightChild == nullptr)
            return;
        while (!stack.empty()) {
            const auto [internalNode, childIdx] = stack.top();
            stack.pop();
            newRightChild = static_cast<internal_node &>(*internalNode).addChildToInternalNode(childIdx,
                                                                                               *newRightChild);
            if (newRightChild == nullptr)
                return;
        }
        ++height;
        roodNode = roodNode->createNewRootNodeFromOld(*newRightChild);
    }

    bool tryRemoveKey(const Key &key) noexcept {
        if (empty())
            return false;
        Stack<std::pair<base_node *, size_type> > stack{height + 1};
        if (!leafSearchWithPath(key, stack))
            return false;
        const auto [leaf, leafIdx] = stack.top();
        stack.pop();
        static_cast<leaf_node &>(*leaf).removeKeyFromLeaf(leafIdx);
        --sz;
        if (roodNode->isLeaf || !static_cast<leaf_node &>(*leaf).needsMerging())
            return true;
        while (!stack.empty()) {
            const auto [internalNode, childIdx] = stack.top();
            stack.pop();
            static_cast<internal_node &>(*internalNode).fixProblemChild(childIdx);
            if (!static_cast<internal_node &>(*internalNode).needsMerging())
                return true;
        }
        if ((!roodNode->isLeaf && !static_cast<internal_node &>(*roodNode).getKeys().empty()) ||
            (roodNode->isLeaf && !static_cast<leaf_node &>(*roodNode).getKeys().empty()))
            return true;
        --height;
        base_node *newRootNode = static_cast<internal_node &>(*roodNode).pullUpChildToRootNode();
        delete roodNode;
        roodNode = newRootNode;
        return true;
    }

    // returns path stack with nodes and child indices, on the top is leaf child with key index
    [[nodiscard]] bool
    leafSearchWithPath(const key_type &key, Stack<std::pair<base_node *, size_type> > &stack) const noexcept {
        base_node *node = roodNode;
        auto [found, idx] = node->isLeaf ? static_cast<leaf_node &>(*node).getKeys().contains(key)
                                         : static_cast<internal_node &>(*node).getKeys().contains(key);
        while (!node->isLeaf) {
            // if key is found in internal node, we need to go to the right child,
            // contains returns lower bound <=, so we need to increment child idx
            if (found)
                ++idx;
            stack.emplace(node, idx);
            node = static_cast<internal_node &>(*node).getChildren()[idx];
            if (node->isLeaf) {
                std::tie(found, idx) = static_cast<leaf_node &>(*node).getKeys().contains(key);
                break;
            } else
                std::tie(found, idx) = static_cast<internal_node &>(*node).getKeys().contains(key);
        }
        stack.emplace(node, idx);
        return found;
    }

    // returns leaf child with key index
    [[nodiscard]] bool leafSearch(const key_type &key, leaf_node *&leaf, size_type &keyIdx) const noexcept {
        base_node *node = roodNode;
        auto [found, idx] = node->isLeaf ? static_cast<leaf_node &>(*node).getKeys().contains(key)
                                         : static_cast<internal_node &>(*node).getKeys().contains(key);
        while (!node->isLeaf) {
            // if key is found in internal node, we need to go to the right child,
            // contains returns lower bound <=, so we need to increment child idx
            if (found)
                ++idx;
            node = static_cast<internal_node &>(*node).getChildren()[idx];
            if (node->isLeaf) {
                std::tie(found, idx) = static_cast<leaf_node &>(*node).getKeys().contains(key);
                break;
            } else
                std::tie(found, idx) = static_cast<internal_node &>(*node).getKeys().contains(key);
        }
        leaf = static_cast<leaf_node *>(node);
        keyIdx = idx;
        return found;
    }

    [[nodiscard]] leaf_node *getFirstLeaf() const noexcept {
        base_node *node = roodNode;
        while (!node->isLeaf) {
            assert(!static_cast<internal_node &>(*node).getChildren().empty());
            node = static_cast<internal_node &>(*node).getChildren().front();
        }
        return static_cast<leaf_node *>(node);
    }
};

// O(1)
template<typename Key, size_t N>
inline void swap(ADS_set<Key, N> &lhs, ADS_set<Key, N> &rhs) noexcept {
    lhs.swap(rhs);
}

#endif
