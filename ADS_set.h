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
    using key_compare = std::less<key_type>;
    using key_equal = std::equal_to<key_type>;

    using size_type = std::uint32_t;

    class View {
    public:
        View(const Type *const parr_, const size_type sz_) : parr{parr_}, sz{sz_} {}

        [[nodiscard]] inline const Type *begin() const noexcept {
            return parr;
        }

        [[nodiscard]] inline const Type *end() const noexcept {
            return parr + sz;
        }

        [[nodiscard]] inline size_type size() const noexcept {
            return sz;
        }

        [[nodiscard]] inline const Type &operator[](const size_type idx) const noexcept {
            assert(idx < sz);
            return parr[idx];
        }

    private:
        const Type *parr = nullptr;
        size_type sz = 0;
    };

    [[nodiscard]] inline operator View() const noexcept {  // NOLINT(google-explicit-constructor)
        return {block, sz};
    }

    [[nodiscard]] inline const Type *begin() const noexcept {
        return block;
    }

    [[nodiscard]] inline const Type *end() const noexcept {
        return block + sz;
    }

    [[nodiscard]] inline const Type &operator[](const size_type idx) const noexcept {
        assert(idx < sz && !empty());
        return block[idx];
    }

    [[nodiscard]] inline Type &operator[](const size_type idx) noexcept {
        assert(idx < sz && !empty());
        return block[idx];
    }

    [[nodiscard]] inline const Type &front() const noexcept {
        assert(!empty());
        return block[0];
    }

    [[nodiscard]] inline const Type &back() const noexcept {
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

    void insertAndShiftRight(const size_type idx, const Type &value) noexcept {
        assert(idx < BlockSize && !full());
        for (size_type i = sz; i > idx; --i)
#pragma GCC diagnostic warning "-Warray-bounds"
                block[i] = std::move(block[i - 1]); // no overflow because now i >= 1
#pragma GCC diagnostic pop
        block[idx] = value;
        ++sz;
    }

    inline void pushFront(const Type &value) noexcept {
        insertAndShiftRight(0, value);
    }

    inline void pushBack(const Type &value) noexcept {
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

    NodeDataBlock splitAddNewValue(const size_type newValueIdx, const Type &newValue) noexcept {
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
    template<class U = Type, class = std::enable_if_t<!std::is_pointer_v<U> > >
    [[nodiscard]] std::pair<bool, size_type> contains(const Type &value) const noexcept {
        if constexpr (BlockSize >= 256) {   // binary search
            const Type *it = binarySearch(begin(), end(), value);  // lower bound, value <= *it
            const size_type idx = it - begin();
            if (it != end() && key_equal{}(value, *it))  // contains value
                return {true, idx};
            return {false, idx};
        } else {   // linear search
            for (size_type idx = 0; idx < sz; ++idx) {
                if (key_compare{}(block[idx], value))   // value >= block[idx]
                    continue;
                if (key_equal{}(value, block[idx]))
                    return {true, idx};
                if (key_compare{}(value, block[idx]))
                    return {false, idx};
            }
            return {false, sz};
        }
    }

    inline void clear() noexcept {
        sz = 0;
    }

private:
    Type block[BlockSize]{};
    size_type sz = 0;

    // only for sorted values, SFINAE, lower bound
    template<class U = Type, class = std::enable_if_t<!std::is_pointer_v<U> > >
    static const Type *binarySearch(const Type *first, const Type *last, const Type &value) {
        const Type *it;
        for (ptrdiff_t count = last - first, step; count > 0;) {
            it = first;
            step = count / 2;
            it += step;
            if (key_compare{}(*it, value)) {
                first = ++it;
                count -= step + 1;
            } else
                count = step;
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

    static_assert(Order >= 1);

    friend class InternalNode<key_type, Order>;

    friend class LeafNode<key_type, Order>;

    const bool isLeaf = true;

    explicit BaseNode(const bool isLeaf_) : isLeaf{isLeaf_} {}

    virtual ~BaseNode() noexcept = default;

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

    [[nodiscard]] base_node *createNewRootNodeFromOld(base_node &newChild) noexcept {
        auto newRootNode = new internal_node;
        newRootNode->children.pushFront(this);
        newRootNode->insertChild(0, newChild);
        return newRootNode;
    }

protected:
    NodeDataBlock<key_type, 2 * Order> keys;
};

template<class KeyType, size_t Order>
class InternalNode : public BaseNode<KeyType, Order> {
public:
    using key_type = KeyType;
    using key_compare = std::less<key_type>;
    using key_equal = std::equal_to<key_type>;

    using size_type = std::size_t;

    using base_node = BaseNode<KeyType, Order>;
    using internal_node = InternalNode<KeyType, Order>;
    using leaf_node = LeafNode<KeyType, Order>;

    friend class BaseNode<KeyType, Order>;

    friend class LeafNode<KeyType, Order>;

    InternalNode() : base_node{false} {}

    ~InternalNode() noexcept override {
        for (const auto child: children)
            delete child;
    }

    friend std::ostream &operator<<(std::ostream &os, const internal_node &node) noexcept {
        node.print(os);
        return os;
    }

    void print(std::ostream &os) const noexcept {
        os << "InternalNode: this = " << this << " keys (size " << this->keys.size() << "): ";
        for (const auto &key: this->keys)
            os << key << ' ';
        os << " children (size " << children.size() << "): ";
        for (const auto child: children)
            os << child << ' ';
    }

    // returns nullptr if was split
    [[nodiscard]] internal_node *addChildToInternalNode(const size_type idx, base_node &newChild) noexcept {
        if (!this->keys.full()) {
            insertChild(idx, newChild);
            return nullptr;
        }
        InternalNode *newRightINode = splitInternalNode(idx, newChild);
        return newRightINode;
    }

    [[nodiscard]] inline const NodeDataBlock<base_node *, 2 * Order + 1> &getChildren() const noexcept {
        return children;
    }

    void fixProblemChild(const size_type problemChildIdx) noexcept {
        assert(!this->keys.empty() && !children.empty());
        if (problemChildIdx > 0 && children[problemChildIdx - 1]->hasEnoughKeysToBorrow())
            borrowKeyFromLeftSibling(problemChildIdx);
        else if (problemChildIdx < children.size() - 1 && children[problemChildIdx + 1]->hasEnoughKeysToBorrow())
            borrowKeyFromRightSibling(problemChildIdx);
        else {
            size_type leftChildIdx = std::min(problemChildIdx, static_cast<size_type>(this->keys.size()) - 1);
            mergeWithRightSibling(leftChildIdx);
        }
    }

    base_node *pullUpChildToRootNode() noexcept {
        assert(this->keys.empty() && children.size() == 1);
        base_node *newRootNode = children.front();
        children.clear();   // nullify children
        return newRootNode;
    }

private:
    NodeDataBlock<base_node *, 2 * Order + 1> children;

    void insertChild(const size_type idx, base_node &newChild) noexcept {
        assert(!this->keys.full() && !children.full() && !children.empty() && !newChild.keys.empty());
        assert(this->keys.contains(newChild.keys.front()).second == idx);
        this->keys.insertAndShiftRight(idx, newChild.keys.front());
        children.insertAndShiftRight(idx + 1, &newChild);
        // deletes the temporary first key of new right internal node child
        if (!newChild.isLeaf)
            newChild.keys.popFront();
    }

    [[nodiscard]] InternalNode *splitInternalNode(const size_type idx, base_node &newChild) noexcept {
        assert(this->keys.full() && children.full());
        auto newRightInternalNode = new internal_node;
        // the first key of new right internal node is saved in order to delete it later in insertChild or splitInternalNode parent method
        newRightInternalNode->keys = this->keys.splitAddNewValue(idx, newChild.keys.front());
        newRightInternalNode->children = children.splitAddNewValue(idx + 1, &newChild);
        // deletes the temporary first key of new right internal node child
        if (!newChild.isLeaf)
            newChild.keys.popFront();
        return newRightInternalNode;
    }

    void borrowKeyFromLeftSibling(const size_type problemChildIdx) noexcept {
        assert(problemChildIdx > 0 && children[problemChildIdx - 1]->hasEnoughKeysToBorrow());
        if (children[problemChildIdx]->isLeaf) {
            auto &leftSibling = static_cast<leaf_node &>(*children[problemChildIdx - 1]);
            auto &rightProblemChild = static_cast<leaf_node &>(*children[problemChildIdx]);
            rightProblemChild.keys.pushFront(leftSibling.keys.back());
            leftSibling.keys.popBack();
            this->keys[problemChildIdx - 1] = rightProblemChild.keys.front();  // update index key
        } else {
            auto &leftSibling = static_cast<internal_node &>(*children[problemChildIdx - 1]);
            auto &rightProblemChild = static_cast<internal_node &>(*children[problemChildIdx]);
            // left sibling last key ->/ left to right problem child parent key \-> right problem child first key
            rightProblemChild.keys.pushFront(this->keys[problemChildIdx - 1]);  // move down index key
            this->keys[problemChildIdx - 1] = leftSibling.keys.back();  // update index key
            rightProblemChild.children.pushFront(leftSibling.children.back());
            leftSibling.keys.popBack();
            leftSibling.children.popBack();
        }
    }

    void borrowKeyFromRightSibling(const size_type problemChildIdx) noexcept {
        assert(problemChildIdx < children.size() - 1 && children[problemChildIdx + 1]->hasEnoughKeysToBorrow());
        if (children[problemChildIdx]->isLeaf) {
            auto &leftProblemChild = static_cast<leaf_node &>(*children[problemChildIdx]);
            auto &rightSibling = static_cast<leaf_node &>(*children[problemChildIdx + 1]);
            leftProblemChild.keys.pushBack(rightSibling.keys.front());
            rightSibling.keys.popFront();
            this->keys[problemChildIdx] = rightSibling.keys.front();  // update index key
        } else {
            auto &leftProblemChild = static_cast<internal_node &>(*children[problemChildIdx]);
            auto &rightSibling = static_cast<internal_node &>(*children[problemChildIdx + 1]);
            // left problem child last key <-/ left to right sibling parent key \<- right sibling first key
            leftProblemChild.keys.pushBack(this->keys[problemChildIdx]);  // move down index key
            this->keys[problemChildIdx] = rightSibling.keys.front();  // update index key
            leftProblemChild.children.pushBack(rightSibling.children.front());
            rightSibling.keys.popFront();
            rightSibling.children.popFront();
        }
    }

    void mergeWithRightSibling(const size_type leftChildIdx) noexcept {
        assert(leftChildIdx < this->keys.size());
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
            leftChild.keys.pushBack(this->keys[leftChildIdx]);
            leftChild.keys.merge(rightChild.keys);
            leftChild.children.merge(rightChild.children);
            rightChild.children.clear();   // nullify children
            delete &rightChild;
        }
        this->keys.eraseAndShiftLeft(leftChildIdx);
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
    using internal_node = InternalNode<KeyType, Order>;
    using leaf_node = LeafNode<KeyType, Order>;

    friend class BaseNode<KeyType, Order>;

    friend class InternalNode<KeyType, Order>;

    LeafNode() : base_node{true} {}

    friend std::ostream &operator<<(std::ostream &os, const leaf_node &node) noexcept {
        node.print(os);
        return os;
    }

    void print(std::ostream &os) const noexcept {
//        os << "LeafNode: prev = " << prev << " this = " << this << " next = " << next << " keys (size "
//           << this->keys.size() << "): ";
        os << "LeafNode: this = " << this << " next = " << next << " keys (size " << this->keys.size() << "): ";
        for (const auto &key: this->keys)
            os << key << ' ';
    }

    // returns nullptr if was split
    [[nodiscard]] leaf_node *addKeyToLeaf(const size_type idx, const key_type &key) noexcept {
        if (!this->keys.full()) {
            this->keys.insertAndShiftRight(idx, key);
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
        this->keys.eraseAndShiftLeft(idx);
    }

private:
    leaf_node *next = nullptr;

    [[nodiscard]] leaf_node *splitLeaf(const size_type idx, const key_type &newKey) noexcept {
        assert(this->keys.full());
        auto newRightLeaf = new leaf_node;
        newRightLeaf->keys = this->keys.splitAddNewValue(idx, newKey);
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
        auto foundStack = leafSearchWithPath(key);
        const auto [leaf, leafIdx] = foundStack.second.top();
        if (foundStack.first) // contains
            return {const_iterator{static_cast<leaf_node *>(leaf), leafIdx}, false};
        tryAddKeyWithPath(key, foundStack);
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
        const auto [found, leaf, leafIdx] = leafSearch(key);
        if (found) // contains
            return {static_cast<leaf_node *>(leaf), leafIdx};
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
        std::pair<bool, Stack<std::pair<base_node *, size_type> > > foundStack = leafSearchWithPath(key);
        return tryAddKeyWithPath(key, foundStack);
    }

    bool tryAddKeyWithPath(const key_type &key,
                           std::pair<bool, Stack<std::pair<base_node *, size_type> > > &foundStack) noexcept {
        auto &[found, stack] = foundStack;
        if (found) // duplicate
            return false;
        assert(!stack.empty());
        const auto [leaf, leafIdx] = stack.top();
        stack.pop();
        ++sz;
        base_node *newRightChild = static_cast<leaf_node *>(leaf)->addKeyToLeaf(leafIdx, key);
        if (newRightChild == nullptr)
            return true;
        while (!stack.empty()) {
            const auto [internalNode, childIdx] = stack.top();
            stack.pop();
            newRightChild = static_cast<internal_node &>(*internalNode).addChildToInternalNode(childIdx,
                                                                                               *newRightChild);
            if (newRightChild == nullptr)
                return true;
        }
        ++height;
        roodNode = roodNode->createNewRootNodeFromOld(*newRightChild);
        return true;
    }

    bool tryRemoveKey(const Key &key) noexcept {
        if (empty())
            return false;
        auto [found, stack] = leafSearchWithPath(key);
        if (!found)
            return false;
        const auto [leaf, leafIdx] = stack.top();
        stack.pop();
        static_cast<leaf_node &>(*leaf).removeKeyFromLeaf(leafIdx);
        --sz;
        if (roodNode->isLeaf || !leaf->needsMerging())
            return true;
        while (!stack.empty()) {
            const auto [internalNode, childIdx] = stack.top();
            stack.pop();
            static_cast<internal_node &>(*internalNode).fixProblemChild(childIdx);
            if (!internalNode->needsMerging())
                return true;
        }
        if (!roodNode->getKeys().empty())
            return true;
        --height;
        base_node *newRootNode = static_cast<internal_node &>(*roodNode).pullUpChildToRootNode();
        delete roodNode;
        roodNode = newRootNode;
        return true;
    }

    // returns path stack with nodes and child indices, on the top is leaf child with key index
    [[nodiscard]] std::pair<bool, Stack<std::pair<base_node *, size_type> > >
    leafSearchWithPath(const key_type &key) const noexcept {
        Stack<std::pair<base_node *, size_type> > stack{height + 1};
        base_node *node = roodNode;
        auto [found, idx] = node->getKeys().contains(key);
        while (!node->isLeaf) {
            // if key is found in internal node, we need to go to the right child,
            // contains returns lower bound <=, so we need to increment child idx
            if (found)
                ++idx;
            stack.emplace(node, idx);
            node = static_cast<internal_node &>(*node).getChildren()[idx];
            std::tie(found, idx) = node->getKeys().contains(key);
        }
        stack.emplace(node, idx);
        return {found, std::move(stack)};
    }

    // returns leaf child with key index
    [[nodiscard]] std::tuple<bool, base_node *, size_type> leafSearch(const key_type &key) const noexcept {
        base_node *node = roodNode;
        auto [found, idx] = node->getKeys().contains(key);
        while (!node->isLeaf) {
            // if key is found in internal node, we need to go to the right child,
            // contains returns lower bound <=, so we need to increment child idx
            if (found)
                ++idx;
            node = static_cast<internal_node &>(*node).getChildren()[idx];
            std::tie(found, idx) = node->getKeys().contains(key);
        }
        return {found, node, idx};
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
