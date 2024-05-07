#ifndef ADS_SET_ADS_SET_H
#define ADS_SET_ADS_SET_H

#include <initializer_list>
#include <iostream>

#define NDEBUG

#include <cassert>

template<class Type, size_t BlockSize>
class NodeDataBlock {
public:
    using key_type = Type;
    using key_compare = std::less<key_type>;
    using key_equal = std::equal_to<key_type>;

    class View {
    public:
        View(const Type *const parr_, const size_t sz_) : parr{parr_}, sz{sz_} {}

        [[nodiscard]] inline const Type *begin() const noexcept {
            return parr;
        }

        [[nodiscard]] inline const Type *end() const noexcept {
            return parr + sz;
        }

        [[nodiscard]] inline size_t size() const noexcept {
            return sz;
        }

        [[nodiscard]] inline const Type &operator[](const size_t idx) const noexcept {
            assert(idx < sz);
            return parr[idx];
        }

    private:
        const Type *parr = nullptr;
        size_t sz = 0;
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

    [[nodiscard]] inline const Type &operator[](const size_t idx) const noexcept {
        assert(idx < sz && !empty());
        return block[idx];
    }

    [[nodiscard]] inline Type &operator[](const size_t idx) noexcept {
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

    [[nodiscard]] inline size_t size() const noexcept {
        return sz;
    }

    [[nodiscard]] inline bool empty() const noexcept {
        return sz == 0;
    }

    [[nodiscard]] inline bool full() const noexcept {
        return sz >= BlockSize;
    }

    void insertAndShiftRight(const size_t idx, const Type &value) noexcept {
        assert(idx < BlockSize && !full());
        for (size_t i = sz; i > idx; --i)
            block[i] = std::move(block[i - 1]); // no overflow because now i >= 1
        block[idx] = value;
        ++sz;
    }

    inline void pushFront(const Type &value) noexcept {
        insertAndShiftRight(0, value);
    }

    inline void pushBack(const Type &value) noexcept {
        insertAndShiftRight(sz, value);
    }

    void eraseAndShiftLeft(const size_t idx) noexcept {
        assert(idx < BlockSize && !empty());
        for (size_t i = idx; i < sz - 1; ++i)   // no overflow because now sz >= 1
            block[i] = std::move(block[i + 1]);
        --sz;
    }

    inline void popFront() noexcept {
        eraseAndShiftLeft(0);
    }

    inline void popBack() noexcept {
        eraseAndShiftLeft(sz - 1);
    }

    NodeDataBlock splitAddNewValue(const size_t newValueIdx, const Type &newValue) noexcept {
        assert(newValueIdx <= BlockSize && full());
        // oldBlockSize >= newBlock_size
        size_t oldBlockSize = (BlockSize + 1) / 2;
        NodeDataBlock newBlock;
        if (newValueIdx < oldBlockSize) {  // insert value in old block
            --oldBlockSize;
            newBlock.sz = sz - oldBlockSize;
            sz = oldBlockSize;
            for (size_t i = 0; i < newBlock.sz; ++i)
                newBlock.block[i] = std::move(block[i + oldBlockSize]);
            insertAndShiftRight(newValueIdx, newValue);
        } else {  // insert value in new block
            newBlock.sz = sz - oldBlockSize + 1;
            sz = oldBlockSize;
            for (size_t i = 0; i < newValueIdx - oldBlockSize; ++i)
                newBlock.block[i] = std::move(block[i + oldBlockSize]);
            newBlock.block[newValueIdx - oldBlockSize] = newValue;
            for (size_t i = newValueIdx - oldBlockSize + 1; i < newBlock.sz; ++i)
                newBlock.block[i] = std::move(block[i + oldBlockSize - 1]);
        }
        return newBlock;
    }

    void merge(View rightBlock) noexcept {
        assert(size() + rightBlock.size() <= BlockSize);
        for (size_t i = 0; i < rightBlock.size(); ++i)
            block[i + sz] = std::move(rightBlock[i]);
        sz += rightBlock.size();
    }

    // only for sorted values, SFINAE
    // returns (true, idx) if contains value, else (false, idx) of first element, so that value <= element
    template<class U = Type, class = std::enable_if_t<!std::is_pointer_v<U> > >
    [[nodiscard]] std::pair<bool, size_t> contains(const Type &value) const noexcept {
        const Type *it = binarySearch(begin(), end(), value);  // lower bound, value <= *it
        const size_t idx = it - begin();
        if (it != end() && key_equal{}(value, *it))  // contains value
            return {true, idx};
        return {false, idx};
    }

private:
    Type block[BlockSize]{};
    size_t sz = 0;

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
    explicit Stack(const size_t max_sz_) noexcept: data{new Type[max_sz_]}, max_sz{max_sz_}, sz{0} {}

    ~Stack() noexcept {
        delete[] data;
    }

    [[nodiscard]] inline size_t size() const noexcept {
        return sz;
    }

    [[nodiscard]] inline bool empty() const noexcept {
        return sz == 0;
    }

    [[nodiscard]] inline size_t capacity() const noexcept {
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
    size_t max_sz = 0;
    size_t sz = 0;
};

template<class KeyType, size_t Order>
class Node {
public:
    using key_type = KeyType;
    using key_compare = std::less<key_type>;
    using key_equal = std::equal_to<key_type>;

    static_assert(Order >= 1);
    static constexpr size_t minNumberOfKeys = Order;
    static constexpr size_t maxNumberOfKeys = 2 * Order;

    enum class NodeType {
        LeafNode,
        InternalNode
    };

    explicit Node(NodeType type_) noexcept: type{type_} {}

    ~Node() {
        if (!isLeaf())
            for (const auto child: children)
                delete child;
    }

    friend std::ostream &operator<<(std::ostream &os, const Node &node) noexcept {
        node.print(os);
        return os;
    }

    void print(std::ostream &os) const noexcept {
        switch (type) {
            case NodeType::LeafNode:
                os << "LeafNode: prev = " << prev << " this = " << this << " next = " << next << " keys (size "
                   << keys.size() << "): ";
                for (const auto &key: keys)
                    os << key << ' ';
                return;
            case NodeType::InternalNode:
                os << "InternalNode: this = " << this << " keys (size " << keys.size() << "): ";
                for (const auto &key: keys)
                    os << key << ' ';
                os << " children (size " << children.size() << "): ";
                for (const auto child: children)
                    os << child << ' ';
        }
    }

    // returns nullptr if was split
    [[nodiscard]] Node *addKeyToLeaf(const size_t idx, const KeyType &key) noexcept {
        assert(isLeaf());
        if (!keys.full()) {
            keys.insertAndShiftRight(idx, key);
            return nullptr;
        }
        Node *newRightLeaf = splitLeaf(idx, key);
        return newRightLeaf;
    }

    // returns nullptr if was split
    [[nodiscard]] Node *addChildToInternalNode(const size_t idx, Node &newChild) noexcept {
        assert(!isLeaf());
        if (!keys.full()) {
            insertChild(idx, newChild);
            return nullptr;
        }
        Node *newRightINode = splitInternalNode(idx, newChild);
        return newRightINode;
    }

    [[nodiscard]] inline bool isLeaf() const noexcept {
        return type == NodeType::LeafNode;
    }

    [[nodiscard]] inline bool needsMerging() const noexcept {
        return keys.size() < minNumberOfKeys;
    }

    [[nodiscard]] inline bool hasEnoughKeysToBorrow() const noexcept {
        assert(keys.size() >= minNumberOfKeys - 1);
        return keys.size() > minNumberOfKeys;
    }

    [[nodiscard]] inline const NodeDataBlock<KeyType, maxNumberOfKeys> &getKeys() const noexcept {
        return keys;
    }

    [[nodiscard]] inline const NodeDataBlock<Node *, maxNumberOfKeys + 1> &getChildren() const noexcept {
        return children;
    }

    [[nodiscard]] inline bool containsKey(const size_t idx, const KeyType &key) const noexcept {
        assert(idx <= keys.size());
        return idx < keys.size() && key_equal{}(keys[idx], key);
    }

    [[nodiscard]] inline const Node *getNext() const noexcept {
        return next;
    }

    [[nodiscard]] inline const Node *getPrev() const noexcept {
        return prev;
    }

    [[nodiscard]] Node *createNewRootNodeFromOld(Node &newChild) noexcept {
        auto newRootNode = new Node<KeyType, Order>{NodeType::InternalNode};
        newRootNode->children.pushFront(this);
        newRootNode->insertChild(0, newChild);
        return newRootNode;
    }

    bool removeKeyFromLeaf(const size_t idx, const KeyType &key) {
        assert(isLeaf());
        if (idx == keys.size() || key_compare{}(key, keys[idx]))  // not found
            return false;
        keys.eraseAndShiftLeft(idx);
        return true;
    }

    void fixProblemChild(const size_t problemChildIdx) noexcept {
        assert(!keys.empty() && !children.empty());
        if (problemChildIdx > 0 && children[problemChildIdx - 1]->hasEnoughKeysToBorrow())
            borrowKeyFromLeftSibling(problemChildIdx);
        else if (problemChildIdx < children.size() - 1 && children[problemChildIdx + 1]->hasEnoughKeysToBorrow())
            borrowKeyFromRightSibling(problemChildIdx);
        else {
            size_t leftChildIdx = std::min(problemChildIdx, keys.size() - 1);
            mergeWithRightSibling(leftChildIdx);
        }
    }

    Node *pullUpChildToRootNode() noexcept {
        assert(keys.empty() && children.size() == 1);
        Node *newRootNode = children.front();
        type = NodeType::LeafNode;  // nullify children
        return newRootNode;
    }

    // returns path stack with nodes and child indices, on the top is leaf child with key index
    [[nodiscard]] Stack<std::pair<Node<KeyType, Order> *, size_t> >
    leafSearchWithPath(const KeyType &key, const size_t height) const noexcept {
        Stack<std::pair<Node<KeyType, Order> *, size_t> > stack{height + 1};
        auto *node = const_cast<Node<KeyType, Order> *>(this);
        auto [found, idx] = node->getKeys().contains(key);
        while (!node->isLeaf()) {
            // if key is found in internal node, we need to go to the right child,
            // contains returns lower bound <=, so we need to increment child idx
            if (found)
                ++idx;
            stack.emplace(node, idx);
            node = node->getChildren()[idx];
            std::tie(found, idx) = node->getKeys().contains(key);
        }
        stack.emplace(node, idx);
        return stack;
    }

    // returns leaf child with key index
    [[nodiscard]] std::pair<Node<KeyType, Order> *, size_t> leafSearch(const KeyType &key) const noexcept {
        auto *node = const_cast<Node<KeyType, Order> *>(this);
        auto [found, idx] = node->getKeys().contains(key);
        while (!node->isLeaf()) {
            // if key is found in internal node, we need to go to the right child,
            // contains returns lower bound <=, so we need to increment child idx
            if (found)
                ++idx;
            node = node->getChildren()[idx];
            std::tie(found, idx) = node->getKeys().contains(key);
        }
        return {node, idx};
    }

    [[nodiscard]] Node *getFirstLeaf() const noexcept {
        auto *node = const_cast<Node<KeyType, Order> *>(this);
        while (!node->isLeaf()) {
            assert(!node->getChildren().empty());
            node = node->children.front();
        }
        return node;
    }

private:
    NodeType type;
    NodeDataBlock<KeyType, maxNumberOfKeys> keys;
    NodeDataBlock<Node *, maxNumberOfKeys + 1> children;
    Node *prev = nullptr;
    Node *next = nullptr;

    [[nodiscard]] Node *splitLeaf(const size_t idx, const KeyType &newKey) noexcept {
        assert(keys.full());
        auto newRightLeaf = new Node<KeyType, Order>(NodeType::LeafNode);
        newRightLeaf->keys = keys.splitAddNewValue(idx, newKey);
        if (next != nullptr)
            next->prev = newRightLeaf;
        newRightLeaf->next = next;
        next = newRightLeaf;
        newRightLeaf->prev = this;
        return newRightLeaf;
    }

    void insertChild(const size_t idx, Node &newChild) noexcept {
        assert(!keys.full() && !children.full() && !children.empty() && !newChild.keys.empty());
        assert(keys.contains(newChild.keys.front()).second == idx);
        keys.insertAndShiftRight(idx, newChild.keys.front());
        children.insertAndShiftRight(idx + 1, &newChild);
        // deletes the temporary first key of new right internal node child
        if (!newChild.isLeaf())
            newChild.keys.popFront();
    }

    [[nodiscard]] Node *splitInternalNode(const size_t idx, Node &newChild) noexcept {
        assert(keys.full() && children.full());
        auto newRightInternalNode = new Node<KeyType, Order>{NodeType::InternalNode};
        // the first key of new right internal node is saved in order to delete it later in insertChild or splitInternalNode parent method
        newRightInternalNode->keys = keys.splitAddNewValue(idx, newChild.keys.front());
        newRightInternalNode->children = children.splitAddNewValue(idx + 1, &newChild);
        // deletes the temporary first key of new right internal node child
        if (!newChild.isLeaf())
            newChild.keys.popFront();
        return newRightInternalNode;
    }

    void borrowKeyFromLeftSibling(const size_t problemChildIdx) noexcept {
        assert(problemChildIdx > 0 && children[problemChildIdx - 1]->hasEnoughKeysToBorrow());
        Node &leftSibling = *children[problemChildIdx - 1];
        Node &rightProblemChild = *children[problemChildIdx];
        switch (rightProblemChild.type) {
            case NodeType::LeafNode:
                rightProblemChild.keys.pushFront(leftSibling.keys.back());
                leftSibling.keys.popBack();
                keys[problemChildIdx - 1] = rightProblemChild.keys.front();  // update index key
                return;
            case NodeType::InternalNode:  // left sibling last key ->/ left to right problem child parent key \-> right problem child first key
                rightProblemChild.keys.pushFront(keys[problemChildIdx - 1]);  // move down index key
                keys[problemChildIdx - 1] = leftSibling.keys.back();  // update index key
                rightProblemChild.children.pushFront(leftSibling.children.back());
                leftSibling.keys.popBack();
                leftSibling.children.popBack();
        }
    }

    void borrowKeyFromRightSibling(const size_t problemChildIdx) noexcept {
        assert(problemChildIdx < children.size() - 1 && children[problemChildIdx + 1]->hasEnoughKeysToBorrow());
        Node &leftProblemChild = *children[problemChildIdx];
        Node &rightSibling = *children[problemChildIdx + 1];
        switch (leftProblemChild.type) {
            case NodeType::LeafNode:
                leftProblemChild.keys.pushBack(rightSibling.keys.front());
                rightSibling.keys.popFront();
                keys[problemChildIdx] = rightSibling.keys.front();  // update index key
                return;
            case NodeType::InternalNode:  // left problem child last key <-/ left to right sibling parent key \<- right sibling first key
                leftProblemChild.keys.pushBack(keys[problemChildIdx]);  // move down index key
                keys[problemChildIdx] = rightSibling.keys.front();  // update index key
                leftProblemChild.children.pushBack(rightSibling.children.front());
                rightSibling.keys.popFront();
                rightSibling.children.popFront();
        }
    }

    void mergeWithRightSibling(const size_t leftChildIdx) noexcept {
        assert(leftChildIdx < keys.size());
        Node &leftChild = *children[leftChildIdx];
        Node &rightChild = *children[leftChildIdx + 1];
        assert(!leftChild.hasEnoughKeysToBorrow() && !rightChild.hasEnoughKeysToBorrow());
        // left child keys + left to right sibling parent key (than deleted in parent) + right sibling keys
        if (!leftChild.isLeaf())
            leftChild.keys.pushBack(keys[leftChildIdx]);
        leftChild.keys.merge(rightChild.keys);
        leftChild.children.merge(rightChild.children);
        if (rightChild.next != nullptr)
            rightChild.next->prev = &leftChild;
        leftChild.next = rightChild.next;
        rightChild.type = NodeType::LeafNode;  // nullify children
        delete &rightChild;
        keys.eraseAndShiftLeft(leftChildIdx);
        children.eraseAndShiftLeft(leftChildIdx + 1);
    }
};

template<class Key, size_t N = 2>
class Iterator {
public:
    using value_type = Key;
    using reference = const value_type &;
    using pointer = const value_type *;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::bidirectional_iterator_tag;

    Iterator() = default;

    Iterator(const Node<Key, N> *const currentLeaf_, const size_t currentKeyIdx_) noexcept:
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

    Iterator &operator--() noexcept {
        if (currentKeyIdx == 0) {
            currentLeaf = currentLeaf->getPrev();
            currentKeyIdx = currentLeaf->getKeys().size() - 1;
        } else
            --currentKeyIdx;
        return *this;
    }

    Iterator operator--(int) noexcept { // NOLINT(*-dcl21-cpp)
        Iterator old = *this;
        operator--();
        return old;
    }

    friend inline bool operator==(const Iterator &lhs, const Iterator &rhs) noexcept {
        return lhs.currentLeaf == rhs.currentLeaf && lhs.currentKeyIdx == rhs.currentKeyIdx;
    }

    friend inline bool operator!=(const Iterator &lhs, const Iterator &rhs) noexcept {
        return lhs.currentLeaf != rhs.currentLeaf || lhs.currentKeyIdx != rhs.currentKeyIdx;
    }

private:
    const Node<Key, N> *currentLeaf = nullptr;
    size_t currentKeyIdx = 0;
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
    using const_iterator = Iterator<Key, N>;
    using iterator = const_iterator;
    using key_compare = std::less<key_type>;   // B+-Tree
    using key_equal = std::equal_to<key_type>; // Hashing
    using hasher = std::hash<key_type>;        // Hashing

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
        Stack<std::pair<Node<key_type, N> *, size_type> > stack = roodNode->leafSearchWithPath(key, height);
        const_iterator it = findWithLeaf(key, stack.top());
        if (it != end())
            return {it, false};
        tryAddKeyWithPath(key, stack);
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
        roodNode = new Node<key_type, N>{Node<key_type, N>::NodeType::LeafNode};
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
        return findWithLeaf(key, roodNode->leafSearch(key));
    }

    [[nodiscard]] inline const_iterator begin() const noexcept {
        if (empty())
            return end();
        return {roodNode->getFirstLeaf(), 0};
    }

    [[nodiscard]] inline const_iterator end() const noexcept {
        return {nullptr, 0};
    }

    void dump(std::ostream &o = std::cerr) const noexcept {
        Stack<std::pair<Node<Key, N> *, size_type> > stack{2 * N * height + 1};
        stack.emplace(roodNode, 0);
        o << "B+ Tree: size = " << sz << " height = " << height << '\n';
        while (!stack.empty()) {
            const auto [node, level] = stack.top();
            stack.pop();
            for (size_type i = 0; i < level; ++i)
                o << '\t';
            if (node->isLeaf())
                o << *node << '\n';
            else {
                o << *node << '\n';
                for (auto child: node->getChildren())
                    stack.emplace(child, level + 1);  // output right to left, "feature"
            }
        }
    }

    // O(size)
    [[nodiscard]] friend bool operator==(const ADS_set<Key, N> &lhs, const ADS_set<Key, N> &rhs) noexcept {
        if (lhs.size() != rhs.size())
            return false;
        auto endIt = lhs.end();
        for (auto lhsIt = lhs.begin(), rhsIt = rhs.begin(); lhsIt != endIt; ++lhsIt, ++rhsIt)
            if (!key_equal{}(*lhsIt, *rhsIt))
                return false;
        return true;
    }

    // O(size)
    [[nodiscard]] friend inline bool operator!=(const ADS_set<Key, N> &lhs, const ADS_set<Key, N> &rhs) noexcept {
        return !operator==(lhs, rhs);
    }

private:
    Node<key_type, N> *roodNode = new Node<key_type, N>{Node<key_type, N>::NodeType::LeafNode};
    size_type sz = 0;
    size_type height = 0;

    [[nodiscard]] inline const_iterator
    findWithLeaf(const key_type &key, const std::pair<Node<key_type, N> *, size_type> &leafAndLeafIdx) const noexcept {
        const auto [leaf, leafIdx] = leafAndLeafIdx;
        if (leaf->containsKey(leafIdx, key)) // contains
            return {leaf, leafIdx};
        return end();
    }

    inline bool tryAddKey(const key_type &key) noexcept {
        Stack<std::pair<Node<key_type, N> *, size_type> > stack = roodNode->leafSearchWithPath(key, height);
        return tryAddKeyWithPath(key, stack);
    }

    bool tryAddKeyWithPath(const key_type &key, Stack<std::pair<Node<key_type, N> *, size_type> > &stack) noexcept {
        assert(!stack.empty());
        const auto [leaf, leafIdx] = stack.top();
        stack.pop();
        if (leaf->containsKey(leafIdx, key)) // duplicate
            return false;
        ++sz;
        Node<key_type, N> *newRightChild = leaf->addKeyToLeaf(leafIdx, key);
        if (newRightChild == nullptr)
            return true;
        while (!stack.empty()) {
            const auto [internalNode, childIdx] = stack.top();
            stack.pop();
            newRightChild = internalNode->addChildToInternalNode(childIdx, *newRightChild);
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
        Stack<std::pair<Node<key_type, N> *, size_t> > stack = roodNode->leafSearchWithPath(key, height);
        const auto [leaf, leafIdx] = stack.top();
        stack.pop();
        bool removed = leaf->removeKeyFromLeaf(leafIdx, key);
        if (!removed)
            return false;
        --sz;
        if (roodNode->isLeaf() || !leaf->needsMerging())
            return true;
        while (!stack.empty()) {
            const auto [internalNode, childIdx] = stack.top();
            stack.pop();
            internalNode->fixProblemChild(childIdx);
            if (!internalNode->needsMerging())
                return true;
        }
        if (!roodNode->getKeys().empty())
            return true;
        --height;
        Node<key_type, N> *newRootNode = roodNode->pullUpChildToRootNode();
        delete roodNode;
        roodNode = newRootNode;
        return true;
    }
};

// O(1)
template<typename Key, size_t N>
inline void swap(ADS_set<Key, N> &lhs, ADS_set<Key, N> &rhs) noexcept {
    lhs.swap(rhs);
}

#endif
