#ifndef UTILITIES_SIMPLE_TREE_HPP
#define UTILITIES_SIMPLE_TREE_HPP

#include <vector>
#include <cassert>
#include <iterator>
#include <type_traits>
#include "range.hpp"

namespace util {

#if 0
/// Defines the order for traversing a tree.
enum TraversalOrder
{
	NativeOrder = 0,
	PreOrder = 1,
	InOrder = 2,
	PostOrder = 3,
	BreadthFirst = 4,
};
#endif

/**
 * Represents a simple tree which only allows adding children to the last
 * node of each depth.
 */
template <class T, class TDepth = size_t>
class simple_tree
{
private:

	//typedef TDepth TSize;

	struct node_t
	{
		T data;     // user data associated with this node
		//TSize size; // 1 + number of children of this node
		TDepth depth;
		//node_t(const T& _data) : data(_data), size(1) { }
		node_t(const T& _data, TDepth _depth) : data(_data), depth(_depth) { }
	};

	std::vector<node_t> _nodes;

	typedef simple_tree<T,TDepth> self_type;

	/// Represents an abstract iterator that points to a specific node 
	/// in a tree, but without the ability to increment itself.
	template <bool IsConst>
	class node_iterator
	{
	protected:

		/// Declares the tree class as friend.
		friend class simple_tree<T,TDepth>;

		/// Type of the tree.
		typedef typename std::conditional<IsConst,
			const self_type, self_type>::type Tree;

		/// Pointer to the underlying tree.
		Tree * _tree;

		/// Index of the node.
		size_t _index;

	public:

		/// Iterator category.
		typedef std::forward_iterator_tag iterator_category;

		/// Type of the data associated with each node.
		typedef T value_type;

		/// Type of iterator difference.
		typedef std::ptrdiff_t difference_type;

		/// Type of reference to data.
		typedef typename std::conditional<IsConst,
			const value_type &, value_type &>::type reference;

		/// Type of pointer to data.
		typedef typename std::conditional<IsConst,
			const value_type *, value_type *>::type pointer;

		/// Constructs an unspecified iterator.
		node_iterator() : _tree(0), _index(0) { }

		/// Copy-constructs an iterator.
		node_iterator(const node_iterator<false> &other)
			: _tree(other._tree), _index(other._index) { }

		/// Constructs an iterator that points to a specific node in a tree.
		node_iterator(Tree *tree, size_t index)
			: _tree(tree), _index(index) { }

		/// Returns the depth of this node (root = 0).
		TDepth depth() const { return _tree->_nodes[_index].depth; }

		/// Returns a reference to the node data.
		reference operator * () const { return _tree->_nodes[_index].data; }

		/// Returns a pointer to the node data.
		pointer operator -> () const { return &_tree->_nodes[_index].data; }

		/// Tests whether the iterator is empty.
		bool operator ! () const { return _tree == 0; }

		/// Tests whether the iterator is not empty.
		operator void* () const { return (_tree == 0)? 0 : (void*)this; }

		/// Tests whether two iterators are equal.
		template <bool IsConst2>
		bool operator == (const node_iterator<IsConst2> &it) const
		{
			return (_index == it._index) && (_tree == it._tree);
		}

		/// Tests whether two iterators are not equal.
		template <bool IsConst2>
		bool operator != (const node_iterator<IsConst2> &it) const
		{
			return ! operator == (it);
		}
	};

	/// Iterator that traverses the siblings of a given node.
	template <bool IsConst>
	class sibling_iterator : public node_iterator<IsConst>
	{
	protected:

		friend class simple_tree<T,TDepth>;

		typedef node_iterator<IsConst> Base;

		/// Constructs an iterator that points to a given node.
		sibling_iterator(typename Base::Tree *tree, size_t index) 
			: node_iterator<IsConst>(tree, index) { }

	public:

		/// Copy- or convert-constructs an iterator.
		template <bool IsConst2>
		sibling_iterator(const node_iterator<IsConst2> &other)
			: node_iterator<IsConst>(other) { }

		/// Advances the iterator to the next sibling of this node. If there
		/// is no next sibling, the advanced iterator points to a position 
		/// where the next sibling would have been inserted.
		sibling_iterator& operator ++ ()
		{
			typename Base::Tree *tree = Base::_tree;
			size_t index = Base::_index;
			TDepth d = tree->_nodes[index].depth;
			size_t i = index + 1;
			while (i < tree->_nodes.size() && tree->_nodes[i].depth > d)
				++i;
			Base::_index = i;
			return *this;
		}
	};

	/// Iterator that traverses a branch of the tree in preorder.
	template <bool IsConst>
	class preorder_iterator : public node_iterator<IsConst>
	{
	protected:

		friend class simple_tree<T,TDepth>;

		typedef node_iterator<IsConst> Base;

		/// Constructs an iterator that points to a given node.
		preorder_iterator(typename Base::Tree *tree, size_t index)
			: node_iterator<IsConst>(tree, index) { }

	public:

		/// Copy- or convert-constructs an iterator.
		preorder_iterator(const node_iterator<false> &other)
			: node_iterator<IsConst>(other) { }

		/// Advances the iterator to the next node in preorder traversal.
		preorder_iterator& operator ++ ()
		{
			++Base::_index;
			return *this;
		}
	};

public:

	/// Type of the data associated with the nodes in the tree.
	typedef T value_type;

	/// Type of the depth of a node.
	typedef TDepth depth_type;

	/// Type of a sibling iterator.
	typedef node_iterator<false> iterator;

	/// Type of a const sibling iterator.
	typedef node_iterator<true> const_iterator;

public:

	/// Creates a tree with the given root data.
	simple_tree(const T& root_data)
	{
		_nodes.push_back(node_t(root_data, 0));
	}

	/// Returns the number of nodes in the tree, including the root node.
	/// The return value is always greater than zero because the root node
	/// always exists.
	size_t size() const { return _nodes.size(); }

	/// Returns a sibling iterator to the root of the tree.
	iterator root() { return iterator(this, 0); }

	/// Returns an iterator to the root of the tree.
	const_iterator root() const { return const_iterator(this, 0); }

	/// Returns a pair of iterators suitable for traversing the children
	/// of a given node.
	range<sibling_iterator<false>> children(const_iterator it)
	{
		return range<sibling_iterator<false>>(
			sibling_iterator<false>(this, it._index + 1),
			++sibling_iterator<false>(this, it._index));
	}

	/// Returns a pair of const iterators suitable for traversing the children
	/// of a given node.
	range<sibling_iterator<true>> children(const_iterator it) const
	{
		return const_cast<self_type*>(this)->children(it);
	}

	/// Returns a pair of iterators suitable for traversing all nodes in 
	/// a given branch in native (preorder) order.
	range<preorder_iterator<false>> traverse(const_iterator root)
	{
		assert(root._tree == this);
		assert(root._index >= 0 && root._index < _nodes.size());
		return range<preorder_iterator<false>>(
			preorder_iterator<false>(this, root._index),
			++sibling_iterator<false>(this, root._index));
	}

	/// Returns a pair of const iterators suitable for traversing all nodes
	/// in a given branch in native (preorder) order.
	range<preorder_iterator<true>> traverse(const_iterator root) const
	{
		return const_cast<self_type*>(this)->traverse(root);
	}

	/**
	 * Inserts a new node as the last child of an existing node.
	 *
	 * @param where Iterator to the parent node.
	 * @param data The node data.
	 * @returns An iterator to the added node.
	 *
	 * @remarks If the child is added to the last sibling of any level, all
	 *      iterators are still valid. Otherwise, only the iterators on the
	 *      parent path of the child is guaranteed to be valid.
	 */
	iterator insert_child(const_iterator where, const T& data)
	{
		assert(where._tree == this);

		// For now, we only support adding to the end of the tree.
		sibling_iterator<true> it(where);
		assert((++it)._index == _nodes.size());

		// Append the child to the end of the tree.
		_nodes.push_back(node_t(data, where.depth() + 1));
		return iterator(this, _nodes.size() - 1);
	}

	/**
	 * Inserts another tree as child of an existing node.
	 *
	 * @param where Iterator to the parent node.
	 * @param subtree A subtree to insert.
	 * @param has_root Flag indicating whether the subtree is rooted. If this
	 *      argument is @c true, the root node of the subtree is inserted as
	 *      a child of the existing node. If this argument is @c false, the
	 *      level-1 nodes of the subtree are inserted as the children.
	 * @returns An iterator to the added node.
	 *
	 * @remarks If the child is added to the last sibling of any level, all
	 *      iterators are still valid. Otherwise, only the iterators on the
	 *      parent path of the child is guaranteed to be valid.
	 */
	iterator insert_child(
		const_iterator where,
		const simple_tree &subtree, 
		bool has_root = true)
	{
		assert(where._tree == this);

		// For now, we only support adding to the end of the tree.
		sibling_iterator<true> it(where);
		assert((++it)._index == _nodes.size());

		iterator ret(this, _nodes.size());

		// Append the subtree to the end of the tree. We need to update the
		// depth field of the tree.
		size_t offset = has_root ? 0 : 1;
		_nodes.reserve(_nodes.size() + subtree._nodes.size() - offset);
		size_t base_depth = where.depth() + 1 - offset;
		size_t n = subtree._nodes.size();
		for (size_t i = offset; i < n; ++i)
		{
			const node_t &child = subtree._nodes[i];
			_nodes.push_back(node_t(child.data, child.depth + (TDepth)base_depth));
		}
		return ret;
	}
};

} // namespace util

#endif // UTILITIES_SIMPLE_TREE_HPP
