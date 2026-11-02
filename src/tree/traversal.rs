//! Lazy, non-recursive traversals in child insertion order.
//!
//! ID and node iterators yield `Result<_, TreeError>`. Invalid node references (or
//! more than two children for inorder) yield one error, then exhaust the iterator.
//! Tree structure is assumed to be acyclic, as for the other tree operations.
use std::collections::VecDeque;
use std::iter::FusedIterator;

use super::{Node, NodeId, Tree, TreeError};

/// Lazy root-before-children traversal. See [`Tree::iter_preorder`].
pub struct PreOrderTreeTraversal<'a> {
    tree: &'a Tree,
    pending: Vec<NodeId>,
}

impl Iterator for PreOrderTreeTraversal<'_> {
    type Item = Result<NodeId, TreeError>;

    fn next(&mut self) -> Option<Self::Item> {
        let id = self.pending.pop()?;
        match self.tree.get(id) {
            Ok(node) => {
                self.pending.extend(node.children.iter().rev().copied());
                Some(Ok(id))
            }
            Err(error) => {
                self.pending.clear();
                Some(Err(error))
            }
        }
    }
}

impl FusedIterator for PreOrderTreeTraversal<'_> {}

/// Lazy children-before-root traversal. See [`Tree::iter_postorder`].
pub struct PostOrderTreeTraversal<'a> {
    tree: &'a Tree,
    // Expanded frames yield their node after all children have been visited.
    pending: Vec<(NodeId, bool)>,
}

impl Iterator for PostOrderTreeTraversal<'_> {
    type Item = Result<NodeId, TreeError>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((id, expanded)) = self.pending.pop() {
            if expanded {
                return Some(Ok(id));
            }
            match self.tree.get(id) {
                Ok(node) => {
                    self.pending.push((id, true));
                    self.pending
                        .extend(node.children.iter().rev().map(|&child| (child, false)));
                }
                Err(error) => {
                    self.pending.clear();
                    return Some(Err(error));
                }
            }
        }
        None
    }
}

impl FusedIterator for PostOrderTreeTraversal<'_> {}

/// Lazy left-child, root, right-child traversal. See [`Tree::iter_inorder`].
pub struct InOrderTreeTraversal<'a> {
    tree: &'a Tree,
    pending: Vec<(NodeId, bool)>,
}

impl Iterator for InOrderTreeTraversal<'_> {
    type Item = Result<NodeId, TreeError>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((id, expanded)) = self.pending.pop() {
            if expanded {
                return Some(Ok(id));
            }
            let node = match self.tree.get(id) {
                Ok(node) => node,
                Err(error) => {
                    self.pending.clear();
                    return Some(Err(error));
                }
            };
            if node.children.len() > 2 {
                self.pending.clear();
                return Some(Err(TreeError::IsNotBinary));
            }
            if let Some(&right) = node.children.get(1) {
                self.pending.push((right, false));
            }
            self.pending.push((id, true));
            if let Some(&left) = node.children.first() {
                self.pending.push((left, false));
            }
        }
        None
    }
}

impl FusedIterator for InOrderTreeTraversal<'_> {}

/// Lazy breadth-first traversal. See [`Tree::iter_levelorder`].
pub struct LevelOrderTreeTraversal<'a> {
    tree: &'a Tree,
    pending: VecDeque<NodeId>,
}

impl Iterator for LevelOrderTreeTraversal<'_> {
    type Item = Result<NodeId, TreeError>;

    fn next(&mut self) -> Option<Self::Item> {
        let id = self.pending.pop_front()?;
        match self.tree.get(id) {
            Ok(node) => {
                self.pending.extend(node.children.iter().copied());
                Some(Ok(id))
            }
            Err(error) => {
                self.pending.clear();
                Some(Err(error))
            }
        }
    }
}

impl FusedIterator for LevelOrderTreeTraversal<'_> {}

impl Tree {
    /// Iterate from the tree root in preorder, preserving child order.
    ///
    /// Construction validates the root; iteration reports invalid descendants.
    /// After yielding an error, the iterator is exhausted.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let ids = tree.iter_preorder()?.collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(ids, tree.preorder(tree.root()?)?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_preorder(&self) -> Result<PreOrderTreeTraversal<'_>, TreeError> {
        self.iter_preorder_from(self.root()?)
    }

    /// Iterate over just the subtree at `root` in preorder.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("named subtree").id;
    /// let ids = tree.iter_preorder_from(root)?.collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(ids.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_preorder_from(&self, root: NodeId) -> Result<PreOrderTreeTraversal<'_>, TreeError> {
        self.get(root)?;
        Ok(PreOrderTreeTraversal {
            tree: self,
            pending: vec![root],
        })
    }

    /// Iterate from the tree root in postorder, preserving child order.
    ///
    /// Construction validates the root; iteration reports invalid descendants.
    /// After yielding an error, the iterator is exhausted.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let ids = tree.iter_postorder()?.collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(ids, tree.postorder(tree.root()?)?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_postorder(&self) -> Result<PostOrderTreeTraversal<'_>, TreeError> {
        self.iter_postorder_from(self.root()?)
    }

    /// Iterate over just the subtree at `root` in postorder.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("named subtree").id;
    /// let ids = tree.iter_postorder_from(root)?.collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(ids.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_postorder_from(
        &self,
        root: NodeId,
    ) -> Result<PostOrderTreeTraversal<'_>, TreeError> {
        self.get(root)?;
        Ok(PostOrderTreeTraversal {
            tree: self,
            pending: vec![(root, false)],
        })
    }

    /// Iterate from the tree root in inorder, preserving child order.
    ///
    /// Construction validates the root; iteration reports invalid descendants.
    /// After yielding an error, the iterator is exhausted. More than two children yields `TreeError::IsNotBinary`.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let ids = tree.iter_inorder()?.collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(ids, tree.inorder(tree.root()?)?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_inorder(&self) -> Result<InOrderTreeTraversal<'_>, TreeError> {
        self.iter_inorder_from(self.root()?)
    }

    /// Iterate over just the subtree at `root` in inorder.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("named subtree").id;
    /// let ids = tree.iter_inorder_from(root)?.collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(ids.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_inorder_from(&self, root: NodeId) -> Result<InOrderTreeTraversal<'_>, TreeError> {
        self.get(root)?;
        Ok(InOrderTreeTraversal {
            tree: self,
            pending: vec![(root, false)],
        })
    }

    /// Iterate from the tree root in levelorder, preserving child order.
    ///
    /// Construction validates the root; iteration reports invalid descendants.
    /// After yielding an error, the iterator is exhausted.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let ids = tree.iter_levelorder()?.collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(ids, tree.levelorder(tree.root()?)?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_levelorder(&self) -> Result<LevelOrderTreeTraversal<'_>, TreeError> {
        self.iter_levelorder_from(self.root()?)
    }

    /// Iterate over just the subtree at `root` in levelorder.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("named subtree").id;
    /// let ids = tree.iter_levelorder_from(root)?.collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(ids.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_levelorder_from(
        &self,
        root: NodeId,
    ) -> Result<LevelOrderTreeTraversal<'_>, TreeError> {
        self.get(root)?;
        Ok(LevelOrderTreeTraversal {
            tree: self,
            pending: VecDeque::from([root]),
        })
    }
}

/// Lazy borrowed nodes in traversal order. See [`Tree::iter_preorder_nodes`].
///
/// Errors are yielded once, followed by exhaustion, as with the ID traversals.
pub struct NodeTraversal<'a, I> {
    tree: &'a Tree,
    ids: I,
    failed: bool,
}

impl<'a, I: Iterator<Item = Result<NodeId, TreeError>>> Iterator for NodeTraversal<'a, I> {
    type Item = Result<&'a Node, TreeError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        let result = self.ids.next()?.and_then(|id| self.tree.get(id));
        self.failed = result.is_err();
        Some(result)
    }
}

impl<I: FusedIterator<Item = Result<NodeId, TreeError>>> FusedIterator for NodeTraversal<'_, I> {}

#[derive(Clone, Copy)]
enum Order {
    Pre,
    Post,
    In,
    Level,
}

enum NodeSlot<'a> {
    Available(&'a mut Node),
    Missing,
    Visited,
}

enum MutableVisit<'a> {
    Enter(NodeId),
    Yield(&'a mut Node),
}

/// Lazy mutable nodes in traversal order. See [`Tree::iter_preorder_nodes_mut`].
///
/// Setup takes O(tree storage size) time and memory to split storage into disjoint
/// references; it does not walk or precompute the traversal. Each reachable node
/// is then visited once. References may coexist, including when collected.
/// Children are captured before a node is yielded, so later edits to that node's
/// links do not redirect iteration. Preserve IDs and links to keep a valid tree.
/// Cycles/shared descendants yield [`TreeError::RepeatedTraversalNode`] once and
/// terminate iteration. All tree caches are invalidated before the first yield.
pub struct NodeTraversalMut<'a> {
    slots: Vec<NodeSlot<'a>>,
    pending: VecDeque<MutableVisit<'a>>,
    order: Order,
}

impl<'a> NodeTraversalMut<'a> {
    fn new(tree: &'a mut Tree, root: NodeId, order: Order) -> Self {
        let slots = tree
            .traversal_nodes_mut()
            .iter_mut()
            .map(|node| {
                if node.deleted {
                    NodeSlot::Missing
                } else {
                    NodeSlot::Available(node)
                }
            })
            .collect();
        Self {
            slots,
            pending: VecDeque::from([MutableVisit::Enter(root)]),
            order,
        }
    }

    fn fail(&mut self, error: TreeError) -> Option<Result<&'a mut Node, TreeError>> {
        self.pending.clear();
        self.slots.clear();
        Some(Err(error))
    }
}

impl<'a> Iterator for NodeTraversalMut<'a> {
    type Item = Result<&'a mut Node, TreeError>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(action) = self.pending.pop_front() {
            let id = match action {
                MutableVisit::Yield(node) => return Some(Ok(node)),
                MutableVisit::Enter(id) => id,
            };
            let Some(slot) = self.slots.get_mut(id.index()) else {
                return self.fail(TreeError::NodeNotFound(id));
            };
            // Taking a slot makes it impossible to yield the same &mut twice.
            let node = match std::mem::replace(slot, NodeSlot::Visited) {
                NodeSlot::Available(node) => node,
                NodeSlot::Missing => return self.fail(TreeError::NodeNotFound(id)),
                NodeSlot::Visited => return self.fail(TreeError::RepeatedTraversalNode(id)),
            };
            match self.order {
                Order::Pre => {
                    push_children(&mut self.pending, &node.children);
                    return Some(Ok(node));
                }
                Order::Level => {
                    self.pending
                        .extend(node.children.iter().copied().map(MutableVisit::Enter));
                    return Some(Ok(node));
                }
                Order::Post => {
                    let children = node.children.clone();
                    self.pending.push_front(MutableVisit::Yield(node));
                    push_children(&mut self.pending, &children);
                }
                Order::In => {
                    if node.children.len() > 2 {
                        return self.fail(TreeError::IsNotBinary);
                    }
                    let left = node.children.first().copied();
                    let right = node.children.get(1).copied();
                    if let Some(right) = right {
                        self.pending.push_front(MutableVisit::Enter(right));
                    }
                    self.pending.push_front(MutableVisit::Yield(node));
                    if let Some(left) = left {
                        self.pending.push_front(MutableVisit::Enter(left));
                    }
                }
            }
        }
        None
    }
}

// Reverse insertion preserves child order for both depth-first traversals.
fn push_children(pending: &mut VecDeque<MutableVisit<'_>>, children: &[NodeId]) {
    for &child in children.iter().rev() {
        pending.push_front(MutableVisit::Enter(child));
    }
}

impl FusedIterator for NodeTraversalMut<'_> {}

impl Tree {
    /// Borrow nodes of the tree in preorder; errors terminate iteration.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let names = tree.iter_preorder_nodes()?
    ///     .map(|node| node.map(|n| n.name.clone()))
    ///     .collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(names.len(), 5);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_preorder_nodes(
        &self,
    ) -> Result<NodeTraversal<'_, PreOrderTreeTraversal<'_>>, TreeError> {
        Ok(NodeTraversal {
            tree: self,
            ids: self.iter_preorder()?,
            failed: false,
        })
    }

    /// Borrow nodes of the subtree at `root` in preorder; errors terminate iteration.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("subtree").id;
    /// let names = tree.iter_preorder_nodes_from(root)?
    ///     .map(|node| node.map(|n| n.name.clone()))
    ///     .collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(names.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_preorder_nodes_from(
        &self,
        root: NodeId,
    ) -> Result<NodeTraversal<'_, PreOrderTreeTraversal<'_>>, TreeError> {
        Ok(NodeTraversal {
            tree: self,
            ids: self.iter_preorder_from(root)?,
            failed: false,
        })
    }

    /// Mutably borrow nodes of the tree in preorder.
    ///
    /// Setup is O(tree storage size); traversal is lazy. See [`NodeTraversalMut`]
    /// for mutation, cache invalidation and error semantics.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// tree.iter_preorder_nodes_mut()?.try_for_each(|node| {
    ///     node?.comment = Some("visited".into());
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_preorder_nodes_mut(&mut self) -> Result<NodeTraversalMut<'_>, TreeError> {
        self.iter_preorder_nodes_mut_from(self.root()?)
    }

    /// Mutably borrow nodes of the subtree at `root` in preorder.
    ///
    /// Setup is O(tree storage size); traversal is lazy. See [`NodeTraversalMut`]
    /// for mutation, cache invalidation and error semantics.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("subtree").id;
    /// tree.iter_preorder_nodes_mut_from(root)?.try_for_each(|node| {
    ///     node?.comment = Some("visited".into());
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_preorder_nodes_mut_from(
        &mut self,
        root: NodeId,
    ) -> Result<NodeTraversalMut<'_>, TreeError> {
        self.get(root)?;
        Ok(NodeTraversalMut::new(self, root, Order::Pre))
    }

    /// Borrow nodes of the tree in postorder; errors terminate iteration.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let names = tree.iter_postorder_nodes()?
    ///     .map(|node| node.map(|n| n.name.clone()))
    ///     .collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(names.len(), 5);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_postorder_nodes(
        &self,
    ) -> Result<NodeTraversal<'_, PostOrderTreeTraversal<'_>>, TreeError> {
        Ok(NodeTraversal {
            tree: self,
            ids: self.iter_postorder()?,
            failed: false,
        })
    }

    /// Borrow nodes of the subtree at `root` in postorder; errors terminate iteration.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("subtree").id;
    /// let names = tree.iter_postorder_nodes_from(root)?
    ///     .map(|node| node.map(|n| n.name.clone()))
    ///     .collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(names.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_postorder_nodes_from(
        &self,
        root: NodeId,
    ) -> Result<NodeTraversal<'_, PostOrderTreeTraversal<'_>>, TreeError> {
        Ok(NodeTraversal {
            tree: self,
            ids: self.iter_postorder_from(root)?,
            failed: false,
        })
    }

    /// Mutably borrow nodes of the tree in postorder.
    ///
    /// Setup is O(tree storage size); traversal is lazy. See [`NodeTraversalMut`]
    /// for mutation, cache invalidation and error semantics.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// tree.iter_postorder_nodes_mut()?.try_for_each(|node| {
    ///     node?.comment = Some("visited".into());
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_postorder_nodes_mut(&mut self) -> Result<NodeTraversalMut<'_>, TreeError> {
        self.iter_postorder_nodes_mut_from(self.root()?)
    }

    /// Mutably borrow nodes of the subtree at `root` in postorder.
    ///
    /// Setup is O(tree storage size); traversal is lazy. See [`NodeTraversalMut`]
    /// for mutation, cache invalidation and error semantics.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("subtree").id;
    /// tree.iter_postorder_nodes_mut_from(root)?.try_for_each(|node| {
    ///     node?.comment = Some("visited".into());
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_postorder_nodes_mut_from(
        &mut self,
        root: NodeId,
    ) -> Result<NodeTraversalMut<'_>, TreeError> {
        self.get(root)?;
        Ok(NodeTraversalMut::new(self, root, Order::Post))
    }

    /// Borrow nodes of the tree in inorder; errors terminate iteration.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let names = tree.iter_inorder_nodes()?
    ///     .map(|node| node.map(|n| n.name.clone()))
    ///     .collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(names.len(), 5);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_inorder_nodes(
        &self,
    ) -> Result<NodeTraversal<'_, InOrderTreeTraversal<'_>>, TreeError> {
        Ok(NodeTraversal {
            tree: self,
            ids: self.iter_inorder()?,
            failed: false,
        })
    }

    /// Borrow nodes of the subtree at `root` in inorder; errors terminate iteration.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("subtree").id;
    /// let names = tree.iter_inorder_nodes_from(root)?
    ///     .map(|node| node.map(|n| n.name.clone()))
    ///     .collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(names.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_inorder_nodes_from(
        &self,
        root: NodeId,
    ) -> Result<NodeTraversal<'_, InOrderTreeTraversal<'_>>, TreeError> {
        Ok(NodeTraversal {
            tree: self,
            ids: self.iter_inorder_from(root)?,
            failed: false,
        })
    }

    /// Mutably borrow nodes of the tree in inorder.
    ///
    /// Setup is O(tree storage size); traversal is lazy. See [`NodeTraversalMut`]
    /// for mutation, cache invalidation and error semantics.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// tree.iter_inorder_nodes_mut()?.try_for_each(|node| {
    ///     node?.comment = Some("visited".into());
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_inorder_nodes_mut(&mut self) -> Result<NodeTraversalMut<'_>, TreeError> {
        self.iter_inorder_nodes_mut_from(self.root()?)
    }

    /// Mutably borrow nodes of the subtree at `root` in inorder.
    ///
    /// Setup is O(tree storage size); traversal is lazy. See [`NodeTraversalMut`]
    /// for mutation, cache invalidation and error semantics.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("subtree").id;
    /// tree.iter_inorder_nodes_mut_from(root)?.try_for_each(|node| {
    ///     node?.comment = Some("visited".into());
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_inorder_nodes_mut_from(
        &mut self,
        root: NodeId,
    ) -> Result<NodeTraversalMut<'_>, TreeError> {
        self.get(root)?;
        Ok(NodeTraversalMut::new(self, root, Order::In))
    }

    /// Borrow nodes of the tree in levelorder; errors terminate iteration.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let names = tree.iter_levelorder_nodes()?
    ///     .map(|node| node.map(|n| n.name.clone()))
    ///     .collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(names.len(), 5);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_levelorder_nodes(
        &self,
    ) -> Result<NodeTraversal<'_, LevelOrderTreeTraversal<'_>>, TreeError> {
        Ok(NodeTraversal {
            tree: self,
            ids: self.iter_levelorder()?,
            failed: false,
        })
    }

    /// Borrow nodes of the subtree at `root` in levelorder; errors terminate iteration.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    /// let tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("subtree").id;
    /// let names = tree.iter_levelorder_nodes_from(root)?
    ///     .map(|node| node.map(|n| n.name.clone()))
    ///     .collect::<Result<Vec<_>, _>>()?;
    /// assert_eq!(names.len(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_levelorder_nodes_from(
        &self,
        root: NodeId,
    ) -> Result<NodeTraversal<'_, LevelOrderTreeTraversal<'_>>, TreeError> {
        Ok(NodeTraversal {
            tree: self,
            ids: self.iter_levelorder_from(root)?,
            failed: false,
        })
    }

    /// Mutably borrow nodes of the tree in levelorder.
    ///
    /// Setup is O(tree storage size); traversal is lazy. See [`NodeTraversalMut`]
    /// for mutation, cache invalidation and error semantics.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// tree.iter_levelorder_nodes_mut()?.try_for_each(|node| {
    ///     node?.comment = Some("visited".into());
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_levelorder_nodes_mut(&mut self) -> Result<NodeTraversalMut<'_>, TreeError> {
        self.iter_levelorder_nodes_mut_from(self.root()?)
    }

    /// Mutably borrow nodes of the subtree at `root` in levelorder.
    ///
    /// Setup is O(tree storage size); traversal is lazy. See [`NodeTraversalMut`]
    /// for mutation, cache invalidation and error semantics.
    ///
    /// ```
    /// use phylotree::tree::{Tree, TreeError};
    /// let mut tree = Tree::from_newick("((A,B)C,D)R;")?;
    /// let root = tree.get_by_name("C").expect("subtree").id;
    /// tree.iter_levelorder_nodes_mut_from(root)?.try_for_each(|node| {
    ///     node?.comment = Some("visited".into());
    ///     Ok::<(), TreeError>(())
    /// })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter_levelorder_nodes_mut_from(
        &mut self,
        root: NodeId,
    ) -> Result<NodeTraversalMut<'_>, TreeError> {
        self.get(root)?;
        Ok(NodeTraversalMut::new(self, root, Order::Level))
    }
}

mod tests {

    use crate::tree::{Node, NodeId, NodeTraversalMut, Tree, TreeError};
    use std::iter::FusedIterator;

    type Traversal<'a> = Box<dyn Iterator<Item = Result<NodeId, TreeError>> + 'a>;

    fn iterators(tree: &Tree, root: NodeId) -> Vec<Traversal<'_>> {
        vec![
            Box::new(tree.iter_preorder_from(root).unwrap()),
            Box::new(tree.iter_postorder_from(root).unwrap()),
            Box::new(tree.iter_inorder_from(root).unwrap()),
            Box::new(tree.iter_levelorder_from(root).unwrap()),
        ]
    }

    fn names(tree: &Tree, ids: impl Iterator<Item = Result<NodeId, TreeError>>) -> Vec<String> {
        ids.map(|id| tree.get(id.unwrap()).unwrap().name.clone().unwrap())
            .collect()
    }

    #[test]
    fn exact_orders_unary_nodes_and_subtrees() {
        let tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
        let root = tree.root().unwrap();
        for (iterator, expected) in iterators(&tree, root).into_iter().zip([
            vec!["F", "B", "A", "D", "C", "E", "G", "I", "H"],
            vec!["A", "C", "E", "D", "B", "H", "I", "G", "F"],
            vec!["A", "B", "C", "D", "E", "F", "H", "I", "G"],
            vec!["F", "B", "G", "A", "D", "I", "C", "E", "H"],
        ]) {
            assert_eq!(names(&tree, iterator), expected);
        }
        let root = tree.get_by_name("B").unwrap().id;
        for (iterator, expected) in iterators(&tree, root).into_iter().zip([
            tree.preorder(root).unwrap(),
            tree.postorder(root).unwrap(),
            tree.inorder(root).unwrap(),
            tree.levelorder(root).unwrap(),
        ]) {
            assert_eq!(iterator.collect::<Result<Vec<_>, _>>().unwrap(), expected);
            assert_eq!(expected.len(), 5);
        }
        let leaf = tree.get_by_name("A").unwrap().id;
        for mut iterator in iterators(&tree, leaf) {
            assert_eq!(iterator.next().unwrap().unwrap(), leaf);
            assert!(iterator.next().is_none());
            assert!(iterator.next().is_none());
        }
    }

    #[test]
    fn polytomies_preserve_order_and_inorder_rejects_them_lazily() {
        let tree = Tree::from_newick("(A,(B,C,D)X,E)R;").unwrap();
        assert_eq!(
            names(&tree, tree.iter_preorder().unwrap()),
            ["R", "A", "X", "B", "C", "D", "E"]
        );
        assert_eq!(
            names(&tree, tree.iter_postorder().unwrap()),
            ["A", "B", "C", "D", "X", "E", "R"]
        );
        assert_eq!(
            names(&tree, tree.iter_levelorder().unwrap()),
            ["R", "A", "X", "E", "B", "C", "D"]
        );
        let mut iterator = tree.iter_inorder().unwrap();
        assert!(matches!(iterator.next(), Some(Err(TreeError::IsNotBinary))));
        assert!(iterator.next().is_none());
        assert!(matches!(
            tree.inorder(tree.root().unwrap()),
            Err(TreeError::IsNotBinary)
        ));
        // Binary subtrees remain traversable even if the whole tree is nonbinary.
        let leaf = tree.get_by_name("A").unwrap().id;
        assert_eq!(tree.iter_inorder_from(leaf).unwrap().count(), 1);
        let tree = Tree::from_newick("(A,(B,C,D)X)R;").unwrap();
        let mut iterator = tree.iter_inorder().unwrap();
        assert!(iterator.next().unwrap().is_ok()); // A
        assert!(iterator.next().unwrap().is_ok()); // R
        assert!(matches!(iterator.next(), Some(Err(TreeError::IsNotBinary))));
        assert!(iterator.next().is_none());
    }

    #[test]
    fn missing_nodes_are_errors_once_then_iteration_is_fused() {
        let mut tree = Tree::from_newick("(A,B)R;").unwrap();
        let root = tree.root().unwrap();
        let missing = NodeId::from(usize::MAX);
        tree.get_mut(root).unwrap().children[1] = missing;
        for mut iterator in iterators(&tree, root) {
            loop {
                match iterator.next() {
                    Some(Ok(id)) => assert_ne!(id, missing),
                    Some(Err(TreeError::NodeNotFound(id))) => {
                        assert_eq!(id, missing);
                        break;
                    }
                    other => panic!("Expected missing node error, got {other:?}"),
                }
            }
            assert!(iterator.next().is_none());
            assert!(iterator.next().is_none());
        }
        for result in [
            tree.preorder(root),
            tree.postorder(root),
            tree.inorder(root),
            tree.levelorder(root),
        ] {
            assert!(matches!(result, Err(TreeError::NodeNotFound(id)) if id == missing));
        }
        // Preorder can return a prefix without looking up any descendant yet.
        assert_eq!(
            tree.iter_preorder()
                .unwrap()
                .take(1)
                .collect::<Result<Vec<_>, _>>()
                .unwrap(),
            [root]
        );
        assert!(matches!(
            tree.iter_preorder_from(missing),
            Err(TreeError::NodeNotFound(_))
        ));
        assert!(matches!(
            tree.iter_postorder_from(missing),
            Err(TreeError::NodeNotFound(_))
        ));
        assert!(matches!(
            tree.iter_inorder_from(missing),
            Err(TreeError::NodeNotFound(_))
        ));
        assert!(matches!(
            tree.iter_levelorder_from(missing),
            Err(TreeError::NodeNotFound(_))
        ));
    }

    #[test]
    fn empty_tree_and_fused_trait() {
        let tree = Tree::new();
        assert!(matches!(tree.iter_preorder(), Err(TreeError::RootNotFound)));
        assert!(matches!(
            tree.iter_postorder(),
            Err(TreeError::RootNotFound)
        ));
        assert!(matches!(tree.iter_inorder(), Err(TreeError::RootNotFound)));
        assert!(matches!(
            tree.iter_levelorder(),
            Err(TreeError::RootNotFound)
        ));
        fn fused(_: impl FusedIterator) {}
        let tree = Tree::from_newick("(A,B)R;").unwrap();
        fused(tree.iter_preorder().unwrap());
        fused(tree.iter_postorder().unwrap());
        fused(tree.iter_inorder().unwrap());
        fused(tree.iter_levelorder().unwrap());
    }

    #[test]
    fn deep_chain_does_not_use_the_call_stack() {
        let mut tree = Tree::new();
        let root = tree.add(Node::new());
        let mut last = root;
        for _ in 1..50_000 {
            last = tree.add_child(Node::new(), last, None).unwrap();
        }
        for (iterator, reverse) in iterators(&tree, root)
            .into_iter()
            .zip([false, true, true, false])
        {
            let ids = iterator.collect::<Result<Vec<_>, _>>().unwrap();
            assert_eq!(ids.len(), 50_000);
            for (i, id) in ids.into_iter().enumerate() {
                assert_eq!(id.index(), if reverse { 49_999 - i } else { i });
            }
        }
        // The compatibility wrappers also avoid recursion.
        assert_eq!(tree.preorder(root).unwrap().len(), 50_000);
        assert_eq!(tree.postorder(root).unwrap().len(), 50_000);
        assert_eq!(tree.inorder(root).unwrap().len(), 50_000);
        assert_eq!(tree.levelorder(root).unwrap().len(), 50_000);
    }

    type Nodes<'a> = Box<dyn Iterator<Item = Result<&'a Node, TreeError>> + 'a>;
    fn borrowed(tree: &Tree, order: usize, root: Option<NodeId>) -> Result<Nodes<'_>, TreeError> {
        Ok(match (order, root) {
            (0, None) => Box::new(tree.iter_preorder_nodes()?),
            (1, None) => Box::new(tree.iter_postorder_nodes()?),
            (2, None) => Box::new(tree.iter_inorder_nodes()?),
            (3, None) => Box::new(tree.iter_levelorder_nodes()?),
            (0, Some(id)) => Box::new(tree.iter_preorder_nodes_from(id)?),
            (1, Some(id)) => Box::new(tree.iter_postorder_nodes_from(id)?),
            (2, Some(id)) => Box::new(tree.iter_inorder_nodes_from(id)?),
            (3, Some(id)) => Box::new(tree.iter_levelorder_nodes_from(id)?),
            _ => unreachable!(),
        })
    }
    fn mutable(
        tree: &mut Tree,
        order: usize,
        root: Option<NodeId>,
    ) -> Result<NodeTraversalMut<'_>, TreeError> {
        match (order, root) {
            (0, None) => tree.iter_preorder_nodes_mut(),
            (1, None) => tree.iter_postorder_nodes_mut(),
            (2, None) => tree.iter_inorder_nodes_mut(),
            (3, None) => tree.iter_levelorder_nodes_mut(),
            (0, Some(id)) => tree.iter_preorder_nodes_mut_from(id),
            (1, Some(id)) => tree.iter_postorder_nodes_mut_from(id),
            (2, Some(id)) => tree.iter_inorder_nodes_mut_from(id),
            (3, Some(id)) => tree.iter_levelorder_nodes_mut_from(id),
            _ => unreachable!(),
        }
    }
    fn ids(tree: &Tree, order: usize, root: NodeId) -> Vec<NodeId> {
        match order {
            0 => tree.preorder(root),
            1 => tree.postorder(root),
            2 => tree.inorder(root),
            3 => tree.levelorder(root),
            _ => unreachable!(),
        }
        .unwrap()
    }

    #[test]
    fn node_orders_match_ids_and_return_original_nodes() {
        for order in 0..4 {
            for subtree in [false, true] {
                let mut tree = Tree::from_newick("((A,(C,E)D)B,((H)I)G)F;").unwrap();
                let root = subtree.then(|| tree.get_by_name("B").unwrap().id);
                let expected = ids(&tree, order, root.unwrap_or(tree.root().unwrap()));
                let nodes = borrowed(&tree, order, root)
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap();
                assert_eq!(
                    nodes.iter().map(|node| node.id).collect::<Vec<_>>(),
                    expected
                );
                for node in nodes {
                    assert!(std::ptr::eq(node, tree.get(node.id).unwrap()));
                }
                // All mutable references may coexist, and mutations affect original storage.
                let mut nodes = mutable(&mut tree, order, root)
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap();
                assert_eq!(
                    nodes.iter().map(|node| node.id).collect::<Vec<_>>(),
                    expected
                );
                for (index, node) in nodes.iter_mut().enumerate() {
                    node.comment = Some(index.to_string());
                }
                drop(nodes);
                for (index, id) in expected.into_iter().enumerate() {
                    assert_eq!(
                        tree.get(id).unwrap().comment.as_deref(),
                        Some(index.to_string().as_str())
                    );
                }
                if subtree {
                    assert!(tree.get(tree.root().unwrap()).unwrap().comment.is_none());
                }
            }
        }
    }

    #[test]
    fn mutable_references_can_be_retained_across_next_and_iterator_drop() {
        for order in 0..4 {
            let mut tree = Tree::from_newick("(A,B)R;").unwrap();
            let mut iter = mutable(&mut tree, order, None).unwrap();
            let first = iter.next().unwrap().unwrap();
            let second = iter.next().unwrap().unwrap();
            assert_ne!(first.id, second.id);
            first.name = Some("first".into());
            second.name = Some("second".into());
            drop(iter);
            first.comment = Some("still borrowed".into());
            assert_eq!(second.name.as_deref(), Some("second"));
            assert_eq!(
                tree.get_by_name("first").unwrap().comment.as_deref(),
                Some("still borrowed")
            );
        }
    }

    #[test]
    fn scaling_matches_tree_rescale_including_both_branch_representations() {
        for order in 0..4 {
            let mut tree = Tree::from_newick("((A:1,B:2)C:3,(D:4,E:5)F:6)R;").unwrap();
            let mut expected = tree.clone();
            expected.rescale(2.0);
            tree.distance_matrix().unwrap(); // Populate caches before mutation.
            tree.rebuild_taxon_index().unwrap();
            mutable(&mut tree, order, None)
                .unwrap()
                .try_for_each(|node| {
                    node?.rescale_edges(2.0);
                    Ok::<(), TreeError>(())
                })
                .unwrap();
            assert!(tree.cached_taxon_index().is_none());
            assert_eq!(tree.to_newick().unwrap(), expected.to_newick().unwrap());
            assert_eq!(tree.diameter().unwrap(), expected.diameter().unwrap());
            assert_eq!(
                tree.distance_matrix().unwrap().as_slice(),
                expected.distance_matrix().unwrap().as_slice()
            );
            for node in tree.iter_preorder_nodes().unwrap() {
                let node = node.unwrap();
                for &child in &node.children {
                    assert_eq!(
                        node.get_child_edge(child),
                        tree.get(child).unwrap().parent_edge
                    );
                }
            }
        }
    }

    #[test]
    fn missing_nodes_errors_are_lazy_and_terminate_iteration() {
        for order in 0..4 {
            let mut tree = Tree::from_newick("(A,B)R;").unwrap();
            let root = tree.root().unwrap();
            let missing = NodeId::from(usize::MAX);
            tree.get_mut(root).unwrap().children[1] = missing;
            let mut read = borrowed(&tree, order, None).unwrap();
            loop {
                match read.next() {
                    Some(Ok(_)) => {}
                    Some(Err(TreeError::NodeNotFound(id))) => {
                        assert_eq!(id, missing);
                        break;
                    }
                    other => panic!("Expected missing-node error, got {other:?}"),
                }
            }
            assert!(read.next().is_none());
            drop(read);
            let mut iter = mutable(&mut tree, order, None).unwrap();
            loop {
                match iter.next() {
                    Some(Ok(_)) => {}
                    Some(Err(TreeError::NodeNotFound(id))) => {
                        assert_eq!(id, missing);
                        break;
                    }
                    other => panic!("Expected missing-node error, got {other:?}"),
                }
            }
            assert!(iter.next().is_none());
            assert!(iter.next().is_none());
            drop(iter);
            assert!(matches!(
                mutable(&mut tree, order, Some(missing)),
                Err(TreeError::NodeNotFound(_))
            ));
            assert!(matches!(
                borrowed(&tree, order, Some(missing)),
                Err(TreeError::NodeNotFound(_))
            ));
            let mut empty = Tree::new();
            assert!(matches!(
                mutable(&mut empty, order, None),
                Err(TreeError::RootNotFound)
            ));
            assert!(matches!(
                borrowed(&empty, order, None),
                Err(TreeError::RootNotFound)
            ));
        }
    }

    #[test]
    fn repeated_nodes_never_produce_aliasing_or_infinite_loops() {
        for order in 0..4 {
            for cycle in [false, true] {
                let mut tree = Tree::from_newick("(A,B)R;").unwrap();
                let root = tree.root().unwrap();
                let a = tree.get_by_name("A").unwrap().id;
                if cycle {
                    tree.get_mut(a).unwrap().children.push(root);
                } else {
                    tree.get_mut(root).unwrap().children[1] = a;
                }
                let mut iter = mutable(&mut tree, order, None).unwrap();
                let mut retained = Vec::new();
                loop {
                    match iter.next() {
                        Some(Ok(node)) => retained.push(node),
                        Some(Err(TreeError::RepeatedTraversalNode(id))) => {
                            assert_eq!(id, if cycle { root } else { a });
                            break;
                        }
                        other => panic!("Expected repeated-node error, got {other:?}"),
                    }
                }
                assert!(iter.next().is_none());
                drop(iter);
                for node in retained {
                    node.comment = Some("retained safely".into());
                }
            }
        }
    }

    #[test]
    fn edits_to_yielded_links_do_not_redirect_traversal() {
        for order in 0..4 {
            let mut tree = Tree::from_newick("((A,B)C,D)R;").unwrap();
            let expected = ids(&tree, order, tree.root().unwrap());
            let mut actual = Vec::new();
            for node in mutable(&mut tree, order, None).unwrap() {
                let node = node.unwrap();
                actual.push(node.id);
                node.children.clear(); // The remaining child work is already captured.
                node.id = NodeId::from(usize::MAX); // Slot identity does not use editable IDs.
            }
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn deleted_nodes_polytomies_and_fused_iterators() {
        for order in 0..4 {
            let mut tree = Tree::from_newick("((A,B)C,D)R;").unwrap();
            let pruned = tree.get_by_name("C").unwrap().id;
            tree.prune(pruned).unwrap();
            let nodes = mutable(&mut tree, order, None)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            assert_eq!(nodes.len(), 2);
            drop(nodes);
            // A stale link to a deleted slot is an error, not a repeated visit.
            let root = tree.root().unwrap();
            tree.get_mut(root).unwrap().children.push(pruned);
            assert!(
                matches!(mutable(&mut tree, order, None).unwrap().collect::<Result<Vec<_>, _>>(), Err(TreeError::NodeNotFound(id)) if id == pruned)
            );
        }
        let mut tree = Tree::from_newick("(A,B,C)R;").unwrap();
        assert!(matches!(
            tree.iter_inorder_nodes().unwrap().next(),
            Some(Err(TreeError::IsNotBinary))
        ));
        let mut iter = tree.iter_inorder_nodes_mut().unwrap();
        assert!(matches!(iter.next(), Some(Err(TreeError::IsNotBinary))));
        assert!(iter.next().is_none());
        drop(iter);
        for order in [0, 1, 3] {
            assert_eq!(mutable(&mut tree, order, None).unwrap().count(), 4);
        }
        fn fused(_: impl FusedIterator) {}
        fused(tree.iter_preorder_nodes().unwrap());
        fused(tree.iter_postorder_nodes().unwrap());
        fused(tree.iter_inorder_nodes().unwrap());
        fused(tree.iter_levelorder_nodes().unwrap());
        fused(tree.iter_preorder_nodes_mut().unwrap());
    }

    #[test]
    fn deep_mutable_walks_remain_non_recursive() {
        let mut tree = Tree::new();
        let mut last = tree.add(Node::new());
        for _ in 1..50_000 {
            last = tree.add_child(Node::new(), last, None).unwrap();
        }
        for order in 0..4 {
            let mut count = 0;
            mutable(&mut tree, order, None)
                .unwrap()
                .try_for_each(|node| {
                    node?.comment = Some("visited".into());
                    count += 1;
                    Ok::<(), TreeError>(())
                })
                .unwrap();
            assert_eq!(count, 50_000);
        }
    }
}
