from re import A
from .pytree import Tree as _rsTree
from typing import Optional


class Tree:
    """A Class representing a Phylogenetic tree"""

    def __init__(self) -> None:
        """Create an empty Tree"""
        self._rs = _rsTree()

    def copy(self) -> "Tree":
        """Create a deep copy of the tree"""
        tree = Tree()
        tree._rs = self._rs.copy()

        return tree

    @staticmethod
    def from_newick(path: str) -> "Tree":
        """Create a Tree from a newick formatted file"""
        tree = Tree()
        tree._rs = _rsTree.from_newick(path)

        return tree

    @staticmethod
    def from_string(newick: str) -> "Tree":
        """Create a Tree from a newick formatted string"""
        tree = Tree()
        tree._rs = _rsTree.from_string(newick)

        return tree

    def to_newick(self) -> str:
        """Represent the Tree as a newick formatted string"""
        return self._rs.to_string()

    def to_file(self, path: str):
        """Write the tree to a newick formatted file"""
        return self._rs.to_file(path)

    @property
    def is_binary(self) -> bool:
        """Returns true if the Tree is binary"""
        return self._rs.is_binary()

    @property
    def is_rooted(self) -> bool:
        """Returns true if the Tree is rooted"""
        return self._rs.is_rooted()

    def height(self) -> float:
        """
        If the tree has branch lengths, this function will
        return the sum of branch lengths in the path between
        the root and the deepest node. If the branches are unweighted
        then it will return the number of branches in that path.
        """
        return self._rs.height()

    def diameter(self) -> float:
        """
        If the tree has branch lengths, this function will
        return the sum of branch lengths in the path between
        the two furthest nodes. If the branches are unweighted
        then it will return the number of branches in that path.
        """
        return self._rs.diameter()

    def n_tips(self) -> int:
        """Returns the number of tips in the Tree"""
        return self._rs.n_tips()

    def n_nodes(self) -> int:
        """Returns the total number of nodes in the Tree"""
        return self._rs.n_nodes()

    def n_cherries(self) -> int:
        """Returns the number of cherries in the Tree"""
        return self._rs.n_cherries()

    def compare(self, other: "Tree") -> dict:
        """
        Compares two trees according to topology and branch lenghts.
        It returns a dict containing:
          - rf: The Robinson-Foulds distance
          - norm_rf: The normalized Robinson-Foulds distance
          - weighted_rf: The Robinson-Foulds distance, weighted by branch lengths
          - branch_score: The Khuner-Felsenstein branch score
        """
        return self._rs.compare(other._rs)

    def colless(self, normalisation: Optional[str] = None) -> float:
        """
        Computes the Colless index of the tree.
        The result can be normalized under different null models
        by specifying the normalisation parameter:
          - "yule": Yule null model
          - "pda": PDA null model
        """
        if normalisation is not None and normalisation not in ["yule", "pda"]:
            raise ValueError("Normalisation must be None, 'yule' or 'pda'")
        return self._rs.colless(normalisation)

    def sackin(self, normalisation: Optional[str] = None) -> float:
        """
        Computes the Sackin index of the tree.
        The result can be normalized under different null models
        by specifying the normalisation parameter:
          - "yule": Yule null model
          - "pda": PDA null model
        """
        if normalisation is not None and normalisation not in ["yule", "pda"]:
            raise ValueError("Normalisation must be None, 'yule' or 'pda'")
        return self._rs.sackin(normalisation)

    def compress(self):
        """
        Compresses the tree by removing nodes with exactly 1 parent and 1 child, while
        preserving branch lengths.
        """
        self._rs.compress()

    def rescale(self, factor: float):
        """
        Apply a multiplicative rescaling factor to all the branch lengths of the tree.
        """
        self._rs.rescale(factor)

    def prune(self, id: Optional[int] = None, name: Optional[str] = None):
        """
        Remove the subtree starting at a given root node.
        You can either choose the root by specifying:
          - The node id of the root
          - The name of the root if it exists
        If the id or name does not exist in the tree this can raise an error
        """
        if id is not None and name is not None:
            raise ValueError("You must specify only one node id or node name")

        if id is not None:
            self._rs.prune(id)
        elif name is not None:
            id = self._rs.get_name_index(name)
            self._rs.prune(id)
        else:
            raise ValueError("You must specify either a node id or a node name")

    def get_distance(
        self,
        ids: Optional[tuple[int, int]] = None,
        names: Optional[tuple[str, str]] = None,
    ) -> tuple[Optional[float], int]:
        """
        Compute the distance between two nodes in the tree.
        You can either specify a pair of node ids, or a pair of node names.
        This returns a tuple containing:
          - the sum of branch lengths on the path between the nodes
            (None if the tree has no branch lengths)
          - the topological distance between the nodes
            (i.e. the number of branches in the path)
        """
        if ids is not None and names is not None:
            raise ValueError("You must specify only one pair of node ids or names")

        if ids is not None:
            brlens, topo = self._rs.get_distance(*ids)
        elif names is not None:
            id1 = self._rs.get_name_index(names[0])
            id2 = self._rs.get_name_index(names[1])
            brlens, topo = self._rs.get_distance(id1, id2)
        else:
            raise ValueError(
                "You must specify either a pair of node ids or a node names"
            )

        return brlens, topo

    def get_id(self, name: str) -> Optional[int]:
        """Get the node id of a node by name"""
        return self._rs.get_name_index(name)

    def get_leaf_names(self) -> list[Optional[str]]:
        """Return a vector containing the names of leaf nodes"""
        return self._rs.get_leaf_names()

    def get_node_attributes(
        self, id: Optional[int] = None, name: Optional[str] = None
    ) -> dict:
        """
        Returns a dictionnary containing attributes of a specified tree node.
        The attribute dictionnary will contain:
          - id: The node id,
          - name: The node name (None if it has no name)
          - parent: The id of the parent node (None if the node is the root)
          - children: a list of node ids corresponding to the children of this node
          - parent_edge: a float denoting the length of the branch between this node and
          its parent (None if the branch has no length or if the node is the root)
          - comment: A comment associated to the node (None if there is no comment)
        """
        if id is not None and name is not None:
            raise ValueError("You cannot specify both a root id and a root name")
        if id is None and name is None:
            raise ValueError(
                "You must specify the node either with a name or a node id"
            )

        id_ = None
        if name is not None:
            id_ = self._rs.get_name_index(name)
        if id is not None:
            id_ = id

        if id_ is None:
            raise ValueError("The specified root was not found")

        return self._rs.get_node_attributes(id_)

    def traversal(
        self,
        root_id: Optional[int] = None,
        root_name: Optional[str] = None,
        order: str = "inorder",
    ) -> list[int]:
        """
        Returns ids of tree nodes ordered by visiting order in a tree traversal.
        You can specify the type of tree traversal with the `order` parameter:
          - inorder
          - preorder
          - postorder
          - levelorder
        Optional you can specify a root node for the traversal either with the root
        id or the root name.
        """

        if root_id is not None and root_name is not None:
            raise ValueError("You cannot specify both a root id and a root name")

        id = self._rs.get_root_id()
        if root_name is not None:
            id = self._rs.get_name_index(root_name)
        if root_id is not None:
            id = root_id

        if id is None:
            raise ValueError("The specified root was not found")

        if order == "inorder":
            return self._rs.inorder(id)
        elif order == "preorder":
            return self._rs.preorder(id)
        elif order == "postorder":
            return self._rs.postorder(id)
        elif order == "levelorder":
            return self._rs.levelorder(id)
        else:
            raise ValueError(f"{order} is not a valid traversal order")

    def to_matrix(self) -> dict[tuple[str, str], float]:
        """
        Returns the distance matrix associated to the tree as
        a dictionnary. The keys are taxon pairs and the values are
        the distance.
        """
        return self._rs.to_matrix()
