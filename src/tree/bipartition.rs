use std::collections::{HashMap, HashSet};

use fixedbitset::FixedBitSet;

use super::taxon_index::TaxonIndex;
use super::tree_impl::{Comparison, Tree, TreeError};
use super::{EdgeDepth, EdgeLength, NodeId};

pub(crate) type Partition = FixedBitSet;
pub(crate) type PartitionMap = HashMap<Partition, (EdgeDepth, EdgeLength)>;
pub(crate) type PartitionSet = HashSet<Partition>;

#[derive(Debug, Clone, Copy)]
struct PartitionData {
    depth: EdgeDepth,
    length: Option<EdgeLength>,
}

/// An immutable, computed bipartition representation of a tree.
///
/// Profiles are useful when a tree participates in repeated comparisons: the
/// tree is traversed once when the profile is created, and subsequent
/// comparisons operate only on hashed [`FixedBitSet`] partitions.
#[derive(Debug, Clone)]
pub struct BipartitionProfile {
    taxon_index: TaxonIndex,
    partitions: HashMap<Partition, PartitionData>,
    root_partitions: PartitionSet,
    rooted: bool,
}

impl BipartitionProfile {
    /// Computes a profile from an immutable tree without modifying or caching
    /// anything in the tree.
    ///
    /// Retain a profile when comparing one tree to many others so its partitions
    /// are computed only once:
    ///
    /// ```
    /// use phylotree::tree::{BipartitionProfile, Tree, TreeList};
    ///
    /// let reference = Tree::from_newick("((A,B),(C,D));").unwrap();
    /// let samples = TreeList::from_newick(
    ///     "((B,A),(D,C)); ((A,C),(B,D));",
    /// ).unwrap();
    /// let reference = BipartitionProfile::from_tree(&reference).unwrap();
    /// let distances: Vec<_> = samples.iter().map(|tree| {
    ///     let sample = BipartitionProfile::from_tree(tree).unwrap();
    ///     reference.robinson_foulds(&sample).unwrap()
    /// }).collect();
    ///
    /// assert_eq!(distances[0], 0);
    /// assert!(distances[1] > 0);
    /// ```
    pub fn from_tree(tree: &Tree) -> Result<Self, TreeError> {
        let taxon_index = tree.taxon_index()?;
        let root = tree.get_root()?;
        let mut partitions: HashMap<Partition, PartitionData> = HashMap::new();
        let mut descendant_taxa: HashMap<NodeId, Partition> = HashMap::new();

        for node_id in tree.postorder(&root)? {
            let node = tree.get(&node_id)?;
            let mut descendants = FixedBitSet::with_capacity(taxon_index.len());
            if node.is_tip() {
                let name = node.name.as_deref().ok_or(TreeError::UnnamedLeaves)?;
                descendants.insert(taxon_index.index_of(name).unwrap());
            } else {
                for child in &node.children {
                    descendants.union_with(descendant_taxa.get(child).unwrap());
                }
            }

            if node.parent.is_some() && !node.is_tip() {
                let partition = Self::canonical_partition(descendants.clone());
                if partition.count_ones(..) > 1 {
                    let old = partitions.get(&partition).copied();
                    let length = match (node.parent_edge, old) {
                        (None, None) => None,
                        (Some(new), Some(old)) => old.length.map(|old| old + new),
                        (Some(new), None) => Some(new),
                        (None, Some(old)) => old.length,
                    };
                    partitions.insert(
                        partition,
                        PartitionData {
                            depth: node.get_depth(),
                            length,
                        },
                    );
                }
            }
            descendant_taxa.insert(node_id, descendants);
        }

        let root_partitions = tree
            .get(&root)?
            .children
            .iter()
            .map(|child| Self::canonical_partition(descendant_taxa[child].clone()))
            .collect();

        Ok(Self {
            taxon_index,
            partitions,
            root_partitions,
            rooted: tree.is_rooted()?,
        })
    }

    fn canonical_partition(partition: Partition) -> Partition {
        let mut complement = partition.clone();
        complement.toggle_range(..);
        complement.min(partition)
    }

    fn require_compatible(&self, other: &Self) -> Result<(), TreeError> {
        if self.taxon_index.is_compatible(&other.taxon_index) {
            Ok(())
        } else {
            Err(TreeError::DifferentTipIndices)
        }
    }

    fn require_branch_lengths(&self) -> Result<(), TreeError> {
        if self.partitions.values().all(|data| data.length.is_some()) {
            Ok(())
        } else {
            Err(TreeError::MissingBranchLengths)
        }
    }

    pub(crate) fn partitions_with_lengths(&self) -> Result<PartitionMap, TreeError> {
        self.partitions
            .iter()
            .map(|(partition, data)| {
                let length = data.length.ok_or(TreeError::MissingBranchLengths)?;
                Ok((partition.clone(), (data.depth, length)))
            })
            .collect()
    }

    /// Returns the taxon-to-bit mapping used by this profile.
    pub fn taxon_index(&self) -> &TaxonIndex {
        &self.taxon_index
    }

    /// Returns whether the source tree was rooted.
    pub fn is_rooted(&self) -> bool {
        self.rooted
    }

    /// Returns the number of nontrivial partitions in the profile.
    pub fn len(&self) -> usize {
        self.partitions.len()
    }

    /// Returns `true` if the profile contains no nontrivial partitions.
    pub fn is_empty(&self) -> bool {
        self.partitions.is_empty()
    }

    /// Counts partitions shared with another compatible profile.
    pub fn common_partition_count(&self, other: &Self) -> Result<usize, TreeError> {
        self.require_compatible(other)?;
        Ok(self
            .partitions
            .keys()
            .filter(|partition| other.partitions.contains_key(*partition))
            .count())
    }

    pub(crate) fn partition_bits(&self) -> impl Iterator<Item = &Partition> {
        self.partitions.keys()
    }

    pub(crate) fn root_partition_bits(&self) -> impl Iterator<Item = &Partition> {
        self.root_partitions.iter()
    }

    pub(crate) fn partition_set(&self) -> PartitionSet {
        self.partitions.keys().cloned().collect()
    }

    /// Computes the Robinson–Foulds distance to another prepared profile.
    pub fn robinson_foulds(&self, other: &Self) -> Result<usize, TreeError> {
        let intersection = self.common_partition_count(other)?;
        let rf = self.partitions.len() + other.partitions.len() - 2 * intersection;
        let same_root = self.root_partitions == other.root_partitions;

        if self.rooted && other.rooted && rf != 0 && !same_root {
            Ok(rf + 2)
        } else {
            Ok(rf)
        }
    }

    /// Computes the normalized Robinson–Foulds distance to another profile.
    pub fn robinson_foulds_norm(&self, other: &Self) -> Result<f64, TreeError> {
        let rf = self.robinson_foulds(other)?;
        let total = self.partitions.len() + other.partitions.len();
        Ok(rf as f64 / total as f64)
    }

    /// Computes the weighted Robinson–Foulds distance to another profile.
    pub fn weighted_robinson_foulds(&self, other: &Self) -> Result<f64, TreeError> {
        self.require_compatible(other)?;
        self.require_branch_lengths()?;
        other.require_branch_lengths()?;
        let mut distance = 0.0;

        for (partition, left) in &self.partitions {
            let left_length = left.length.unwrap();
            distance += match other.partitions.get(partition) {
                Some(right) => (left_length - right.length.unwrap()).abs(),
                None => left_length,
            };
        }
        for (partition, right) in &other.partitions {
            if !self.partitions.contains_key(partition) {
                distance += right.length.unwrap();
            }
        }
        Ok(distance)
    }

    /// Computes the Kuhner–Felsenstein branch score to another profile.
    pub fn kuhner_felsenstein(&self, other: &Self) -> Result<f64, TreeError> {
        self.require_compatible(other)?;
        self.require_branch_lengths()?;
        other.require_branch_lengths()?;
        let mut distance = 0.0;

        for (partition, left) in &self.partitions {
            let left_length = left.length.unwrap();
            distance += match other.partitions.get(partition) {
                Some(right) => f64::powi(left_length - right.length.unwrap(), 2),
                None => f64::powi(left_length, 2),
            };
        }
        for (partition, right) in &other.partitions {
            if !self.partitions.contains_key(partition) {
                distance += f64::powi(right.length.unwrap(), 2);
            }
        }
        Ok(distance.sqrt())
    }

    /// Computes topology and branch-length comparison metrics in one pass.
    pub fn compare(&self, other: &Self) -> Result<Comparison, TreeError> {
        self.require_compatible(other)?;
        self.require_branch_lengths()?;
        other.require_branch_lengths()?;
        let total = self.partitions.len() + other.partitions.len();
        let mut intersection = 0.0;
        let mut weighted_rf = 0.0;
        let mut branch_score = 0.0;

        for (partition, left) in &self.partitions {
            let left_length = left.length.unwrap();
            if let Some(right) = other.partitions.get(partition) {
                let right_length = right.length.unwrap();
                weighted_rf += (left_length - right_length).abs();
                branch_score += f64::powi(left_length - right_length, 2);
                intersection += 1.0;
            } else {
                weighted_rf += left_length;
                branch_score += f64::powi(left_length, 2);
            }
        }
        for (partition, right) in &other.partitions {
            if !self.partitions.contains_key(partition) {
                let right_length = right.length.unwrap();
                weighted_rf += right_length;
                branch_score += f64::powi(right_length, 2);
            }
        }

        let mut rf = total as f64 - 2.0 * intersection;
        if self.rooted && other.rooted && rf != 0.0 && self.root_partitions != other.root_partitions
        {
            rf += 2.0;
        }

        Ok(Comparison {
            rf,
            norm_rf: rf / total as f64,
            weighted_rf,
            branch_score: branch_score.sqrt(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profiles_compare_without_mutating_trees() {
        let first = Tree::from_newick("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let second = Tree::from_newick("((B:1,A:1):1,(D:1,C:1):1);").unwrap();
        let first_profile = BipartitionProfile::from_tree(&first).unwrap();
        let second_profile = BipartitionProfile::from_tree(&second).unwrap();

        assert_eq!(first_profile.robinson_foulds(&second_profile).unwrap(), 0);
        assert_eq!(
            first_profile
                .weighted_robinson_foulds(&second_profile)
                .unwrap(),
            0.0
        );
    }

    #[test]
    fn profile_creation_does_not_populate_an_invalidated_tree_cache() {
        let mut tree = Tree::from_newick("((A,B),C);").unwrap();
        let root = tree.get_root().unwrap();
        tree.get_mut(&root).unwrap().support = Some(1.0);
        assert!(tree.cached_taxon_index().is_none());

        BipartitionProfile::from_tree(&tree).unwrap();
        assert!(tree.cached_taxon_index().is_none());

        tree.rebuild_taxon_index().unwrap();
        assert!(tree.cached_taxon_index().is_some());
    }

    #[test]
    fn profiles_reject_incompatible_taxa() {
        let first = Tree::from_newick("((A,B),C);").unwrap();
        let second = Tree::from_newick("((A,B),D);").unwrap();
        let first_profile = BipartitionProfile::from_tree(&first).unwrap();
        let second_profile = BipartitionProfile::from_tree(&second).unwrap();

        assert!(matches!(
            first_profile.robinson_foulds(&second_profile),
            Err(TreeError::DifferentTipIndices)
        ));
    }
}
