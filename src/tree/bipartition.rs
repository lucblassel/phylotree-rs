use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;

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

type ProfileMap = HashMap<Partition, PartitionData>;

/// Selects the rooted or unrooted interpretation of tree partitions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobinsonFouldsMode {
    /// Compare directed descendant clades. Both source trees must be rooted.
    Rooted,
    /// Compare unordered canonical splits, ignoring root placement.
    Unrooted,
}

/// An immutable, computed clade representation of a tree.
///
/// Profiles always store nontrivial descendant clades and lazily retain canonical
/// unrooted splits after the first unrooted operation.
#[derive(Debug, Clone)]
pub struct BipartitionProfile {
    taxon_index: TaxonIndex,
    rooted_clades: ProfileMap,
    unrooted_splits: OnceLock<ProfileMap>,
    rooted: bool,
}

impl BipartitionProfile {
    /// Computes a profile from an immutable tree.
    ///
    /// Retain a profile when comparing one tree to many others so its lazily
    /// computed unrooted splits can be reused:
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
    ///     reference.unrooted_robinson_foulds(&sample).unwrap()
    /// }).collect();
    ///
    /// assert_eq!(distances[0], 0);
    /// assert!(distances[1] > 0);
    /// ```
    pub fn from_tree(tree: &Tree) -> Result<Self, TreeError> {
        let taxon_index = tree.taxon_index()?;
        let root = tree.get_root()?;
        let mut clades = ProfileMap::new();
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

            let clade_size = descendants.count_ones(..);
            if node.parent.is_some() && !node.is_tip() && clade_size > 1 {
                clades.insert(
                    descendants.clone(),
                    PartitionData {
                        depth: node.get_depth(),
                        length: node.parent_edge,
                    },
                );
            }
            descendant_taxa.insert(node_id, descendants);
        }

        Ok(Self {
            taxon_index,
            rooted_clades: clades,
            unrooted_splits: OnceLock::new(),
            rooted: tree.is_rooted()?,
        })
    }

    fn canonical_partition(partition: Partition) -> Partition {
        let mut complement = partition.clone();
        complement.toggle_range(..);
        complement.min(partition)
    }

    fn merge_split(splits: &mut ProfileMap, partition: Partition, data: PartitionData) {
        use std::collections::hash_map::Entry;

        match splits.entry(partition) {
            Entry::Vacant(entry) => {
                entry.insert(data);
            }
            Entry::Occupied(mut entry) => {
                let old = entry.get_mut();
                old.depth = old.depth.min(data.depth);
                old.length = match (old.length, data.length) {
                    (Some(left), Some(right)) => Some(left + right),
                    _ => None,
                };
            }
        }
    }

    fn project_unrooted(clades: &ProfileMap) -> ProfileMap {
        let mut splits = ProfileMap::with_capacity(clades.len());
        for (clade, data) in clades {
            let partition = Self::canonical_partition(clade.clone());
            if partition.count_ones(..) > 1 {
                Self::merge_split(&mut splits, partition, *data);
            }
        }
        splits
    }

    fn require_compatible(&self, other: &Self) -> Result<(), TreeError> {
        if self.taxon_index.is_compatible(&other.taxon_index) {
            Ok(())
        } else {
            Err(TreeError::DifferentTipIndices)
        }
    }

    fn partitions(&self, mode: RobinsonFouldsMode) -> Result<&ProfileMap, TreeError> {
        match mode {
            RobinsonFouldsMode::Rooted => {
                if !self.rooted {
                    return Err(TreeError::RootedComparisonRequiresRootedTrees);
                }
                Ok(&self.rooted_clades)
            }
            RobinsonFouldsMode::Unrooted => Ok(self
                .unrooted_splits
                .get_or_init(|| Self::project_unrooted(&self.rooted_clades))),
        }
    }

    fn require_branch_lengths(partitions: &ProfileMap) -> Result<(), TreeError> {
        if partitions.values().all(|data| data.length.is_some()) {
            Ok(())
        } else {
            Err(TreeError::MissingBranchLengths)
        }
    }

    pub(crate) fn partitions_with_lengths(&self) -> Result<PartitionMap, TreeError> {
        self.partitions(RobinsonFouldsMode::Unrooted)?
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

    /// Returns the number of stored nontrivial descendant clades.
    pub fn len(&self) -> usize {
        self.rooted_clades.len()
    }

    /// Returns `true` if the profile contains no nontrivial descendant clades.
    pub fn is_empty(&self) -> bool {
        self.rooted_clades.is_empty()
    }

    /// Counts groupings shared with another compatible profile in the selected mode.
    pub fn common_partition_count(
        &self,
        other: &Self,
        mode: RobinsonFouldsMode,
    ) -> Result<usize, TreeError> {
        self.require_compatible(other)?;
        let left = self.partitions(mode)?;
        let right = other.partitions(mode)?;
        Ok(left
            .keys()
            .filter(|partition| right.contains_key(*partition))
            .count())
    }

    pub(crate) fn partition_set(
        &self,
        mode: RobinsonFouldsMode,
    ) -> Result<PartitionSet, TreeError> {
        Ok(self.partitions(mode)?.keys().cloned().collect())
    }

    /// Computes the Robinson–Foulds distance in the selected mode.
    pub fn robinson_foulds(
        &self,
        other: &Self,
        mode: RobinsonFouldsMode,
    ) -> Result<usize, TreeError> {
        self.require_compatible(other)?;
        let left = self.partitions(mode)?;
        let right = other.partitions(mode)?;
        let intersection = left
            .keys()
            .filter(|partition| right.contains_key(*partition))
            .count();
        Ok(left.len() + right.len() - 2 * intersection)
    }

    /// Computes rooted Robinson–Foulds distance from descendant clades.
    pub fn rooted_robinson_foulds(&self, other: &Self) -> Result<usize, TreeError> {
        self.robinson_foulds(other, RobinsonFouldsMode::Rooted)
    }

    /// Computes unrooted Robinson–Foulds distance from canonical splits.
    pub fn unrooted_robinson_foulds(&self, other: &Self) -> Result<usize, TreeError> {
        self.robinson_foulds(other, RobinsonFouldsMode::Unrooted)
    }

    /// Computes normalized Robinson–Foulds distance in the selected mode.
    pub fn robinson_foulds_norm(
        &self,
        other: &Self,
        mode: RobinsonFouldsMode,
    ) -> Result<f64, TreeError> {
        self.require_compatible(other)?;
        let left = self.partitions(mode)?;
        let right = other.partitions(mode)?;
        let intersection = left
            .keys()
            .filter(|partition| right.contains_key(*partition))
            .count();
        let total = left.len() + right.len();
        if total == 0 {
            Ok(0.0)
        } else {
            Ok((total - 2 * intersection) as f64 / total as f64)
        }
    }

    /// Computes weighted Robinson–Foulds distance in the selected mode.
    pub fn weighted_robinson_foulds(
        &self,
        other: &Self,
        mode: RobinsonFouldsMode,
    ) -> Result<f64, TreeError> {
        self.require_compatible(other)?;
        let left = self.partitions(mode)?;
        let right = other.partitions(mode)?;
        Self::require_branch_lengths(left)?;
        Self::require_branch_lengths(right)?;
        let mut distance = 0.0;

        for (partition, left_data) in left.iter() {
            let left_length = left_data.length.unwrap();
            distance += match right.get(partition) {
                Some(right_data) => (left_length - right_data.length.unwrap()).abs(),
                None => left_length,
            };
        }
        for (partition, right_data) in right.iter() {
            if !left.contains_key(partition) {
                distance += right_data.length.unwrap();
            }
        }
        Ok(distance)
    }

    /// Computes the Kuhner–Felsenstein branch score in the selected mode.
    pub fn kuhner_felsenstein(
        &self,
        other: &Self,
        mode: RobinsonFouldsMode,
    ) -> Result<f64, TreeError> {
        self.require_compatible(other)?;
        let left = self.partitions(mode)?;
        let right = other.partitions(mode)?;
        Self::require_branch_lengths(left)?;
        Self::require_branch_lengths(right)?;
        let mut distance = 0.0;

        for (partition, left_data) in left.iter() {
            let left_length = left_data.length.unwrap();
            distance += match right.get(partition) {
                Some(right_data) => f64::powi(left_length - right_data.length.unwrap(), 2),
                None => f64::powi(left_length, 2),
            };
        }
        for (partition, right_data) in right.iter() {
            if !left.contains_key(partition) {
                distance += f64::powi(right_data.length.unwrap(), 2);
            }
        }
        Ok(distance.sqrt())
    }

    /// Computes topology and branch-length comparison metrics in one pass.
    pub fn compare(&self, other: &Self, mode: RobinsonFouldsMode) -> Result<Comparison, TreeError> {
        self.require_compatible(other)?;
        let left = self.partitions(mode)?;
        let right = other.partitions(mode)?;
        Self::require_branch_lengths(left)?;
        Self::require_branch_lengths(right)?;
        let total = left.len() + right.len();
        let mut intersection = 0.0;
        let mut weighted_rf = 0.0;
        let mut branch_score = 0.0;

        for (partition, left_data) in left.iter() {
            let left_length = left_data.length.unwrap();
            if let Some(right_data) = right.get(partition) {
                let right_length = right_data.length.unwrap();
                weighted_rf += (left_length - right_length).abs();
                branch_score += f64::powi(left_length - right_length, 2);
                intersection += 1.0;
            } else {
                weighted_rf += left_length;
                branch_score += f64::powi(left_length, 2);
            }
        }
        for (partition, right_data) in right.iter() {
            if !left.contains_key(partition) {
                let right_length = right_data.length.unwrap();
                weighted_rf += right_length;
                branch_score += f64::powi(right_length, 2);
            }
        }

        let rf = total as f64 - 2.0 * intersection;
        Ok(Comparison {
            rf,
            norm_rf: if total == 0 { 0.0 } else { rf / total as f64 },
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

        assert_eq!(
            first_profile
                .robinson_foulds(&second_profile, RobinsonFouldsMode::Rooted)
                .unwrap(),
            0
        );
        assert_eq!(
            first_profile
                .weighted_robinson_foulds(&second_profile, RobinsonFouldsMode::Unrooted)
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
            first_profile.robinson_foulds(&second_profile, RobinsonFouldsMode::Unrooted),
            Err(TreeError::DifferentTipIndices)
        ));
    }

    #[test]
    fn storage_strategies_produce_identical_unrooted_results() {
        let first = Tree::from_newick("((A,B),(C,D));").unwrap();
        let second = Tree::from_newick("((A,C),(B,D));").unwrap();
        let first_clades = BipartitionProfile::from_tree(&first).unwrap();
        let first_cached = BipartitionProfile::from_tree(&first).unwrap();
        let second_clades = BipartitionProfile::from_tree(&second).unwrap();
        let second_cached = BipartitionProfile::from_tree(&second).unwrap();

        let expected = first_clades
            .unrooted_robinson_foulds(&second_clades)
            .unwrap();
        assert_eq!(
            first_clades
                .unrooted_robinson_foulds(&second_cached)
                .unwrap(),
            expected
        );
        assert_eq!(
            first_cached
                .unrooted_robinson_foulds(&second_clades)
                .unwrap(),
            expected
        );
        assert_eq!(
            first_cached
                .unrooted_robinson_foulds(&second_cached)
                .unwrap(),
            expected
        );
    }

    #[test]
    fn rooted_rf_detects_root_placement_while_unrooted_rf_ignores_it() {
        let first = Tree::from_newick("((A,B),(C,D));").unwrap();
        let second = Tree::from_newick("(A,(B,(C,D)));").unwrap();
        let first = BipartitionProfile::from_tree(&first).unwrap();
        let second = BipartitionProfile::from_tree(&second).unwrap();

        assert_eq!(first.unrooted_robinson_foulds(&second).unwrap(), 0);
        assert_eq!(first.rooted_robinson_foulds(&second).unwrap(), 2);
    }

    #[test]
    fn rooted_rf_rejects_an_unrooted_tree() {
        let rooted = Tree::from_newick("((A,B),(C,D));").unwrap();
        let unrooted = Tree::from_newick("(A,B,(C,D));").unwrap();
        let rooted = BipartitionProfile::from_tree(&rooted).unwrap();
        let unrooted = BipartitionProfile::from_tree(&unrooted).unwrap();

        assert!(matches!(
            rooted.rooted_robinson_foulds(&unrooted),
            Err(TreeError::RootedComparisonRequiresRootedTrees)
        ));
        assert_eq!(rooted.unrooted_robinson_foulds(&unrooted).unwrap(), 0);
    }
}
