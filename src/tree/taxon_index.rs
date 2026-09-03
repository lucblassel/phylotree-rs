use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

#[derive(Debug)]
struct TaxonIndexInner {
    taxa: Vec<String>,
    indices: HashMap<String, usize>,
    mapping_hash: u64,
}

/// A shared bidirectional mapping between taxon labels and bipartition bit positions.
///
/// Taxon indices have value semantics: independently constructed indices are
/// compatible and compare equal when they assign the same labels to the same
/// positions. Clones share their underlying allocation.
#[derive(Debug, Clone)]
pub struct TaxonIndex {
    inner: Arc<TaxonIndexInner>,
}

impl TaxonIndex {
    pub(crate) fn from_sorted_taxa(taxa: Vec<String>) -> Self {
        let indices = taxa
            .iter()
            .enumerate()
            .map(|(index, taxon)| (taxon.clone(), index))
            .collect();
        let mut hasher = DefaultHasher::new();
        taxa.hash(&mut hasher);
        let mapping_hash = hasher.finish();

        Self {
            inner: Arc::new(TaxonIndexInner {
                taxa,
                indices,
                mapping_hash,
            }),
        }
    }

    pub(crate) fn empty() -> Self {
        Self::from_sorted_taxa(Vec::new())
    }

    /// Returns the number of taxa in the index.
    pub fn len(&self) -> usize {
        self.inner.taxa.len()
    }

    /// Returns `true` if the index contains no taxa.
    pub fn is_empty(&self) -> bool {
        self.inner.taxa.is_empty()
    }

    /// Returns all taxa in bit-index order.
    pub fn taxa(&self) -> &[String] {
        &self.inner.taxa
    }

    /// Returns the taxon assigned to a bit index.
    pub fn taxon(&self, index: usize) -> Option<&str> {
        self.inner.taxa.get(index).map(String::as_str)
    }

    /// Returns the bit index assigned to a taxon label.
    pub fn index_of(&self, taxon: &str) -> Option<usize> {
        self.inner.indices.get(taxon).copied()
    }

    /// Returns whether both indices assign the same taxon to every bit position.
    ///
    /// Independently parsed trees are compatible when they have the same taxa,
    /// regardless of their leaf or child order.
    ///
    /// ```
    /// use phylotree::tree::Tree;
    ///
    /// let first = Tree::from_newick("((A,B),C);").unwrap();
    /// let second = Tree::from_newick("(C,(B,A));").unwrap();
    ///
    /// assert!(first.taxon_index().unwrap().is_compatible(&second.taxon_index().unwrap()));
    /// ```
    pub fn is_compatible(&self, other: &Self) -> bool {
        self == other
    }
}

impl PartialEq for TaxonIndex {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner)
            || (self.inner.mapping_hash == other.inner.mapping_hash
                && self.inner.taxa == other.inner.taxa)
    }
}

impl Eq for TaxonIndex {}

impl Hash for TaxonIndex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.inner.mapping_hash.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn independently_constructed_equal_mappings_are_compatible() {
        let first = TaxonIndex::from_sorted_taxa(vec!["A".into(), "B".into()]);
        let second = TaxonIndex::from_sorted_taxa(vec!["A".into(), "B".into()]);
        let incompatible = TaxonIndex::from_sorted_taxa(vec!["A".into(), "C".into()]);

        assert_eq!(first, second);
        assert!(first.is_compatible(&second));
        assert!(!first.is_compatible(&incompatible));
    }

    #[test]
    fn maps_labels_and_positions_in_both_directions() {
        let index = TaxonIndex::from_sorted_taxa(vec!["A".into(), "B".into()]);

        assert_eq!(index.index_of("B"), Some(1));
        assert_eq!(index.taxon(0), Some("A"));
        assert_eq!(index.index_of("C"), None);
    }
}
