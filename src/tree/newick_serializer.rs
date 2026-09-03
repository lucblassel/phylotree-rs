use super::{Node, NodeId, Tree, TreeError};

/// Selects which node fields are included in Newick output.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum NewickFormat {
    /// Output all supported and available fields.
    AllFields,
    /// Output only the tree topology.
    Topology,
    /// Output all fields except comments.
    NoComments,
    /// Output all node names.
    OnlyNames,
    /// Output all branch lengths.
    OnlyLengths,
    /// Output leaf branch lengths and all node names.
    LeafLengthsAllNames,
    /// Output leaf branch lengths and leaf node names.
    LeafLengthsLeafNames,
    /// Output internal branch lengths and leaf node names.
    InternalLengthsLeafNames,
    /// Output all branch lengths and leaf node names.
    AllLengthsLeafNames,
}

/// Serializes a tree to Newick using a selected output format.
pub(crate) struct NewickSerializer<'a> {
    tree: &'a Tree,
    format: NewickFormat,
}

impl<'a> NewickSerializer<'a> {
    /// Creates a serializer for `tree` using `format`.
    pub(crate) fn new(tree: &'a Tree, format: NewickFormat) -> Self {
        Self { tree, format }
    }

    /// Serializes the complete tree, including the terminating semicolon.
    pub(crate) fn serialize(&self) -> Result<String, TreeError> {
        let root = self.tree.get_root()?;
        let mut output = self.serialize_subtree(&root)?;
        output.push(';');
        Ok(output)
    }

    fn serialize_subtree(&self, node_id: &NodeId) -> Result<String, TreeError> {
        let node = self.tree.get(node_id)?;
        let mut output = String::new();

        if !node.children.is_empty() {
            output.push('(');
            for (index, child_id) in node.children.iter().enumerate() {
                if index > 0 {
                    output.push(',');
                }
                output.push_str(&self.serialize_subtree(child_id)?);
            }
            output.push(')');
        }

        output.push_str(&serialize_node(node, self.format));
        Ok(output)
    }
}

/// Serializes the fields belonging to one node, without its children.
pub(crate) fn serialize_node(node: &Node, format: NewickFormat) -> String {
    let mut output = String::new();

    if includes_name(format, node.is_tip()) {
        if let Some(name) = &node.name {
            output.push_str(&serialize_name(name));
        }
    }

    if includes_branch_length(format, node.is_tip()) {
        if let Some(length) = node.parent_edge {
            output.push(':');
            output.push_str(&length.to_string());
        }
    }

    if format == NewickFormat::AllFields {
        if let Some(comment) = &node.comment {
            output.push('[');
            output.push_str(&serialize_comment(comment));
            output.push(']');
        }
    }

    output
}

fn includes_name(format: NewickFormat, is_tip: bool) -> bool {
    match format {
        NewickFormat::AllFields
        | NewickFormat::NoComments
        | NewickFormat::OnlyNames
        | NewickFormat::LeafLengthsAllNames => true,
        NewickFormat::LeafLengthsLeafNames
        | NewickFormat::InternalLengthsLeafNames
        | NewickFormat::AllLengthsLeafNames => is_tip,
        NewickFormat::Topology | NewickFormat::OnlyLengths => false,
    }
}

fn includes_branch_length(format: NewickFormat, is_tip: bool) -> bool {
    match format {
        NewickFormat::AllFields
        | NewickFormat::NoComments
        | NewickFormat::OnlyLengths
        | NewickFormat::AllLengthsLeafNames => true,
        NewickFormat::InternalLengthsLeafNames => !is_tip,
        NewickFormat::LeafLengthsLeafNames | NewickFormat::LeafLengthsAllNames => is_tip,
        NewickFormat::Topology | NewickFormat::OnlyNames => false,
    }
}

fn serialize_name(name: &str) -> String {
    if !name.is_empty() && name.chars().all(is_safe_unquoted_name_char) {
        return name.to_owned();
    }

    format!("'{}'", name.replace('\'', "''"))
}

fn is_safe_unquoted_name_char(c: char) -> bool {
    !(c.is_whitespace() || matches!(c, '(' | ')' | '[' | ']' | ',' | ':' | ';' | '\''))
}

fn serialize_comment(comment: &str) -> String {
    comment.replace(']', "']")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn node(name: Option<&str>, length: Option<f64>, comment: Option<&str>, is_tip: bool) -> Node {
        let mut node = Node::new();
        node.name = name.map(str::to_owned);
        node.parent_edge = length;
        node.comment = comment.map(str::to_owned);
        if !is_tip {
            node.children.push(1);
        }
        node
    }

    #[test]
    fn serializes_safe_and_quoted_names() {
        assert_eq!(serialize_name("A_B-123"), "A_B-123");
        assert_eq!(serialize_name("hungarian dog"), "'hungarian dog'");
        assert_eq!(serialize_name("A,B"), "'A,B'");
        assert_eq!(serialize_name("Darwin's finch"), "'Darwin''s finch'");
        assert_eq!(serialize_name(""), "''");
    }

    #[test]
    fn serializes_all_node_fields() {
        let node = node(Some("A B"), Some(0.1), Some("metadata"), true);
        assert_eq!(
            serialize_node(&node, NewickFormat::AllFields),
            "'A B':0.1[metadata]"
        );
        assert_eq!(serialize_node(&node, NewickFormat::NoComments), "'A B':0.1");
    }

    #[test]
    fn serializes_leaf_and_internal_field_selections() {
        let leaf = node(Some("leaf"), Some(0.1), None, true);
        let internal = node(Some("internal"), Some(0.2), None, false);

        assert_eq!(
            serialize_node(&leaf, NewickFormat::LeafLengthsLeafNames),
            "leaf:0.1"
        );
        assert_eq!(
            serialize_node(&internal, NewickFormat::LeafLengthsLeafNames),
            ""
        );
        assert_eq!(
            serialize_node(&leaf, NewickFormat::InternalLengthsLeafNames),
            "leaf"
        );
        assert_eq!(
            serialize_node(&internal, NewickFormat::InternalLengthsLeafNames),
            ":0.2"
        );
    }

    #[test]
    fn serializes_topology_without_attributes() {
        let node = node(Some("A"), Some(0.1), Some("metadata"), true);
        assert_eq!(serialize_node(&node, NewickFormat::Topology), "");
    }

    #[test]
    fn escapes_closing_brackets_in_comments() {
        let node = node(None, None, Some("left]right"), true);
        assert_eq!(
            serialize_node(&node, NewickFormat::AllFields),
            "[left']right]"
        );
    }
}

#[cfg(test)]
// These tests are lifted from the ETE3 test suite to ensure compatible Newick serialization.
mod tests_ete3_serialization {
    use super::*;

    // These Newick trees come from the ETE3 test suite:
    // https://github.com/etetoolkit/ete/blob/439c7efd8668408fb813e78c3febce1dc627fb43/ete3/test/datasets.py
    const NW_SIMPLE5: &str = "(H,(A,(B,C,D)),D,T,S,(U,Y));";
    const NW_SIMPLE6: &str = "(H,(A,(B,(C),(T))),D);";
    const NW_FULL: &str = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/src/tree/fixtures/ete3/nw_full.nwk"
    ));
    const NW2_FULL: &str = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/src/tree/fixtures/ete3/nw2_full.nwk"
    ));

    #[test]
    fn roundtrips_full_nhx() {
        let tree = Tree::from_newick(NW_FULL).expect("Could not parse NW_FULL");
        assert_eq!(NW_FULL, tree.to_newick().expect("Could not write NW_FULL"));
    }

    #[test]
    fn roundtrips_complex_newick() {
        let tree = Tree::from_newick(NW2_FULL).expect("Could not parse NW2_FULL");
        assert_eq!(
            NW2_FULL,
            tree.to_newick().expect("Could not write NW2_FULL")
        );
    }

    #[test]
    fn roundtrips_weird_topologies() {
        let cases = [(NW_SIMPLE5, "NW_SIMPLE5"), (NW_SIMPLE6, "NW_SIMPLE6")];
        for (newick, name) in cases {
            let tree =
                Tree::from_newick(newick).unwrap_or_else(|_| panic!("Could not parse {name}"));
            assert_eq!(
                newick,
                tree.to_newick()
                    .unwrap_or_else(|_| panic!("Could not write {name}"))
            );
        }
    }

    #[ignore]
    #[test]
    fn roundtrips_single_node_trees() {
        for newick in ["(hola);", "hola;"] {
            let tree = Tree::from_newick(newick).expect("Could not read tree");
            assert_eq!(newick, tree.to_newick().expect("Could not write tree"));
        }
    }

    #[test]
    fn roundtrips_with_features() {
        let newick = "(((A[&&NHX:name=A],B[&&NHX:name=B])[&&NHX:name=NoName],C[&&NHX:name=C])[&&NHX:name=I],(D[&&NHX:name=D],F[&&NHX:name=F])[&&NHX:name=J])[&&NHX:name=root];";
        let tree = Tree::from_newick(newick).expect("Could not parse tree");

        assert_eq!(newick, tree.to_newick().expect("Could not write tree"));
    }
}
