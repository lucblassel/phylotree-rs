use super::{Node, NodeId, Tree, TreeError};

/// Selects which node names are included in Newick output.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub enum NewickNameFormat {
    /// Omit all node names.
    None,
    /// Output names only for leaf nodes.
    Leaves,
    /// Output names only for internal nodes.
    Internal,
    /// Output names for all nodes.
    #[default]
    All,
}

impl NewickNameFormat {
    fn includes(self, is_tip: bool) -> bool {
        match self {
            Self::None => false,
            Self::Leaves => is_tip,
            Self::Internal => !is_tip,
            Self::All => true,
        }
    }
}

/// Selects which branch lengths are included in Newick output.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub enum NewickBranchLengthFormat {
    /// Omit all branch lengths.
    None,
    /// Output branch lengths only for leaf nodes.
    Leaves,
    /// Output branch lengths only for internal nodes.
    Internal,
    /// Output branch lengths for all nodes.
    #[default]
    All,
}

impl NewickBranchLengthFormat {
    fn includes(self, is_tip: bool) -> bool {
        match self {
            Self::None => false,
            Self::Leaves => is_tip,
            Self::Internal => !is_tip,
            Self::All => true,
        }
    }
}

/// Selects which raw node comments are included in Newick output.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub enum NewickCommentFormat {
    /// Omit all raw comments.
    None,
    /// Output comments only for leaf nodes.
    Leaves,
    /// Output comments only for internal nodes.
    Internal,
    /// Output comments for all nodes.
    #[default]
    All,
}

impl NewickCommentFormat {
    fn includes(self, is_tip: bool) -> bool {
        match self {
            Self::None => false,
            Self::Leaves => is_tip,
            Self::Internal => !is_tip,
            Self::All => true,
        }
    }
}

/// Selects how support values on internal branches are represented.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum NewickSupportFormat {
    /// Omit dedicated support values.
    #[default]
    None,
    /// Serialize support in the internal-node label position, such as `(A,B)95`.
    InternalNodeLabel,
    /// Serialize support as an NHX `B` annotation, such as `[&&NHX:B=95]`.
    NhxComment,
    /// Serialize support as a BEAST annotation, such as `[&support=0.95]`.
    BeastComment,
    /// Serialize support as a plain comment containing an exact custom key.
    CustomComment(String),
}

impl NewickSupportFormat {
    /// Creates a custom comment format using `key` in `[key=<support>]`.
    ///
    /// Both borrowed string slices and owned strings are accepted.
    pub fn custom_comment(key: impl Into<String>) -> Self {
        Self::CustomComment(key.into())
    }
}

/// Independent options controlling each supported Newick node attribute.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NewickSerializeOptions {
    /// Which node names to output.
    pub names: NewickNameFormat,
    /// Which branch lengths to output.
    pub branch_lengths: NewickBranchLengthFormat,
    /// Which raw node comments to output.
    pub comments: NewickCommentFormat,
    /// How to output support values attached to internal branches.
    pub support: NewickSupportFormat,
}

impl Default for NewickSerializeOptions {
    fn default() -> Self {
        Self {
            names: NewickNameFormat::All,
            branch_lengths: NewickBranchLengthFormat::All,
            comments: NewickCommentFormat::All,
            support: NewickSupportFormat::None,
        }
    }
}

/// A compatibility preset for commonly used combinations of serialization options.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NewickFormat {
    /// Output all names, branch lengths, and comments.
    ///
    AllFields,
    /// Output only the tree topology.
    Topology,
    /// Output all names and branch lengths, but no comments.
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
    /// Output attributes using an explicit set of independent options.
    Custom(NewickSerializeOptions),
}

impl NewickFormat {
    /// Expands this compatibility preset into independent serialization options.
    pub fn options(&self) -> NewickSerializeOptions {
        let support = NewickSupportFormat::None;
        match self {
            Self::AllFields => NewickSerializeOptions {
                names: NewickNameFormat::All,
                branch_lengths: NewickBranchLengthFormat::All,
                comments: NewickCommentFormat::All,
                support,
            },
            Self::Topology => NewickSerializeOptions {
                names: NewickNameFormat::None,
                branch_lengths: NewickBranchLengthFormat::None,
                comments: NewickCommentFormat::None,
                support,
            },
            Self::NoComments => NewickSerializeOptions {
                names: NewickNameFormat::All,
                branch_lengths: NewickBranchLengthFormat::All,
                comments: NewickCommentFormat::None,
                support,
            },
            Self::OnlyNames => NewickSerializeOptions {
                names: NewickNameFormat::All,
                branch_lengths: NewickBranchLengthFormat::None,
                comments: NewickCommentFormat::None,
                support,
            },
            Self::OnlyLengths => NewickSerializeOptions {
                names: NewickNameFormat::None,
                branch_lengths: NewickBranchLengthFormat::All,
                comments: NewickCommentFormat::None,
                support,
            },
            Self::LeafLengthsAllNames => NewickSerializeOptions {
                names: NewickNameFormat::All,
                branch_lengths: NewickBranchLengthFormat::Leaves,
                comments: NewickCommentFormat::None,
                support,
            },
            Self::LeafLengthsLeafNames => NewickSerializeOptions {
                names: NewickNameFormat::Leaves,
                branch_lengths: NewickBranchLengthFormat::Leaves,
                comments: NewickCommentFormat::None,
                support,
            },
            Self::InternalLengthsLeafNames => NewickSerializeOptions {
                names: NewickNameFormat::Leaves,
                branch_lengths: NewickBranchLengthFormat::Internal,
                comments: NewickCommentFormat::None,
                support,
            },
            Self::AllLengthsLeafNames => NewickSerializeOptions {
                names: NewickNameFormat::Leaves,
                branch_lengths: NewickBranchLengthFormat::All,
                comments: NewickCommentFormat::None,
                support,
            },
            Self::Custom(options) => options.clone(),
        }
    }
}

impl From<NewickFormat> for NewickSerializeOptions {
    fn from(format: NewickFormat) -> Self {
        format.options()
    }
}

/// Serializes a tree to Newick using independent output options.
pub(crate) struct NewickSerializer<'a> {
    tree: &'a Tree,
    options: NewickSerializeOptions,
}

impl<'a> NewickSerializer<'a> {
    /// Creates a serializer for `tree` using `options`.
    pub(crate) fn new(tree: &'a Tree, options: NewickSerializeOptions) -> Self {
        Self { tree, options }
    }

    /// Serializes the complete tree, including the terminating semicolon.
    pub(crate) fn serialize(&self) -> Result<String, TreeError> {
        self.validate()?;
        let root = self.tree.get_root()?;
        let mut output = self.serialize_subtree(&root)?;
        output.push(';');
        Ok(output)
    }

    fn validate(&self) -> Result<(), TreeError> {
        if let NewickSupportFormat::CustomComment(key) = &self.options.support {
            if !is_valid_comment_key(key) {
                return Err(TreeError::InvalidNewickSupportKey(key.clone()));
            }
        }
        Ok(())
    }

    fn serialize_subtree(&self, node_id: &NodeId) -> Result<String, TreeError> {
        let node = self.tree.get(node_id)?;
        validate_node(node, &self.options)?;
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

        output.push_str(&serialize_node_fields(node, &self.options));
        Ok(output)
    }
}

fn validate_node(node: &Node, options: &NewickSerializeOptions) -> Result<(), TreeError> {
    if node.is_tip() || node.support.is_none() {
        return Ok(());
    }

    if options.support == NewickSupportFormat::InternalNodeLabel
        && options.names.includes(false)
        && node.name.is_some()
    {
        return Err(TreeError::ConflictingNewickNodeLabel(node.id));
    }

    if options.comments.includes(false) {
        if let Some(encoded) = node
            .comment
            .as_deref()
            .and_then(|comment| support_value_from_comment(comment, &options.support))
        {
            let encoded = encoded
                .parse::<f64>()
                .map_err(|_| TreeError::InvalidNewickSupportComment(node.id))?;
            if Some(encoded) != node.support {
                return Err(TreeError::ConflictingNewickSupport(node.id));
            }
        }
    }

    Ok(())
}

/// Serializes the fields belonging to one node with a compatibility preset.
pub(crate) fn serialize_node(node: &Node, format: NewickFormat) -> Result<String, TreeError> {
    serialize_node_with_options(node, format.options())
}

/// Serializes the fields belonging to one node using independent options.
pub(crate) fn serialize_node_with_options(
    node: &Node,
    options: NewickSerializeOptions,
) -> Result<String, TreeError> {
    if let NewickSupportFormat::CustomComment(key) = &options.support {
        if !is_valid_comment_key(key) {
            return Err(TreeError::InvalidNewickSupportKey(key.clone()));
        }
    }
    validate_node(node, &options)?;
    Ok(serialize_node_fields(node, &options))
}

fn serialize_node_fields(node: &Node, options: &NewickSerializeOptions) -> String {
    let mut output = String::new();
    let is_tip = node.is_tip();

    if options.names.includes(is_tip) {
        if let Some(name) = &node.name {
            output.push_str(&serialize_name(name));
        }
    }

    if !is_tip && options.support == NewickSupportFormat::InternalNodeLabel {
        if let Some(support) = node.support {
            output.push_str(&support.to_string());
        }
    }

    if options.branch_lengths.includes(is_tip) {
        if let Some(length) = node.parent_edge {
            output.push(':');
            output.push_str(&length.to_string());
        }
    }

    if options.comments.includes(is_tip) {
        if let Some(comment) = &node.comment {
            output.push('[');
            output.push_str(&serialize_comment(comment));
            output.push(']');
        }
    }

    if !is_tip {
        if let Some(support) = node.support {
            if is_comment_support_format(&options.support)
                && !(options.comments.includes(is_tip)
                    && comment_already_contains_support(node, &options.support))
            {
                output.push('[');
                output.push_str(&serialize_support_comment(support, &options.support));
                output.push(']');
            }
        }
    }

    output
}

fn is_comment_support_format(format: &NewickSupportFormat) -> bool {
    matches!(
        format,
        NewickSupportFormat::NhxComment
            | NewickSupportFormat::BeastComment
            | NewickSupportFormat::CustomComment(_)
    )
}

fn comment_already_contains_support(node: &Node, format: &NewickSupportFormat) -> bool {
    let Some(comment) = node.comment.as_deref() else {
        return false;
    };
    let Some(value) = support_value_from_comment(comment, format) else {
        return false;
    };
    value.parse::<f64>().ok() == node.support
}

fn support_value_from_comment<'a>(
    comment: &'a str,
    format: &NewickSupportFormat,
) -> Option<&'a str> {
    match format {
        NewickSupportFormat::None | NewickSupportFormat::InternalNodeLabel => None,
        NewickSupportFormat::NhxComment => comment
            .strip_prefix("&&NHX:")
            .and_then(|body| find_annotation_value(body, "B", &[':'])),
        NewickSupportFormat::BeastComment => comment
            .strip_prefix('&')
            .filter(|body| !body.starts_with('&'))
            .and_then(|body| find_annotation_value(body, "support", &[','])),
        NewickSupportFormat::CustomComment(key) => {
            let body = comment
                .strip_prefix("&&NHX:")
                .or_else(|| comment.strip_prefix('&'))
                .unwrap_or(comment);
            find_annotation_value(body, key, &[':', ','])
        }
    }
}

fn find_annotation_value<'a>(body: &'a str, key: &str, separators: &[char]) -> Option<&'a str> {
    body.split(|c| separators.contains(&c)).find_map(|field| {
        let (field_key, value) = field.split_once('=')?;
        (field_key.trim() == key).then(|| value.trim())
    })
}

fn serialize_support_comment(support: f64, format: &NewickSupportFormat) -> String {
    match format {
        NewickSupportFormat::NhxComment => format!("&&NHX:B={support}"),
        NewickSupportFormat::BeastComment => format!("&support={support}"),
        NewickSupportFormat::CustomComment(key) => format!("{key}={support}"),
        NewickSupportFormat::None | NewickSupportFormat::InternalNodeLabel => String::new(),
    }
}

fn is_valid_comment_key(key: &str) -> bool {
    !key.is_empty()
        && key
            .chars()
            .all(|c| !c.is_whitespace() && !matches!(c, '[' | ']' | ':' | ',' | '='))
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
    use crate::tree::{InternalNodeLabelMode, NewickParseOptions, SupportCommentMode};

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
    fn compatibility_formats_expand_to_the_expected_options() {
        assert_eq!(
            NewickFormat::Topology.options(),
            NewickSerializeOptions {
                names: NewickNameFormat::None,
                branch_lengths: NewickBranchLengthFormat::None,
                comments: NewickCommentFormat::None,
                support: NewickSupportFormat::None,
            }
        );
        assert_eq!(
            NewickFormat::InternalLengthsLeafNames.options(),
            NewickSerializeOptions {
                names: NewickNameFormat::Leaves,
                branch_lengths: NewickBranchLengthFormat::Internal,
                comments: NewickCommentFormat::None,
                support: NewickSupportFormat::None,
            }
        );

        let custom = NewickSerializeOptions {
            support: NewickSupportFormat::NhxComment,
            ..Default::default()
        };
        assert_eq!(NewickFormat::Custom(custom.clone()).options(), custom);
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
    fn serializes_all_legacy_node_fields() {
        let node = node(Some("A B"), Some(0.1), Some("metadata"), true);
        assert_eq!(
            serialize_node(&node, NewickFormat::AllFields).unwrap(),
            "'A B':0.1[metadata]"
        );
        assert_eq!(
            serialize_node(&node, NewickFormat::NoComments).unwrap(),
            "'A B':0.1"
        );
    }

    #[test]
    fn serializes_independent_field_selections() {
        let leaf = node(Some("leaf"), Some(0.1), Some("leaf-comment"), true);
        let internal = node(Some("internal"), Some(0.2), Some("internal-comment"), false);
        let options = NewickSerializeOptions {
            names: NewickNameFormat::Leaves,
            branch_lengths: NewickBranchLengthFormat::Internal,
            comments: NewickCommentFormat::Leaves,
            support: NewickSupportFormat::None,
        };

        assert_eq!(
            serialize_node_with_options(&leaf, options.clone()).unwrap(),
            "leaf[leaf-comment]"
        );
        assert_eq!(
            serialize_node_with_options(&internal, options).unwrap(),
            ":0.2"
        );
    }

    #[test]
    fn serializes_support_in_each_supported_representation() {
        let mut internal = node(None, Some(0.2), None, false);
        internal.support = Some(95.0);

        for (support, expected) in [
            (NewickSupportFormat::InternalNodeLabel, "95:0.2"),
            (NewickSupportFormat::NhxComment, ":0.2[&&NHX:B=95]"),
            (NewickSupportFormat::BeastComment, ":0.2[&support=95]"),
            (
                NewickSupportFormat::custom_comment("bootstrap"),
                ":0.2[bootstrap=95]",
            ),
        ] {
            let options = NewickSerializeOptions {
                support,
                ..Default::default()
            };
            assert_eq!(
                serialize_node_with_options(&internal, options).unwrap(),
                expected
            );
        }
    }

    #[test]
    fn rejects_name_and_support_competing_for_the_internal_label() {
        let mut internal = node(Some("Clade"), None, None, false);
        internal.support = Some(95.0);
        let options = NewickSerializeOptions {
            support: NewickSupportFormat::InternalNodeLabel,
            ..Default::default()
        };

        assert!(matches!(
            serialize_node_with_options(&internal, options),
            Err(TreeError::ConflictingNewickNodeLabel(0))
        ));
    }

    #[test]
    fn validates_custom_support_keys() {
        let mut internal = node(None, None, None, false);
        internal.support = Some(95.0);
        let options = NewickSerializeOptions {
            support: NewickSupportFormat::custom_comment("invalid:key"),
            ..Default::default()
        };

        assert!(matches!(
            serialize_node_with_options(&internal, options),
            Err(TreeError::InvalidNewickSupportKey(key)) if key == "invalid:key"
        ));
    }

    #[test]
    fn avoids_duplicating_preserved_support_comments() {
        let mut internal = node(None, None, Some("&&NHX:B=95"), false);
        internal.support = Some(95.0);
        let options = NewickSerializeOptions {
            support: NewickSupportFormat::NhxComment,
            ..Default::default()
        };

        assert_eq!(
            serialize_node_with_options(&internal, options).unwrap(),
            "[&&NHX:B=95]"
        );
    }

    #[test]
    fn roundtrips_internal_label_support_with_matching_options() {
        let parse_options = NewickParseOptions {
            internal_node_labels: InternalNodeLabelMode::NumericSupport,
            ..Default::default()
        };
        let tree = Tree::from_newick_with_options("((A,B)95,C)99;", parse_options).unwrap();
        let serialize_options = NewickSerializeOptions {
            support: NewickSupportFormat::InternalNodeLabel,
            ..Default::default()
        };

        assert_eq!(
            tree.to_formatted_newick(NewickFormat::Custom(serialize_options))
                .unwrap(),
            "((A,B)95,C)99;"
        );
    }

    #[test]
    fn roundtrips_comment_support_without_duplication() {
        let cases = [
            (
                SupportCommentMode::Nhx,
                NewickSupportFormat::NhxComment,
                "(A,B)[&&NHX:B=95];",
            ),
            (
                SupportCommentMode::Beast,
                NewickSupportFormat::BeastComment,
                "(A,B)[&support=0.95];",
            ),
            (
                SupportCommentMode::custom("bootstrap"),
                NewickSupportFormat::custom_comment("bootstrap"),
                "(A,B)[bootstrap=95];",
            ),
        ];

        for (parse_support, serialize_support, newick) in cases {
            let parse_options = NewickParseOptions {
                support_comments: parse_support,
                ..Default::default()
            };
            let tree = Tree::from_newick_with_options(newick, parse_options).unwrap();
            let serialize_options = NewickSerializeOptions {
                support: serialize_support,
                ..Default::default()
            };
            assert_eq!(
                tree.to_formatted_newick(NewickFormat::Custom(serialize_options))
                    .unwrap(),
                newick
            );
        }
    }

    #[test]
    fn serializes_topology_without_attributes() {
        let node = node(Some("A"), Some(0.1), Some("metadata"), true);
        assert_eq!(serialize_node(&node, NewickFormat::Topology).unwrap(), "");
    }

    #[test]
    fn escapes_closing_brackets_in_comments() {
        let node = node(None, None, Some("left]right"), true);
        assert_eq!(
            serialize_node(&node, NewickFormat::AllFields).unwrap(),
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
