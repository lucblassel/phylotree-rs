use thiserror::Error;

use super::{
    tokenizer::{NewickToken, Span, SpannedToken, Tokenizer, TokenizerError},
    Node, NodeId,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ParserState {
    Child,        // After '(' or ','
    Attribute,    // After a child, ')' or at root
    BranchLength, // After ':'
}

/// Controls how unquoted labels on internal nodes are interpreted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InternalNodeLabelMode {
    /// Preserve every internal-node label as a node name.
    #[default]
    Name,
    /// Interpret numeric internal-node labels as support and preserve other labels as names.
    NumericSupport,
    /// Require every unquoted internal-node label to be a numeric support value.
    ///
    /// Quoted internal-node labels remain names, providing an explicit way to
    /// represent a numeric-looking or nonnumeric name.
    Support,
}

/// Controls whether branch support is extracted from node comments.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum SupportCommentMode {
    /// Preserve comments without interpreting them as support annotations.
    #[default]
    None,
    /// Read the `B` key from comments using the `&&NHX:B=<value>` convention.
    Nhx,
    /// Read the `support` key from comments using the `&support=<value>` convention.
    Beast,
    /// Read an exact key from plain, NHX, or BEAST-style key-value comments.
    Custom(String),
}

impl SupportCommentMode {
    /// Creates a custom comment mode matching `key` exactly.
    ///
    /// Both borrowed string slices and owned strings are accepted.
    pub fn custom(key: impl Into<String>) -> Self {
        Self::Custom(key.into())
    }

    fn support_value<'a>(&self, comment: &'a str) -> Option<&'a str> {
        match self {
            Self::None => None,
            Self::Nhx => comment
                .strip_prefix("&&NHX:")
                .and_then(|body| find_annotation_value(body, "B", &[':'])),
            Self::Beast => comment
                .strip_prefix('&')
                .filter(|body| !body.starts_with('&'))
                .and_then(|body| find_annotation_value(body, "support", &[','])),
            Self::Custom(key) => {
                let body = comment
                    .strip_prefix("&&NHX:")
                    .or_else(|| comment.strip_prefix('&'))
                    .unwrap_or(comment);
                find_annotation_value(body, key, &[':', ','])
            }
        }
    }
}

fn find_annotation_value<'a>(body: &'a str, key: &str, separators: &[char]) -> Option<&'a str> {
    body.split(|c| separators.contains(&c)).find_map(|field| {
        let (field_key, value) = field.split_once('=')?;
        (field_key.trim() == key).then(|| value.trim())
    })
}

/// Options controlling convention-dependent Newick parsing behavior.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct NewickParseOptions {
    /// How labels attached to internal nodes are interpreted.
    pub internal_node_labels: InternalNodeLabelMode,
    /// Whether comments are inspected for branch-support annotations.
    pub support_comments: SupportCommentMode,
}

/// An error encountered while parsing a Newick tree.
#[derive(Debug, Error)]
pub enum NewickParseError {
    /// Tokenization failed before a complete syntax token could be produced.
    #[error("Tokenizer error: {err}")]
    TokenizerError {
        /// The underlying tokenizer error.
        #[source]
        err: TokenizerError,
        /// Byte range associated with the error.
        span: Span,
    },
    /// The input ended before the tree was terminated by a semicolon.
    #[error("Unexpected end of file")]
    UnexpectedEOF,
    /// A valid token appeared where the Newick grammar did not permit it.
    #[error("Unexpected token {token:?}; expected {expected}")]
    UnexpectedToken {
        /// The token encountered by the parser.
        token: NewickToken,
        /// Half-open byte range containing the token.
        span: Span,
        /// Description of the token or construct expected at this position.
        expected: &'static str,
    },
    /// A node attribute appeared when there was no node to receive it.
    #[error("Cannot parse node attributes without a current node")]
    MissingNodeForAttribute {
        /// Half-open byte range containing the attribute.
        span: Span,
    },
    /// A branch-length field was not a valid floating-point number.
    #[error("Invalid branch length '{input}'")]
    InvalidBranchLength {
        /// Text found in the branch-length field.
        input: String,
        /// Half-open byte range containing the invalid value.
        span: Span,
    },
    /// A label or recognized comment annotation contained an invalid support value.
    #[error("Invalid branch support value '{input}'")]
    InvalidSupportValue {
        /// Text found where a support value was expected.
        input: String,
        /// Half-open byte range containing the label or comment.
        span: Span,
    },
    /// Multiple support encodings assigned different values to the same branch.
    #[error("Conflicting branch support values {existing} and {new}")]
    ConflictingSupportValues {
        /// Support value assigned earlier in the node's attributes.
        existing: f64,
        /// Conflicting support value encountered later.
        new: f64,
        /// Half-open byte range containing the conflicting label or comment.
        span: Span,
    },
}

impl From<TokenizerError> for NewickParseError {
    fn from(err: TokenizerError) -> Self {
        Self::TokenizerError {
            span: err.span(),
            err,
        }
    }
}

impl From<std::io::Error> for NewickParseError {
    fn from(err: std::io::Error) -> Self {
        TokenizerError::from(err).into()
    }
}

pub struct NewickParser<'a, T: Tokenizer> {
    tokenizer: &'a mut T,
    nodes: Vec<Node>,
    stack: Vec<NodeId>,
    current_node: Option<NodeId>,
    state: ParserState,
    options: NewickParseOptions,
}

impl<'a, T: Tokenizer> NewickParser<'a, T> {
    pub fn new(tokenizer: &'a mut T) -> Self {
        Self::with_options(tokenizer, NewickParseOptions::default())
    }

    /// Creates a parser using explicit convention-dependent parsing options.
    pub fn with_options(tokenizer: &'a mut T, options: NewickParseOptions) -> Self {
        Self {
            tokenizer,
            nodes: Vec::new(),
            stack: Vec::new(),
            current_node: None,
            state: ParserState::Child,
            options,
        }
    }

    pub fn parse(&mut self) -> Result<Vec<Node>, NewickParseError> {
        loop {
            if self.tokenizer.peek()?.is_none() {
                return Err(NewickParseError::UnexpectedEOF);
            }
            let token = self
                .tokenizer
                .next_token()?
                .ok_or(NewickParseError::UnexpectedEOF)?;

            let finished = match self.state {
                ParserState::Child => self.handle_expect_child(token)?,
                ParserState::Attribute => self.handle_expect_attribute(token)?,
                ParserState::BranchLength => self.handle_expect_branch_length(token)?,
            };

            if finished {
                if let Some(token) = self.tokenizer.next_token()? {
                    return Err(Self::unexpected(token, "end of input"));
                }
                return Ok(std::mem::take(&mut self.nodes));
            }
        }
    }

    fn add_node(&mut self, parent: Option<NodeId>) -> NodeId {
        let id = self.nodes.len();
        let mut node = Node::new();
        node.set_id(id);

        if let Some(parent_id) = parent {
            node.set_parent(parent_id, None);
            node.set_depth(self.nodes[parent_id].get_depth() + 1);
            self.nodes[parent_id].add_child(id, None);
        }

        self.nodes.push(node);
        id
    }

    fn add_anonymous_child(&mut self, token: SpannedToken) -> Result<NodeId, NewickParseError> {
        let parent = self
            .stack
            .last()
            .copied()
            .ok_or_else(|| Self::unexpected(token, "an opening parenthesis"))?;
        Ok(self.add_node(Some(parent)))
    }

    fn add_leaf_or_root(&mut self) -> NodeId {
        self.add_node(self.stack.last().copied())
    }

    fn set_name(&mut self, node: NodeId, name: String) {
        self.nodes[node].set_name(name);
    }

    fn set_unquoted_label(
        &mut self,
        node: NodeId,
        label: String,
        span: Span,
    ) -> Result<(), NewickParseError> {
        if self.nodes[node].is_tip()
            || self.options.internal_node_labels == InternalNodeLabelMode::Name
        {
            self.set_name(node, label);
            return Ok(());
        }

        match label.parse::<f64>() {
            Ok(support) => self.set_support(node, support, span),
            Err(_)
                if self.options.internal_node_labels == InternalNodeLabelMode::NumericSupport =>
            {
                self.set_name(node, label);
                Ok(())
            }
            Err(_) => Err(NewickParseError::InvalidSupportValue { input: label, span }),
        }
    }

    fn set_support(
        &mut self,
        node: NodeId,
        support: f64,
        span: Span,
    ) -> Result<(), NewickParseError> {
        if let Some(existing) = self.nodes[node].support {
            if existing != support {
                return Err(NewickParseError::ConflictingSupportValues {
                    existing,
                    new: support,
                    span,
                });
            }
        } else {
            self.nodes[node].support = Some(support);
        }
        Ok(())
    }

    fn set_branch_length(&mut self, node: NodeId, length: f64) {
        self.nodes[node].parent_edge = Some(length);
        if let Some(parent) = self.nodes[node].parent {
            self.nodes[parent].set_child_edge(&node, Some(length));
        }
    }

    fn add_comment(
        &mut self,
        node: NodeId,
        comment: &str,
        span: Span,
    ) -> Result<(), NewickParseError> {
        if !self.nodes[node].is_tip() {
            if let Some(input) = self.options.support_comments.support_value(comment) {
                let support =
                    input
                        .parse::<f64>()
                        .map_err(|_| NewickParseError::InvalidSupportValue {
                            input: input.to_owned(),
                            span,
                        })?;
                self.set_support(node, support, span)?;
            }
        }

        match &mut self.nodes[node].comment {
            Some(current) => current.push_str(comment),
            None => self.nodes[node].comment = Some(comment.to_owned()),
        }
        Ok(())
    }

    fn unexpected(token: SpannedToken, expected: &'static str) -> NewickParseError {
        NewickParseError::UnexpectedToken {
            token: token.token,
            span: token.span,
            expected,
        }
    }

    fn current_node_for_attribute(&self, span: Span) -> Result<NodeId, NewickParseError> {
        self.current_node
            .ok_or(NewickParseError::MissingNodeForAttribute { span })
    }

    fn handle_expect_child(&mut self, token: SpannedToken) -> Result<bool, NewickParseError> {
        match token.token {
            NewickToken::ParenOpen => {
                let parent = self.stack.last().copied();
                let node = self.add_node(parent);
                self.stack.push(node);
                self.current_node = None;
                Ok(false)
            }
            NewickToken::Comma => {
                self.add_anonymous_child(token)?;
                self.current_node = None;
                Ok(false)
            }
            NewickToken::ParenClose => {
                self.add_anonymous_child(token)?;
                self.current_node = self.stack.pop();
                self.state = ParserState::Attribute;
                Ok(false)
            }
            NewickToken::UnquotedString(name) | NewickToken::QuotedString(name) => {
                let node = self.add_leaf_or_root();
                self.set_name(node, name);
                self.current_node = Some(node);
                self.state = ParserState::Attribute;
                Ok(false)
            }
            NewickToken::Colon => {
                let node = self.add_leaf_or_root();
                self.current_node = Some(node);
                self.state = ParserState::BranchLength;
                Ok(false)
            }
            NewickToken::Comment(comment) => {
                let node = self.add_leaf_or_root();
                self.add_comment(node, &comment, token.span)?;
                self.current_node = Some(node);
                self.state = ParserState::Attribute;
                Ok(false)
            }
            _ => Err(Self::unexpected(
                token,
                "an opening parenthesis, node attribute, or anonymous child",
            )),
        }
    }

    fn handle_expect_attribute(&mut self, token: SpannedToken) -> Result<bool, NewickParseError> {
        let node = self.current_node_for_attribute(token.span)?;
        match &token.token {
            NewickToken::UnquotedString(label) => {
                if self.nodes[node].name.is_some() || self.nodes[node].parent_edge.is_some() {
                    return Err(Self::unexpected(token, "a delimiter or branch length"));
                }
                self.set_unquoted_label(node, label.clone(), token.span)?;
                Ok(false)
            }
            NewickToken::QuotedString(name) => {
                if self.nodes[node].name.is_some() || self.nodes[node].parent_edge.is_some() {
                    return Err(Self::unexpected(token, "a delimiter or branch length"));
                }
                self.set_name(node, name.clone());
                Ok(false)
            }
            NewickToken::Colon => {
                if self.nodes[node].parent_edge.is_some() {
                    return Err(Self::unexpected(token, "a delimiter"));
                }
                self.state = ParserState::BranchLength;
                Ok(false)
            }
            NewickToken::Comment(comment) => {
                self.add_comment(node, comment, token.span)?;
                Ok(false)
            }
            NewickToken::Comma => {
                if self.stack.is_empty() {
                    return Err(Self::unexpected(token, "a terminating semicolon"));
                }
                self.current_node = None;
                self.state = ParserState::Child;
                Ok(false)
            }
            NewickToken::ParenClose => {
                let closed = self
                    .stack
                    .pop()
                    .ok_or_else(|| Self::unexpected(token, "an open subtree"))?;
                self.current_node = Some(closed);
                Ok(false)
            }
            NewickToken::SemiColon if self.stack.is_empty() && self.current_node.is_some() => {
                Ok(true)
            }
            _ => Err(Self::unexpected(token, "a node attribute or delimiter")),
        }
    }

    fn handle_expect_branch_length(
        &mut self,
        token: SpannedToken,
    ) -> Result<bool, NewickParseError> {
        let node = self.current_node_for_attribute(token.span)?;
        match token.token {
            NewickToken::Comment(comment) => {
                self.add_comment(node, &comment, token.span)?;
                Ok(false)
            }
            NewickToken::UnquotedString(input) => {
                let length =
                    input
                        .parse::<f64>()
                        .map_err(|_| NewickParseError::InvalidBranchLength {
                            input,
                            span: token.span,
                        })?;
                self.set_branch_length(node, length);
                self.state = ParserState::Attribute;
                Ok(false)
            }
            _ => Err(Self::unexpected(
                token,
                "a comment or unquoted branch length",
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::tokenizer::NewickTokenizer;

    pub(super) fn parse(newick: &str) -> Result<Vec<Node>, NewickParseError> {
        parse_with_options(newick, NewickParseOptions::default())
    }

    fn parse_with_options(
        newick: &str,
        options: NewickParseOptions,
    ) -> Result<Vec<Node>, NewickParseError> {
        let mut tokenizer = NewickTokenizer::new(newick.as_bytes());
        NewickParser::with_options(&mut tokenizer, options).parse()
    }

    pub(super) fn node_named<'a>(nodes: &'a [Node], name: &str) -> &'a Node {
        nodes
            .iter()
            .find(|node| node.name.as_deref() == Some(name))
            .unwrap_or_else(|| panic!("node named '{name}' was not parsed"))
    }

    fn node_by_id(nodes: &[Node], id: NodeId) -> &Node {
        nodes
            .iter()
            .find(|node| node.id == id)
            .unwrap_or_else(|| panic!("node with id {id} was not parsed"))
    }

    fn assert_is_child(nodes: &[Node], parent: &Node, child: &Node) {
        assert_eq!(child.parent, Some(parent.id));
        assert!(parent.children.contains(&child.id));
        assert_eq!(node_by_id(nodes, child.id).parent, Some(parent.id));
    }

    pub(super) fn assert_node_graph_is_consistent(nodes: &[Node]) {
        assert!(!nodes.is_empty());
        assert_eq!(
            nodes.iter().filter(|node| node.parent.is_none()).count(),
            1,
            "a parsed tree must have exactly one root"
        );

        for (expected_id, node) in nodes.iter().enumerate() {
            assert_eq!(node.id, expected_id);

            if let Some(parent_id) = node.parent {
                let parent = node_by_id(nodes, parent_id);
                assert!(parent.children.contains(&node.id));
                assert_eq!(node.get_depth(), parent.get_depth() + 1);
                assert_eq!(parent.get_child_edge(&node.id), node.parent_edge);
            } else {
                assert_eq!(node.get_depth(), 0);
            }

            for child_id in &node.children {
                assert_eq!(node_by_id(nodes, *child_id).parent, Some(node.id));
            }
        }
    }

    #[test]
    fn parses_a_single_named_node() {
        let nodes = parse("A;").unwrap();

        assert_eq!(nodes.len(), 1);
        assert_eq!(nodes[0].name.as_deref(), Some("A"));
        assert_eq!(nodes[0].parent, None);
        assert!(nodes[0].children.is_empty());
        assert_eq!(nodes[0].parent_edge, None);
        assert_eq!(nodes[0].comment, None);
    }

    #[test]
    fn parses_nested_topology_and_internal_names() {
        let nodes = parse("(A,(B,C)Inner)Root;").unwrap();

        assert_eq!(nodes.len(), 5);
        let root = node_named(&nodes, "Root");
        let inner = node_named(&nodes, "Inner");
        let a = node_named(&nodes, "A");
        let b = node_named(&nodes, "B");
        let c = node_named(&nodes, "C");

        assert_eq!(root.parent, None);
        assert_eq!(root.children.len(), 2);
        assert_is_child(&nodes, root, a);
        assert_is_child(&nodes, root, inner);
        assert_eq!(inner.children.len(), 2);
        assert_is_child(&nodes, inner, b);
        assert_is_child(&nodes, inner, c);
    }

    #[test]
    fn parses_unnamed_nodes() {
        let nodes = parse("(,);").unwrap();

        assert_eq!(nodes.len(), 3);
        let root = nodes.iter().find(|node| node.parent.is_none()).unwrap();
        assert_eq!(root.name, None);
        assert_eq!(root.children.len(), 2);
        for child_id in &root.children {
            let child = node_by_id(&nodes, *child_id);
            assert_eq!(child.name, None);
            assert_eq!(child.parent, Some(root.id));
            assert!(child.children.is_empty());
            assert_eq!(child.get_depth(), 1);
        }
    }

    #[test]
    fn parses_nested_unnamed_topology() {
        let nodes = parse("((,),);").unwrap();

        assert_eq!(nodes.len(), 5);
        let root = node_by_id(&nodes, 0);
        let inner = node_by_id(&nodes, 1);
        assert_eq!(root.parent, None);
        assert_eq!(root.children, vec![inner.id, 4]);
        assert_eq!(inner.parent, Some(root.id));
        assert_eq!(inner.children, vec![2, 3]);
        assert_eq!(inner.get_depth(), 1);
        assert_eq!(node_by_id(&nodes, 2).get_depth(), 2);
        assert_eq!(node_by_id(&nodes, 3).get_depth(), 2);
        assert_eq!(node_by_id(&nodes, 4).get_depth(), 1);
    }

    #[test]
    fn requires_a_semicolon_for_unnamed_topology() {
        assert!(matches!(parse("(,)"), Err(NewickParseError::UnexpectedEOF)));
    }

    #[test]
    fn rejects_unclosed_or_extra_parentheses() {
        assert!(matches!(
            parse("(,;"),
            Err(NewickParseError::UnexpectedToken {
                token: NewickToken::SemiColon,
                span: Span { start: 2, end: 3 },
                ..
            })
        ));
        assert!(matches!(
            parse("(,));"),
            Err(NewickParseError::UnexpectedToken {
                token: NewickToken::ParenClose,
                span: Span { start: 3, end: 4 },
                ..
            })
        ));
    }

    #[test]
    fn rejects_tokens_after_the_semicolon() {
        assert!(matches!(
            parse("(,);(,);"),
            Err(NewickParseError::UnexpectedToken {
                token: NewickToken::ParenOpen,
                span: Span { start: 4, end: 5 },
                ..
            })
        ));
    }

    #[test]
    fn preserves_quoted_and_numeric_looking_names() {
        let nodes = parse("('A B',123,123abc)'Root name';").unwrap();

        for expected_name in ["A B", "123", "123abc", "Root name"] {
            node_named(&nodes, expected_name);
        }
    }

    #[test]
    fn parses_branch_lengths_and_updates_both_ends_of_each_edge() {
        let nodes = parse("(A:0.1,B:-2.5,(C:3e-4)Inner:4)Root;").unwrap();

        let root = node_named(&nodes, "Root");
        let inner = node_named(&nodes, "Inner");
        let expected = [("A", 0.1), ("B", -2.5), ("C", 3e-4), ("Inner", 4.0)];

        for (name, length) in expected {
            let node = node_named(&nodes, name);
            assert_eq!(node.parent_edge, Some(length));
            let parent = node_by_id(&nodes, node.parent.unwrap());
            assert_eq!(parent.get_child_edge(&node.id), Some(length));
        }

        assert_eq!(root.parent_edge, None);
        assert_is_child(&nodes, root, inner);
    }

    #[test]
    fn parses_branch_lengths_on_unnamed_nodes_and_the_root() {
        let nodes = parse("(:0.1,:0.2,(:0.3,:0.4):0.5):0;").unwrap();

        assert_eq!(nodes.len(), 6);
        let root = node_by_id(&nodes, 0);
        let inner = node_by_id(&nodes, 3);
        assert_eq!(root.parent, None);
        assert_eq!(root.parent_edge, Some(0.0));
        assert_eq!(root.children, vec![1, 2, 3]);
        assert_eq!(inner.parent_edge, Some(0.5));
        assert_eq!(inner.children, vec![4, 5]);

        for (id, length) in [(1, 0.1), (2, 0.2), (3, 0.5), (4, 0.3), (5, 0.4)] {
            let node = node_by_id(&nodes, id);
            assert_eq!(node.parent_edge, Some(length));
            let parent = node_by_id(&nodes, node.parent.unwrap());
            assert_eq!(parent.get_child_edge(&id), Some(length));
        }
    }

    // Cases originally covered by tree_impl::tests::read_newick.
    #[test]
    fn parses_newick_validator_cases_from_tree_impl() {
        let cases = [
            ("((D,E)B,(F,G)C)A;", 7),
            ("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;", 6),
            ("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);", 6),
            ("(dog:20,(elephant:30,horse:60):20):50;", 5),
            ("(A,B,(C,D));", 6),
            ("(A,B,(C,D)E)F;", 6),
            (
                "(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0;",
                9,
            ),
            ("(:0.1,:0.2,(:0.3,:0.4):0.5);", 6),
            ("(:0.1,:0.2,(:0.3,:0.4):0.5):0;", 6),
            ("((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;", 6),
            ("(,,(,));", 6),
            (
                "('hungarian dog':20,('indian elephant':30,'swedish horse':60):20):50;",
                5,
            ),
            (
                "('hungarian dog':20[Comment_1],('indian elephant':30,'swedish horse':60[Another interesting comment]):20):50;",
                5,
            ),
        ];

        for (newick, expected_nodes) in cases {
            let nodes = parse(newick).unwrap_or_else(|err| {
                panic!("failed to parse {newick:?}: {err}");
            });
            assert_eq!(nodes.len(), expected_nodes, "wrong node count for {newick}");
            assert_node_graph_is_consistent(&nodes);
        }
    }

    // Cases originally covered by tree_impl::tests::binary_from_newick.
    #[test]
    fn parses_binary_and_polytomous_cases_from_tree_impl() {
        for newick in ["((A,B,C)D,E)F;", "(A,B,(C,D)E)F;", "((D,E)B,(F,G)C)A;"] {
            let nodes = parse(newick).unwrap_or_else(|err| {
                panic!("failed to parse {newick:?}: {err}");
            });
            assert_node_graph_is_consistent(&nodes);
        }
    }

    // Cases originally covered by tree_impl::tests::read_newick_fails.
    #[test]
    fn rejects_malformed_cases_from_tree_impl() {
        for newick in ["((D,E)B,(F,G,C)A;", "((D,E)B,(F,G)C)A"] {
            assert!(parse(newick).is_err(), "unexpectedly parsed {newick:?}");
        }
    }

    #[test]
    fn attaches_comments_to_the_preceding_node() {
        let nodes = parse("(A[first]:0.1,B:0.2[second])Root[root];").unwrap();

        assert_eq!(node_named(&nodes, "A").comment.as_deref(), Some("first"));
        assert_eq!(node_named(&nodes, "B").comment.as_deref(), Some("second"));
        assert_eq!(node_named(&nodes, "Root").comment.as_deref(), Some("root"));
    }

    #[test]
    fn attaches_comments_to_anonymous_nodes() {
        let nodes = parse("([leaf],([left],[right])[inner])[root];").unwrap();

        for (id, expected) in [
            (0, "root"),
            (1, "leaf"),
            (2, "inner"),
            (3, "left"),
            (4, "right"),
        ] {
            assert_eq!(node_by_id(&nodes, id).comment.as_deref(), Some(expected));
        }
        assert_node_graph_is_consistent(&nodes);
    }

    #[test]
    fn concatenates_multiple_comments_on_a_node() {
        let nodes = parse("(A[first][second])Root;").unwrap();

        assert_eq!(
            node_named(&nodes, "A").comment.as_deref(),
            Some("firstsecond")
        );
    }

    #[test]
    fn accepts_comments_before_branch_lengths() {
        let nodes = parse("(A:[metadata]0.1)Root;").unwrap();
        let a = node_named(&nodes, "A");

        assert_eq!(a.comment.as_deref(), Some("metadata"));
        assert_eq!(a.parent_edge, Some(0.1));
    }

    #[test]
    fn rejects_invalid_branch_lengths() {
        assert!(matches!(
            parse("(A:not-a-number)Root;"),
            Err(NewickParseError::InvalidBranchLength {
                input,
                span: Span { start: 3, end: 15 },
            }) if input == "not-a-number"
        ));
    }

    #[test]
    fn rejects_quoted_branch_lengths() {
        assert!(matches!(
            parse("(A:'0.1')Root;"),
            Err(NewickParseError::UnexpectedToken {
                token: NewickToken::QuotedString(value),
                span: Span { start: 3, end: 8 },
                ..
            }) if value == "0.1"
        ));
    }

    #[test]
    fn rejects_empty_input() {
        assert!(matches!(parse(""), Err(NewickParseError::UnexpectedEOF)));
    }

    #[test]
    fn requires_a_terminating_semicolon() {
        assert!(matches!(
            parse("(A,B)Root"),
            Err(NewickParseError::UnexpectedEOF)
        ));
    }

    #[test]
    fn reports_a_missing_node_in_attribute_state() {
        let mut tokenizer = NewickTokenizer::new("Name".as_bytes());
        let mut parser = NewickParser::new(&mut tokenizer);
        parser.state = ParserState::Attribute;

        assert!(matches!(
            parser.parse(),
            Err(NewickParseError::MissingNodeForAttribute {
                span: Span { start: 0, end: 4 }
            })
        ));
    }

    #[test]
    fn reports_a_missing_node_in_branch_length_state() {
        let mut tokenizer = NewickTokenizer::new("0.1".as_bytes());
        let mut parser = NewickParser::new(&mut tokenizer);
        parser.state = ParserState::BranchLength;

        assert!(matches!(
            parser.parse(),
            Err(NewickParseError::MissingNodeForAttribute {
                span: Span { start: 0, end: 3 }
            })
        ));
    }

    #[test]
    fn preserves_internal_numeric_labels_and_comments_by_default() {
        let nodes = parse("(A,B)95[&&NHX:B=90];").unwrap();
        let root = &nodes[0];

        assert_eq!(root.name.as_deref(), Some("95"));
        assert_eq!(root.support, None);
        assert_eq!(root.comment.as_deref(), Some("&&NHX:B=90"));
    }

    #[test]
    fn interprets_only_numeric_internal_labels_in_numeric_support_mode() {
        let options = NewickParseOptions {
            internal_node_labels: InternalNodeLabelMode::NumericSupport,
            ..Default::default()
        };
        let nodes = parse_with_options("((1,2)95,C)Root;", options).unwrap();
        let inner = node_by_id(&nodes, 1);

        assert_eq!(inner.name, None);
        assert_eq!(inner.support, Some(95.0));
        assert_eq!(node_by_id(&nodes, 2).name.as_deref(), Some("1"));
        assert_eq!(node_by_id(&nodes, 3).name.as_deref(), Some("2"));
        assert_eq!(node_by_id(&nodes, 0).name.as_deref(), Some("Root"));
    }

    #[test]
    fn strict_support_mode_rejects_nonnumeric_unquoted_internal_labels() {
        let options = NewickParseOptions {
            internal_node_labels: InternalNodeLabelMode::Support,
            ..Default::default()
        };

        assert!(matches!(
            parse_with_options("(A,B)Clade;", options),
            Err(NewickParseError::InvalidSupportValue {
                input,
                span: Span { start: 5, end: 10 },
            }) if input == "Clade"
        ));
    }

    #[test]
    fn strict_support_mode_preserves_explicitly_quoted_internal_names() {
        let options = NewickParseOptions {
            internal_node_labels: InternalNodeLabelMode::Support,
            ..Default::default()
        };
        let nodes = parse_with_options("(A,B)'95';", options).unwrap();

        assert_eq!(nodes[0].name.as_deref(), Some("95"));
        assert_eq!(nodes[0].support, None);
    }

    #[test]
    fn extracts_nhx_support_and_preserves_the_comment() {
        let options = NewickParseOptions {
            support_comments: SupportCommentMode::Nhx,
            ..Default::default()
        };
        let nodes = parse_with_options("(A,B)[&&NHX:S=human:B=98.5];", options).unwrap();

        assert_eq!(nodes[0].support, Some(98.5));
        assert_eq!(nodes[0].comment.as_deref(), Some("&&NHX:S=human:B=98.5"));
    }

    #[test]
    fn extracts_beast_support_and_preserves_the_comment() {
        let options = NewickParseOptions {
            support_comments: SupportCommentMode::Beast,
            ..Default::default()
        };
        let nodes = parse_with_options("(A,B)[&height=1.2,support=0.97];", options).unwrap();

        assert_eq!(nodes[0].support, Some(0.97));
        assert_eq!(
            nodes[0].comment.as_deref(),
            Some("&height=1.2,support=0.97")
        );
    }

    #[test]
    fn extracts_support_using_a_custom_literal_key() {
        let borrowed_key = NewickParseOptions {
            support_comments: SupportCommentMode::custom("bootstrap"),
            ..Default::default()
        };
        let owned_key = NewickParseOptions {
            support_comments: SupportCommentMode::custom(String::from("bootstrap")),
            ..Default::default()
        };

        let plain = parse_with_options("(A,B)[bootstrap=87];", borrowed_key).unwrap();
        let prefixed = parse_with_options("(A,B)[&bootstrap=88];", owned_key).unwrap();
        assert_eq!(plain[0].support, Some(87.0));
        assert_eq!(prefixed[0].support, Some(88.0));
    }

    #[test]
    fn ignores_support_annotations_on_leaf_nodes() {
        let options = NewickParseOptions {
            support_comments: SupportCommentMode::Nhx,
            ..Default::default()
        };
        let nodes = parse_with_options("(A[&&NHX:B=100],B)Root;", options).unwrap();

        let a = node_named(&nodes, "A");
        assert_eq!(a.support, None);
        assert_eq!(a.comment.as_deref(), Some("&&NHX:B=100"));
    }

    #[test]
    fn rejects_invalid_or_conflicting_support_values() {
        let comments = NewickParseOptions {
            support_comments: SupportCommentMode::Nhx,
            ..Default::default()
        };
        assert!(matches!(
            parse_with_options("(A,B)[&&NHX:B=high];", comments),
            Err(NewickParseError::InvalidSupportValue { input, .. }) if input == "high"
        ));

        let combined = NewickParseOptions {
            internal_node_labels: InternalNodeLabelMode::NumericSupport,
            support_comments: SupportCommentMode::Nhx,
        };
        assert!(matches!(
            parse_with_options("(A,B)95[&&NHX:B=90];", combined),
            Err(NewickParseError::ConflictingSupportValues {
                existing: 95.0,
                new: 90.0,
                ..
            })
        ));
    }

    #[test]
    fn accepts_matching_support_values_from_multiple_encodings() {
        let options = NewickParseOptions {
            internal_node_labels: InternalNodeLabelMode::NumericSupport,
            support_comments: SupportCommentMode::Nhx,
        };
        let nodes = parse_with_options("(A,B)95[&&NHX:B=95.0];", options).unwrap();

        assert_eq!(nodes[0].support, Some(95.0));
    }

    #[test]
    fn propagates_tokenizer_errors_with_their_spans() {
        assert!(matches!(
            parse("'unterminated"),
            Err(NewickParseError::TokenizerError {
                err: TokenizerError::UnterminatedQuoteString { .. },
                span: Span { start: 0, end: 13 },
            })
        ));
    }
}

#[cfg(test)]
mod tests_ete3_parsing {
    use super::tests::{assert_node_graph_is_consistent, node_named, parse};

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
    fn parses_full_nhx() {
        let nodes = parse(NW_FULL).expect("could not parse NW_FULL");

        assert_node_graph_is_consistent(&nodes);
        assert_eq!(
            node_named(&nodes, "Hsa0010730").comment.as_deref(),
            Some("&&NHX:flag=Red")
        );
        assert!(nodes
            .iter()
            .any(|node| node.comment.as_deref() == Some("&&NHX:flag=Black")));
        assert!(nodes
            .iter()
            .any(|node| node.comment.as_deref() == Some("&&NHX:flag=White")));
    }

    #[test]
    fn parses_complex_newick() {
        let nodes = parse(NW2_FULL).expect("could not parse NW2_FULL");

        assert_node_graph_is_consistent(&nodes);
        assert!(nodes.len() > 500);
        node_named(&nodes, "YGR138C");
        node_named(&nodes, "YPR156C");
    }

    #[test]
    fn parses_weird_topologies() {
        for (newick, name) in [(NW_SIMPLE5, "NW_SIMPLE5"), (NW_SIMPLE6, "NW_SIMPLE6")] {
            let nodes = parse(newick).unwrap_or_else(|err| panic!("could not parse {name}: {err}"));
            assert_node_graph_is_consistent(&nodes);
        }
    }

    #[test]
    fn parses_single_node_trees() {
        let parenthesized = parse("(hola);").expect("could not parse parenthesized node");
        assert_node_graph_is_consistent(&parenthesized);
        assert_eq!(parenthesized.len(), 2);
        assert_eq!(node_named(&parenthesized, "hola").parent, Some(0));

        let bare = parse("hola;").expect("could not parse bare node");
        assert_node_graph_is_consistent(&bare);
        assert_eq!(bare.len(), 1);
        assert_eq!(node_named(&bare, "hola").parent, None);
    }

    #[test]
    fn parses_nhx_features_as_comments() {
        let newick = "(((A[&&NHX:name=A],B[&&NHX:name=B])[&&NHX:name=NoName],C[&&NHX:name=C])[&&NHX:name=I],(D[&&NHX:name=D],F[&&NHX:name=F])[&&NHX:name=J])[&&NHX:name=root];";
        let nodes = parse(newick).expect("could not parse NHX features");

        assert_node_graph_is_consistent(&nodes);
        assert_eq!(
            node_named(&nodes, "A").comment.as_deref(),
            Some("&&NHX:name=A")
        );
        assert!(nodes
            .iter()
            .any(|node| node.comment.as_deref() == Some("&&NHX:name=NoName")));
        assert!(nodes
            .iter()
            .any(|node| node.comment.as_deref() == Some("&&NHX:name=root")));
    }
}
