use ariadne::{Color, Config, Label, Report, ReportKind, Source};
use phylotree::tree::{NewickParseError, TokenizerError};
use std::io::{self, Write};

pub(crate) fn write_newick_error(
    writer: impl Write,
    source_name: &str,
    source: &str,
    error: &NewickParseError,
    color: bool,
) -> io::Result<()> {
    let span = error
        .span()
        .map(|span| span.range())
        .unwrap_or(source.len()..source.len());

    Report::build(ReportKind::Error, (source_name, span.clone()))
        .with_config(Config::default().with_color(color))
        .with_message(error.to_string())
        .with_label(
            Label::new((source_name, span))
                .with_color(Color::Red)
                .with_message(label_message(error)),
        )
        .finish()
        .write((source_name, Source::from(source)), writer)
}

fn label_message(error: &NewickParseError) -> String {
    match error {
        NewickParseError::TokenizerError { err, .. } => match err {
            TokenizerError::IOError { .. } => "the input stream could not be read".into(),
            TokenizerError::UnterminatedQuoteString { .. } => {
                "this quoted label is missing its closing quote".into()
            }
            TokenizerError::UnterminatedComment { .. } => {
                "this comment is missing its closing bracket".into()
            }
            TokenizerError::InvalidCharacter { char, .. } => {
                format!("'{char}' is not valid at this position")
            }
        },
        NewickParseError::UnexpectedEOF => "expected more input before the end of the tree".into(),
        NewickParseError::UnexpectedToken { expected, .. } => {
            format!("expected {expected}")
        }
        NewickParseError::MissingNodeForAttribute { .. } => {
            "this attribute has no preceding node".into()
        }
        NewickParseError::InvalidBranchLength { .. } => {
            "expected a floating-point branch length".into()
        }
        NewickParseError::InvalidSupportValue { .. } => {
            "expected a numeric branch support value".into()
        }
        NewickParseError::ConflictingSupportValues { existing, .. } => {
            format!("this conflicts with the previous support value {existing}")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use phylotree::tree::Tree;

    fn render(source: &str) -> String {
        let error = Tree::from_newick(source).unwrap_err();
        let mut output = Vec::new();
        write_newick_error(&mut output, "tree.nwk", source, &error, false).unwrap();
        String::from_utf8(output).unwrap()
    }

    #[test]
    fn renders_a_spanned_parse_error() {
        let rendered = render("(A:not-a-number,B);");

        assert!(rendered.contains("tree.nwk"));
        assert!(rendered.contains("not-a-number"));
        assert!(rendered.contains("expected a floating-point branch length"));
    }

    #[test]
    fn renders_unexpected_eof_at_the_end_of_the_source() {
        let rendered = render("(A,B)");

        assert!(rendered.contains("Unexpected end of file"));
        assert!(rendered.contains("expected more input before the end of the tree"));
    }
}
