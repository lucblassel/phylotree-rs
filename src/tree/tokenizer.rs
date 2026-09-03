use std::io::{BufReader, Read};
use thiserror::Error;
use utf8_chars::BufReadCharsExt;

/// A lexical token recognized in a Newick/NHX input stream.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum NewickToken {
    /// `(`
    ParenOpen,
    /// `)`
    ParenClose,
    /// `,`
    Comma,
    /// `:`
    Colon,
    /// `;`
    SemiColon,
    /// A single-quoted name, with doubled quotes decoded as one quote.
    QuotedString(String),
    /// An atom that contains neither whitespace nor Newick delimiters.
    ///
    /// Numeric-looking atoms are intentionally left as strings so the parser
    /// can interpret them according to their grammatical context.
    UnquotedString(String),
    /// The contents of a bracket-delimited comment, excluding the brackets.
    Comment(String),
}

/// A half-open byte range (`start..end`) in the original UTF-8 input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Span {
    /// Byte offset of the first byte in the range.
    pub start: usize,
    /// Byte offset immediately after the last byte in the range.
    pub end: usize,
}

/// A token together with its location in the original input.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SpannedToken {
    /// The token parsed from the input.
    pub token: NewickToken,
    /// The byte range containing the complete token, including delimiters.
    pub span: Span,
}

/// An error encountered while reading or tokenizing Newick/NHX input.
#[derive(Debug, Error)]
pub enum TokenizerError {
    #[error("IO Error: {source}")]
    IOError { source: std::io::Error },
    #[error("Unterminated quoted string")]
    UnterminatedQuoteString { span: Span },
    #[error("Unterminated newick comment")]
    UnterminatedComment { span: Span },
    #[error("Invalid character: '{char}'")]
    InvalidCharacter { char: char, span: Span },
}

impl TokenizerError {
    /// Returns the byte range associated with this error.
    ///
    /// I/O errors are not tied to a token and therefore return `0..0`.
    pub fn span(&self) -> Span {
        match self {
            TokenizerError::IOError { .. } => Span { start: 0, end: 0 },
            TokenizerError::UnterminatedQuoteString { span } => *span,
            TokenizerError::UnterminatedComment { span } => *span,
            TokenizerError::InvalidCharacter { span, .. } => *span,
        }
    }
}

impl From<std::io::Error> for TokenizerError {
    fn from(err: std::io::Error) -> Self {
        Self::IOError { source: err }
    }
}

/// A token source consumed by the Newick parser.
pub trait Tokenizer {
    /// Returns and consumes the next token, or `None` at the end of the input.
    fn next_token(&mut self) -> Result<Option<SpannedToken>, TokenizerError>;

    /// Returns the next token without consuming it.
    ///
    /// Repeated calls return the same token until [`Tokenizer::next_token`] is
    /// called.
    fn peek(&mut self) -> Result<Option<&SpannedToken>, TokenizerError>;
}

/// A streaming Newick/NHX tokenizer over any [`Read`] implementation.
///
/// Input is decoded incrementally as UTF-8 and token spans are reported as
/// half-open byte ranges into the original stream.
struct NewickTokenizer<R: Read> {
    reader: BufReader<R>,
    byte_offset: usize,
    peeked_token: Option<SpannedToken>,
    current_char: Option<char>,
    current_start: usize,
}

// Characters that are not allowed in unquoted strings
// in addition to whitespace.
fn is_special_char(c: char) -> bool {
    matches!(c, '(' | ')' | '[' | ']' | ',' | ':' | ';')
}

// Character that is allowd in an unquoted string
fn is_safe_char(c: char) -> bool {
    !(c.is_whitespace() || is_special_char(c))
}

impl<R: Read> NewickTokenizer<R> {
    /// Creates a tokenizer that reads incrementally from `reader`.
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
            byte_offset: 0,
            peeked_token: None,
            current_char: None,
            current_start: 0,
        }
    }

    fn get_current_span(&self) -> Span {
        Span {
            start: self.current_start,
            end: self.byte_offset,
        }
    }

    // Peek at the next character without advancing the consumed byte offset.
    fn peek_char(&mut self) -> std::io::Result<Option<char>> {
        if self.current_char.is_none() {
            self.current_char = self.reader.read_char()?;
        }

        Ok(self.current_char)
    }

    // Consume the next character and advance the byte offset by its UTF-8 width.
    fn consume_char(&mut self) -> std::io::Result<Option<char>> {
        let c = match self.current_char.take() {
            Some(c) => Some(c),
            None => self.reader.read_char()?,
        };
        if let Some(c) = c {
            self.byte_offset += c.len_utf8();
        }
        Ok(c)
    }

    // Skip whitespace between newick tokens
    fn skip_whitespace(&mut self) -> std::io::Result<()> {
        while let Some(c) = self.peek_char()? {
            if c.is_whitespace() {
                self.consume_char()?;
            } else {
                break;
            }
        }
        Ok(())
    }

    fn scan_unquoted_string(&mut self, c: char) -> Result<NewickToken, TokenizerError> {
        let mut s = String::from(c);
        loop {
            match self.peek_char()? {
                Some(c) if is_safe_char(c) => {
                    let c = self.consume_char()?.unwrap();
                    s.push(c);
                }
                _ => break,
            }
        }

        Ok(NewickToken::UnquotedString(s))
    }

    fn scan_quoted_string(&mut self) -> Result<NewickToken, TokenizerError> {
        let mut s = String::new();
        loop {
            let c = match self.consume_char()? {
                Some(c) => c,
                None => {
                    return Err(TokenizerError::UnterminatedQuoteString {
                        span: self.get_current_span(),
                    })
                }
            };

            if c == '\'' {
                // check if escaped
                if let Some('\'') = self.peek_char()? {
                    self.consume_char()?; // Consume escaping '
                    s.push('\'');
                } else {
                    // closing quote
                    break;
                }
            } else {
                s.push(c);
            }
        }

        Ok(NewickToken::QuotedString(s))
    }

    fn scan_comment(&mut self) -> Result<NewickToken, TokenizerError> {
        let mut s = String::new();
        loop {
            match self.consume_char()? {
                Some(']') => break,
                Some('\'') => match self.consume_char()? {
                    Some('[') => s.push('['),
                    Some(']') => s.push(']'),
                    Some(c) => {
                        s.push('\'');
                        s.push(c);
                    }
                    None => {
                        return Err(TokenizerError::UnterminatedComment {
                            span: self.get_current_span(),
                        })
                    }
                },
                Some(c) => s.push(c),
                None => {
                    return Err(TokenizerError::UnterminatedComment {
                        span: self.get_current_span(),
                    })
                }
            }
        }

        Ok(NewickToken::Comment(s))
    }

    /// Parses and consumes one token, bypassing the token-level peek cache.
    fn next_token_uncached(&mut self) -> Result<Option<SpannedToken>, TokenizerError> {
        self.skip_whitespace()?;

        self.current_start = self.byte_offset;
        let c = match self.consume_char()? {
            Some(c) => c,
            None => return Ok(None),
        };

        let token = match c {
            '(' => NewickToken::ParenOpen,
            ')' => NewickToken::ParenClose,
            ',' => NewickToken::Comma,
            ':' => NewickToken::Colon,
            ';' => NewickToken::SemiColon,
            '[' => self.scan_comment()?,
            ']' => {
                return Err(TokenizerError::InvalidCharacter {
                    char: c,
                    span: self.get_current_span(),
                })
            }
            '\'' => self.scan_quoted_string()?,
            _ => self.scan_unquoted_string(c)?,
        };

        let end = self.byte_offset;
        Ok(Some(SpannedToken {
            token,
            span: Span {
                start: self.current_start,
                end,
            },
        }))
    }
}

impl<R: Read> Tokenizer for NewickTokenizer<R> {
    fn next_token(&mut self) -> Result<Option<SpannedToken>, TokenizerError> {
        if let Some(token) = self.peeked_token.take() {
            return Ok(Some(token));
        }
        self.next_token_uncached()
    }

    fn peek(&mut self) -> Result<Option<&SpannedToken>, TokenizerError> {
        if self.peeked_token.is_none() {
            let token = self.next_token()?;
            self.peeked_token = token;
        }
        Ok(self.peeked_token.as_ref())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    pub fn get_token_vec(newick: &str) -> Vec<NewickToken> {
        let mut tokenizer = NewickTokenizer::new(newick.as_bytes());
        let mut tokens = vec![];
        while let Some(next) = tokenizer.next_token().unwrap() {
            tokens.push(next.token)
        }
        tokens
    }

    #[test]
    pub fn test_basic_newick() {
        let newick = "(,);";
        let expected = vec![
            NewickToken::ParenOpen,
            NewickToken::Comma,
            NewickToken::ParenClose,
            NewickToken::SemiColon,
        ];
        let tokens = get_token_vec(newick);
        assert_eq!(tokens, expected)
    }

    #[test]
    pub fn test_standard_newick() {
        let newick = "(A_B_C:123,'this is a name':456);";
        let expected = vec![
            NewickToken::ParenOpen,
            NewickToken::UnquotedString("A_B_C".into()),
            NewickToken::Colon,
            NewickToken::UnquotedString("123".into()),
            NewickToken::Comma,
            NewickToken::QuotedString("this is a name".into()),
            NewickToken::Colon,
            NewickToken::UnquotedString("456".into()),
            NewickToken::ParenClose,
            NewickToken::SemiColon,
        ];
        let tokens = get_token_vec(newick);
        assert_eq!(tokens, expected)
    }

    #[test]
    pub fn test_number_like_branch_lengths() {
        let newick = "(A:-1.23,B:1e+34);";
        let expected = vec![
            NewickToken::ParenOpen,
            NewickToken::UnquotedString("A".into()),
            NewickToken::Colon,
            NewickToken::UnquotedString("-1.23".into()),
            NewickToken::Comma,
            NewickToken::UnquotedString("B".into()),
            NewickToken::Colon,
            NewickToken::UnquotedString("1e+34".into()),
            NewickToken::ParenClose,
            NewickToken::SemiColon,
        ];
        let tokens = get_token_vec(newick);
        assert_eq!(tokens, expected)
    }

    #[test]
    pub fn test_comments() {
        let newick = "(A[comment&^GDSA]:['[']],B);";
        let expected = vec![
            NewickToken::ParenOpen,
            NewickToken::UnquotedString("A".into()),
            NewickToken::Comment("comment&^GDSA".into()),
            NewickToken::Colon,
            NewickToken::Comment("[]".into()),
            NewickToken::Comma,
            NewickToken::UnquotedString("B".into()),
            NewickToken::ParenClose,
            NewickToken::SemiColon,
        ];
        let tokens = get_token_vec(newick);
        assert_eq!(tokens, expected)
    }

    #[test]
    fn skips_whitespace_between_tokens() {
        let tokens = get_token_vec(" \n ( A ,\tB ) ; \r\n");
        assert_eq!(
            tokens,
            vec![
                NewickToken::ParenOpen,
                NewickToken::UnquotedString("A".into()),
                NewickToken::Comma,
                NewickToken::UnquotedString("B".into()),
                NewickToken::ParenClose,
                NewickToken::SemiColon,
            ]
        );
    }

    #[test]
    fn decodes_doubled_quotes_in_quoted_strings() {
        assert_eq!(
            get_token_vec("'Darwin''s finch';"),
            vec![
                NewickToken::QuotedString("Darwin's finch".into()),
                NewickToken::SemiColon,
            ]
        );
    }

    #[test]
    fn emits_number_like_atoms_as_unquoted_strings() {
        assert_eq!(
            get_token_vec("0.5 -2 3E-4 6e+7 -.5 +2 .5 1e+"),
            vec![
                NewickToken::UnquotedString("0.5".into()),
                NewickToken::UnquotedString("-2".into()),
                NewickToken::UnquotedString("3E-4".into()),
                NewickToken::UnquotedString("6e+7".into()),
                NewickToken::UnquotedString("-.5".into()),
                NewickToken::UnquotedString("+2".into()),
                NewickToken::UnquotedString(".5".into()),
                NewickToken::UnquotedString("1e+".into()),
            ]
        );
    }

    #[test]
    fn keeps_numeric_prefixes_attached_to_unquoted_strings() {
        assert_eq!(
            get_token_vec("123abc 1e+invalid"),
            vec![
                NewickToken::UnquotedString("123abc".into()),
                NewickToken::UnquotedString("1e+invalid".into()),
            ]
        );
    }

    #[test]
    fn peek_does_not_consume_the_next_token() {
        let mut tokenizer = NewickTokenizer::new("(A)".as_bytes());

        let first_peek = tokenizer.peek().unwrap().cloned();
        let second_peek = tokenizer.peek().unwrap().cloned();

        assert_eq!(first_peek, second_peek);
        assert_eq!(tokenizer.next_token().unwrap(), first_peek);
        assert_eq!(
            tokenizer.next_token().unwrap().unwrap().token,
            NewickToken::UnquotedString("A".into())
        );
    }

    #[test]
    fn spans_are_half_open_byte_ranges_and_exclude_whitespace() {
        let mut tokenizer = NewickTokenizer::new(" \té:'a''b';".as_bytes());
        let mut tokens = Vec::new();
        while let Some(token) = tokenizer.next_token().unwrap() {
            tokens.push(token);
        }

        assert_eq!(
            tokens,
            vec![
                SpannedToken {
                    token: NewickToken::UnquotedString("é".into()),
                    span: Span { start: 2, end: 4 },
                },
                SpannedToken {
                    token: NewickToken::Colon,
                    span: Span { start: 4, end: 5 },
                },
                SpannedToken {
                    token: NewickToken::QuotedString("a'b".into()),
                    span: Span { start: 5, end: 11 },
                },
                SpannedToken {
                    token: NewickToken::SemiColon,
                    span: Span { start: 11, end: 12 },
                },
            ]
        );
    }

    #[test]
    fn reports_unterminated_constructs_with_their_spans() {
        let mut quoted = NewickTokenizer::new("'abc".as_bytes());
        assert!(matches!(
            quoted.next_token(),
            Err(TokenizerError::UnterminatedQuoteString {
                span: Span { start: 0, end: 4 }
            })
        ));

        let mut comment = NewickTokenizer::new("[abc".as_bytes());
        assert!(matches!(
            comment.next_token(),
            Err(TokenizerError::UnterminatedComment {
                span: Span { start: 0, end: 4 }
            })
        ));
    }

    #[test]
    fn rejects_an_unmatched_closing_bracket() {
        let mut tokenizer = NewickTokenizer::new("]".as_bytes());
        assert!(matches!(
            tokenizer.next_token(),
            Err(TokenizerError::InvalidCharacter {
                char: ']',
                span: Span { start: 0, end: 1 }
            })
        ));
    }
}
