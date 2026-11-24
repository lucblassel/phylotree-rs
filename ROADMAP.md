# Roadmap to the First Rust Beta Release

This roadmap focuses on tightening the Rust API, fixing confirmed correctness issues, enforcing tree and matrix invariants, and preparing a reliable first beta release.

Python bindings, Python packaging, performance-only optimizations, and new domain functionality are intentionally out of scope.

## Release target

The unreleased branch contains breaking changes relative to `0.1.3`. Publish it as **`0.2.0`** rather than another `0.1.x` release.

A beta is ready when:

- [ ] All release-blocking correctness issues below are fixed.
- [ ] Public safe APIs do not panic on malformed input or invalid operations.
- [x] Every successful tree mutation preserves documented invariants.
- [x] Every failed tree mutation leaves the tree unchanged.
- [x] Public constructors cannot create invalid `DistanceMatrix` values.
- [ ] Unit tests, doctests, Clippy, formatting, and documentation checks pass in CI.
- [ ] The README examples compile and match the released API.
- [ ] Breaking changes from `0.1.3` are documented in the changelog.

---

## 1. Fix release-blocking correctness bugs

### 1.1 Fix `DistanceMatrix::to_map`

`DistanceMatrix::to_map` currently iterates over diagonal pairs and unwraps an error from `get_pair_index`, so it panics for every non-empty matrix.

Relevant code: `src/distance.rs`.

- [x] Decide whether `to_map` returns all ordered pairs, including diagonal zeroes, or only unique off-diagonal pairs.
- [x] Document that representation explicitly.
- [x] Remove the internal `unwrap` and implement the selected behavior without panicking.
- [x] Add tests for empty, one-taxon, and multi-taxon matrices.
- [x] Test diagonal, symmetric, and off-diagonal entries explicitly.

**Acceptance criteria**

- `to_map` succeeds for every valid matrix size, including zero and one.
- The returned keys and diagonal behavior match the documentation.
- No safe caller input can trigger a panic.

### 1.2 Make strict PHYLIP parsing panic-free and genuinely strict

`from_phylip_strict` can index outside a row when the input contains more rows than declared. Row parsing also truncates extra fields with `.take(...)`, allowing malformed input to be silently accepted.

Relevant code: `src/distance.rs`.

- [x] Reject extra rows before indexing by row number.
- [x] Parse complete rows rather than silently truncating extra columns.
- [x] Require the exact expected number of values per row.
- [x] Reject too few and too many rows.
- [x] Reject duplicate taxon labels.
- [x] Validate diagonal values in square matrices without unchecked indexing.
- [x] Reject non-finite values in strict PHYLIP input.
- [x] Audit parser code for indexing and `unwrap` calls reachable from malformed input.

**Regression tests**

- [x] Empty input.
- [x] Declared size zero.
- [x] Missing rows.
- [x] Extra rows.
- [x] Missing columns.
- [x] Extra columns.
- [x] Malformed trailing columns.
- [x] Non-zero diagonal.
- [x] Asymmetric square matrix.
- [x] Duplicate taxon labels.
- [x] Non-finite distances.

**Acceptance criteria**

- Every malformed input returns a typed parse error rather than panicking.
- The strict parser accepts only inputs with exactly the documented shape.

### 1.3 Validate random-tree generator sizes

The generators currently handle zero and one leaf inconsistently. `generate_tree(0, ...)` underflows, and `generate_yule(0, ...)` can continue growing indefinitely.

Relevant code: `src/lib.rs` and the CLI argument definitions.

- [x] Define the minimum supported leaf count.
- [x] Return a dedicated typed error for unsupported sizes.
- [x] Apply the same rules to every generator.
- [x] Make one-leaf behavior consistent, including tip naming.
- [x] Make tip numbering consistent across generators.
- [x] Enforce the same range in CLI argument validation.
- [x] Remove arithmetic that can underflow for invalid sizes.

**Regression tests**

- [x] Test sizes 0, 1, and 2 for every generator.
- [x] Confirm valid generated trees have the requested number of live leaves.
- [x] Confirm all generated leaves have unique names.

### 1.4 Remove public unconditional panics

`DistanceMatrix::neighbor_joining` is public but implemented with `todo!()`.

Relevant code: `src/distance.rs`.

- [x] Remove the method from the public API for the beta, hide it behind an explicitly experimental feature, or return a typed unsupported-operation error.
- [x] Do not leave `todo!`, `unimplemented!`, or unconditional `panic!` in public safe library methods.
- [x] Audit public library methods for input-reachable `unwrap`, `expect`, `panic`, `todo`, `unimplemented`, and `unreachable` calls.

Implementing neighbor joining is not required for this beta roadmap; removing or explicitly hiding the incomplete API is sufficient.

**Acceptance criteria**

- Every public safe library method either succeeds or returns a documented error for invalid caller input.

---

## 2. Restore and enforce tree invariants

### 2.1 Define the invariants of a valid `Tree`

Document the internal guarantees that all constructors and mutation methods must preserve.

1. [x] Exactly one live root exists for every non-empty connected tree.
2. [x] Every live non-root node has one live parent.
    - [x] Made the structural set_parent, add_child, remove_child node methods crate-private.
    - [x] Added a Tree::remove_child method that uses the Node methods. Right now it deletes the node. 
    - [x] ~We need to decide if the remove_child method can work on internal nodes and therefore prune the subtree or if we restrict it to leaves.~ `remove_child` only works on leaves, and the error messages tells users to use `prune` for subtrees.
3. [x] Every child ID refers to a live node in the same tree.
4. [x] Parent and child relationships are reciprocal.
5. [x] A parent contains each child at most once.
6. [x] Branch lengths agree between the child and parent representations.
7. [x] Every live node is reachable from the root.
8. [x] The topology is acyclic.
9. [x] Every stored node depth equals its actual distance in edges from the root.
10. [x] Deleted storage slots are ignored by semantic tree operations.
11. [x] Cached or derived indexes are either valid or explicitly invalidated.
12. [x] Add an internal invariant-checking helper for tests.
13. [x] Consider exposing a public `Tree::validate()` method if users can construct or deserialize trees through low-level APIs.

### 2.2 Make tombstone handling consistent

Deleted nodes are reset to empty nodes and retained in the backing vector. Some operations filter them while others count them as leaves or roots.

Relevant code: `src/tree/node.rs` and `src/tree/tree_impl.rs`.

- [x] Preserve the storage-slot ID when marking a node deleted, or otherwise ensure deleted nodes cannot masquerade as node zero internally.
- [x] Make `n_leaves` count only live leaves.
- [x] Make root lookup ignore deleted nodes.
- [x] Make node searches ignore deleted nodes.
- [x] Make tree length and whole-tree metrics ignore deleted nodes.
- [x] Audit every direct iteration over `self.nodes` and decide whether deleted nodes must be filtered.
- [x] Ensure taxon-index construction and unique-name validation use only live leaves.
- [x] Ensure serialization and NEXUS counts use live node and leaf counts.
- [x] Define whether `Tree::size` means live nodes or allocated storage slots.
- [x] Prefer a user-facing `len` or `node_count` for live nodes; keep storage length private unless it has a clear use.

**Regression tests**

- [x] Counts and names after pruning a leaf.
- [x] Counts and names after pruning an internal subtree.
- [x] Counts and names after compression.
- [x] Root lookup after attempting to prune the root.
- [x] Taxon-index rebuilding after prune and compression.
- [x] Serialization and NEXUS taxon counts after deletion.

**Acceptance criteria**

- All public queries agree on which nodes are live.
- Pruning and compression do not cause false duplicate-name or unnamed-leaf errors.

### 2.3 Keep depths correct after every topology mutation

`compress_node`, `resolve`, and `merge_children` currently move nodes without updating stored depths. This can produce incorrect Sackin values and incorrect bipartition depth metadata.

Relevant code: `src/tree/tree_impl.rs`.

- [x] Update depths after `compress`.
- [x] Update depths after `resolve`.
- [x] Update depths after `merge_children`.
- [x] Audit `prune` and future topology-changing helpers for the same requirement.
- [x] Decide whether to update only affected subtrees or call `reset_depths` once after a successful compound mutation.
- [x] Keep depth repair inside mutation methods rather than requiring callers to remember it.

**Regression tests**

- [x] For every live parent-child pair, assert `child.depth == parent.depth + 1` after each mutation.
- [x] Confirm Sackin values after compressing unary nodes.
- [x] Confirm depths after resolving a polytomy.
- [x] Confirm depths after merging siblings.
- [x] Confirm bipartition depth metadata after topology changes.

### 2.4 Make mutations atomic on error

`merge_children` can partially remove a child and then return an error, leaving the tree inconsistent.

Relevant code: `src/tree/tree_impl.rs`.

- [x] Reject attempts to merge a node with itself before mutation.
- [x] Validate both node IDs before mutation.
- [x] Validate both parent relationships before mutation.
- [x] Validate that the parent contains both children exactly once.
- [x] Validate all other fallible preconditions before making structural changes.
- [x] Perform the commit only after validation succeeds.
- [x] Audit other compound mutation methods for partial updates before possible errors.

**Regression tests**

- [x] Merging a node with itself returns an error and leaves the tree unchanged.
- [x] Merging non-siblings returns an error and leaves the tree unchanged.
- [x] Merging missing or deleted nodes returns an error and leaves the tree unchanged.
- [x] Failed compression leaves the tree unchanged.
- [x] Compare serialized form and invariant state before and after each failed operation.

**Acceptance criteria**

- Every mutation returning `Err` has no externally observable effect.

### 2.5 Clarify and correct binary-tree validation

`is_binary` currently rejects only nodes with more than two children, so unary internal nodes can be classified as binary.

Relevant code: `src/tree/tree_impl.rs`.

- [x] Define rooted binary-tree degree requirements.
- [x] Define unrooted binary-tree degree requirements, including the virtual root convention.
- [x] Decide how zero-node and one-node trees are classified.
- [x] Make `is_binary` enforce the documented definition.
- [x] Ensure Colless, Sackin, and cherry calculations use the same definition.
- [x] Remove metric logic that silently treats a missing second child as an empty subtree unless that is explicitly part of the definition.

**Regression tests**

- [x] Empty tree.
- [x] Single-node tree.
- [x] Unary root.
- [x] Unary internal node.
- [x] Proper rooted binary tree.
- [x] Rooted polytomy.
- [x] Proper unrooted binary tree under the crate's virtual-root convention.
- [x] Unrooted multifurcation.

### 2.6 Make branch-length semantics explicit

`height` and `diameter` can fall back from branch-length sums to edge counts on individual paths. In a partially weighted tree this can compare values expressed in different units.

Relevant code: `src/tree/tree_impl.rs`.

- [x] Split topological and branch-length variants, or require an explicit metric argument.
- [x] Make branch-length variants return `MissingBranchLengths` when required path lengths are incomplete.
- [x] Keep edge-count variants independent of branch-length presence.
- [x] Define whether negative branch lengths are valid.
- [x] Reject non-finite branch lengths and support values during parsing unless there is a documented reason to preserve them.
- [x] Apply the same semantics to distance (✅), height (✅), diameter (✅), and comparison APIs.

A possible API shape is:

```rust
pub enum TreeDistanceMetric {
    EdgeCount,
    BranchLength,
}
```

or separate methods such as `height_by_edges` and `height_by_length`.

**Acceptance criteria**

- No result combines or compares topological counts with branch-length units implicitly.
- Partially weighted trees have deterministic, documented behavior.

### 2.7 Fix Newick comment round trips

Comment tokenization and serialization currently use inconsistent apostrophe and closing-bracket escaping rules.

Relevant code: `src/tree/newick/tokenizer.rs` and `src/tree/newick/serializer.rs`.

- [x] Define one documented escaping convention for comments.
- [x] Make tokenization and serialization exact inverses for supported comments.
- [x] Verify that trailing apostrophes do not consume the closing bracket.
- [x] Add round-trip tests for apostrophes, brackets, repeated comments, Unicode, and empty comments.
- [x] Add property-based or table-driven parser/serializer round-trip tests.

---

## 3. Tighten the public Rust API

These changes are intentionally grouped before `0.2.0` because they are easier to make before downstream users depend on the current representation.

### 3.1 Protect `Tree` structure from accidental corruption

`Node` exposes structural fields publicly, and `Tree::get_mut` allows callers to modify them without maintaining reciprocal links, depths, or indexes.

Relevant code: `src/tree/node.rs` and `src/tree/tree_impl.rs`.

- [x] Make `Node::id`, `parent`, `children`, and structural edge bookkeeping private.
- [x] Provide read-only accessors for structural data.
- [x] Route topology changes through invariant-preserving `Tree` methods.
- [x] Decide whether metadata fields such as name, support, and comment remain public or receive setters.
- [x] Prevent `Tree::add` from importing stale IDs, parents, children, deletion flags, or caches from a detached `Node`.
- [x] Clearly distinguish low-level disconnected construction from normal rooted-tree construction.
- [x] Ensure setting a branch length updates both ends of the edge.
- [x] Ensure clearing a branch length is supported and updates both ends; `set_child_edge(None)` must not silently leave an old value in place.

**Acceptance criteria**

- Safe external code cannot directly create cycles or mismatched parent-child links by editing a `Node`.

### 3.2 Introduce a `NodeId` newtype

`NodeId` is currently a `usize` alias, allowing accidental mixing with arbitrary indexes.

Relevant code: `src/tree/mod.rs` and all node-ID call sites.

- [x] Replace the alias with a small `Copy` newtype.
- [x] Keep the inner storage index private where practical.
- [x] Implement ordering, hashing, formatting, and lossless `usize` conversion traits, plus an explicit storage-index accessor.
- [x] Accept `NodeId` by value rather than `&NodeId` in public methods.
- [x] Document that IDs are tree-local.
- [x] Document whether IDs remain stable after prune, compression, and other mutations.

### 3.3 Make `DistanceMatrix` states valid by construction

Relevant code: `src/distance.rs`.

- [x] Remove the redundant stored `size` field and make `taxa` private.
- [x] Expose `len`, `is_empty`, and `taxa` accessors.
- [x] Replace unchecked public construction with a checked constructor returning `Result`.
- [x] Validate compact matrix length against `n * (n - 1) / 2`.
- [x] Handle sizes zero and one without underflow.
- [x] Use checked arithmetic for compact storage size.
- [x] Reject duplicate taxon labels.
- [x] Validate taxon existence before returning a same-taxon distance of zero.
- [x] Add index-based accessors so algorithms do not repeatedly scan names.
- [x] Keep compact construction checked internally; no unchecked constructor is exposed.

**Acceptance criteria**

- Every publicly constructible `DistanceMatrix` satisfies its shape and taxon invariants.
- Public accessors cannot panic because of a malformed matrix created through safe public APIs.

### 3.4 Remove unsound `Eq` semantics from `Node`

`Node::PartialEq` uses epsilon-based floating-point comparison and then implements `Eq`. Approximate equality is not an equivalence relation, and NaN is not reflexive.

Relevant code: `src/tree/node.rs`.

- [x] Remove `Eq`, or change `PartialEq` to exact equivalence semantics.
- [x] If approximate comparison is useful, expose an explicitly named method with a documented tolerance.
- [x] Define whether node equality is metadata equality, detached-node equality, or structural identity.
- [x] Avoid implying whole-tree structural equality from a weak node comparison.

### 3.5 Stop leaking private partition representations

Legacy public methods expose private aliases backed by `FixedBitSet` and collection implementation details.

Relevant code: `src/tree/bipartition.rs` and `src/tree/tree_impl.rs`.

- [x] Deprecate or remove public methods returning raw `Partition` or `PartitionSet` aliases.
- [x] Prefer the opaque public `Bipartition` type.
- [x] Ensure every public return and argument type is intentionally public and documented.
- [x] Add any small accessors required to make `Bipartition` sufficient without exposing its bitset.

### 3.6 Normalize naming and result types

Perform a focused consistency pass rather than a broad redesign.

- [ ] Choose consistent count names: `len`, `node_count`, `leaf_count`, and `cherry_count`, or clearly document the existing convention.
- [ ] Consider dropping `get_` from computed/query names such as root, leaves, descendants, distance, and common ancestor.
- [ ] Rename or clarify `compare_topologies`, since it also computes branch-length metrics.
- [ ] Make direct RF methods and `Comparison::rf` use consistent numeric types where possible.
- [ ] Replace ambiguous tuples with named structs where field meaning is not obvious, especially tree distance results.
- [ ] Use `impl AsRef<Path>` for public path-based APIs.
- [ ] Add symmetric reader/writer APIs where existing string and file APIs are inconsistent.

Compatibility aliases with deprecation notices may be retained where they make migration from `0.1.3` straightforward.

### 3.7 Make public enums evolvable

- [ ] Review public error and configuration enums for `#[non_exhaustive]`.
- [ ] Add `#[non_exhaustive]` where downstream exhaustive matching should not be part of the stability contract.
- [ ] Ensure errors preserve useful sources rather than collapsing failures into general string variants.
- [ ] Replace `GeneralError(&'static str)` with specific variants for known failure conditions.

### 3.8 Decouple core API types from CLI dependencies

This is API and dependency tightening rather than a new feature.

Relevant code: `Cargo.toml`, `src/lib.rs`, and `src/distr.rs`.

- [ ] Remove unconditional `clap::ValueEnum` coupling from core enums, or gate it behind the `cli` feature.
- [ ] Make `clap` optional if it is only required by the binary.
- [ ] Review presentation-only dependencies such as `ptree` for feature gating.
- [ ] Confirm `cargo check --no-default-features` builds a lean core library.
- [ ] Confirm the default feature set matches the intended experience for library users.

---

## 4. Strengthen regression and invariant testing

### 4.1 Add mutation invariant tests

For `prune`, `compress`, `resolve`, and `merge_children`, test after success that:

- [x] There is exactly one expected live root.
- [x] Every live node is reachable.
- [x] Every parent-child relationship is reciprocal.
- [x] No child occurs twice under one parent.
- [x] Edge lengths agree in both representations.
- [x] Stored depths are correct.
- [x] Live node and leaf counts are correct.
- [x] Taxon-index rebuilding succeeds when leaf names are valid.
- [x] Newick serialization can be parsed back into an equivalent valid tree.

For every documented error path:

- [x] Snapshot or clone the tree before the operation.
- [x] Assert that an error is returned.
- [x] Assert that the tree is unchanged.
- [x] Re-run the invariant checker.

### 4.2 Add panic-resistance tests

- [ ] Exercise every public constructor with boundary sizes.
- [ ] Exercise matrix parsers with malformed dimensions and rows.
- [ ] Exercise public methods with missing, deleted, and out-of-range node IDs.
- [ ] Exercise empty and single-node trees.
- [ ] Exercise partially weighted trees.
- [ ] Exercise duplicate and missing leaf names.
- [ ] Ensure safe public APIs return errors rather than unwinding.

Consider using `proptest` or a similar tool for compact matrix dimensions, Newick round trips, and mutation sequences, but do not make property testing a release blocker if focused regression tests cover the confirmed failures.

### 4.3 Improve test effectiveness

- [ ] Replace tests that compare a method with the same method call rather than an independent expected result.
- [ ] Add small hand-calculated reference trees for every metric.
- [ ] Verify rooted and unrooted comparison expectations independently.
- [ ] Test public behavior rather than only internal intermediate values.
- [ ] Keep regression tests next to the bug they prevent and name them after the failure mode.

### 4.4 Add deep-tree robustness tests

Even if iterative traversal refactoring is deferred, establish expected behavior for highly unbalanced trees.

- [ ] Add a reasonably deep caterpillar-tree traversal test.
- [ ] Test parse/serialize round trips on deep trees.
- [ ] Test pruning and depth rebuilding on deep trees.
- [ ] Document any practical depth limit if stack-recursive implementations remain in the beta.

---

## 5. Documentation, CI, and beta release gates

### 5.1 Repair public documentation

Relevant code: `README.md` and Rustdoc throughout `src`.

- [ ] Import every type used in README examples.
- [ ] Remove undefined variables from examples.
- [ ] Update comparison examples to the current rooted/unrooted API.
- [ ] Update the advertised crate version.
- [ ] Correct stale method names and signatures.
- [ ] Document node-ID stability and tree-locality.
- [ ] Document binary-tree definitions.
- [ ] Document behavior for empty and single-node trees.
- [ ] Document branch-length requirements for every metric.
- [ ] Document deletion/live-node count semantics.
- [ ] Fix the broken `from_phylip_strict` Rustdoc link in `src/distance.rs`.
- [ ] Review spelling and terminology in public errors and docs.

### 5.2 Modernize Rust CI sufficiently for the beta

Relevant code: `.github/workflows/rustcheck.yaml`.

- [ ] Replace archived `actions-rs` actions with maintained Rust toolchain and cache actions.
- [ ] Update `actions/checkout` to a maintained version.
- [ ] Run `cargo fmt --all -- --check`.
- [ ] Run `cargo check --all-targets`.
- [ ] Run `cargo test --all-targets`.
- [ ] Run `cargo clippy --all-targets -- -D warnings`.
- [ ] Run `cargo doc --no-deps` with broken links denied.
- [ ] Run `cargo check --no-default-features`.
- [ ] Run the intended supported feature combinations.
- [ ] Add an MSRV to `Cargo.toml` and test it if the project intends to promise one.

A useful documentation gate is:

```shell
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps
```

### 5.3 Prepare the breaking release

- [ ] Set the package version to `0.2.0`.
- [ ] Add a `0.2.0` changelog section.
- [ ] List breaking API changes from `0.1.3`.
- [ ] Provide migration examples for renamed methods and changed comparison signatures.
- [ ] Confirm Cargo package metadata points to the correct repository and homepage.
- [ ] Clarify the SPDX license expression as `GPL-3.0-only` or `GPL-3.0-or-later` according to project intent.
- [ ] Run `cargo package` and inspect the packaged file list.
- [ ] Run tests against the packaged crate where practical.
- [ ] Confirm no generated artifacts, benchmark outputs, or local scratch files are included.

---

## Recommended execution order

1. [x] Fix `DistanceMatrix::to_map` and matrix construction invariants.
2. [x] Make PHYLIP parsing strict and panic-free.
3. [x] Define tree invariants and add a reusable test checker.
4. [x] Fix tombstone filtering and live node counts.
5. [x] Repair depths after topology mutations.
6. [x] Make compound mutations atomic on error.
7. [x] Correct binary-tree validation.
8. [x] Make branch-length versus edge-count semantics explicit.
9. [x] Validate random-tree generator sizes.
10. [x] Fix Newick comment round trips and non-finite numeric handling.
11. [x] Remove incomplete or panic-only public APIs.
12. [x] Make `Tree` structure and `DistanceMatrix` fields invariant-preserving.
13. [x] Introduce the `NodeId` newtype and normalize node-ID signatures.
14. [ ] Remove unsound equality and private representation leaks.
15. [ ] Perform the focused naming and error-enum cleanup.
16. [ ] Add regression, panic-resistance, and deep-tree tests.
17. [ ] Repair README and Rustdoc.
18. [ ] Strengthen CI and complete the `0.2.0` release checklist.

## Explicitly deferred

The following are intentionally not part of this beta roadmap:

- Python API or packaging work
- New tree-building algorithms
- Consensus-tree functionality
- New file formats
- Rerooting and other new tree operations
- Performance-only rewrites that are not required for correctness

These can be planned after the core Rust API and invariants are stable.
