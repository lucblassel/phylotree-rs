use phylotree::*;

fn main() {
    let random = generate_tree(100, true);
    let newick = random.to_newick();
    println!("{newick}")
}
