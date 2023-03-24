use phylotree::*;
fn main() {
    println!("Generating random binary tree");
    let random = generate_tree(10000, true);
    random.print();
}
