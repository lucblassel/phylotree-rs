import phylotree as pt
from phylodm import PhyloDM
import sys

TREE = "test_trees/bench_trees/5000.test.nwk"

if __name__ == "__main__":
    if sys.argv[1] == "pdm":
        for _ in range(10):
            tree = PhyloDM.load_from_newick_path(TREE)
            dm = tree.dm()
            taxa = tree.taxa()
        print("Done PDM")
    else:
        for i in range(10):
            pdm = pt.DistanceMatrix.from_newick(TREE)
            dm = pdm.distances
            taxa = pdm.taxa
        print("Done Phylotree")

