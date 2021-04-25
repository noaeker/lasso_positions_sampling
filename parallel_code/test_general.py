from spr_prune_and_regraft import *
from config import *
import numpy as np
import matplotlib.pyplot as plt

trees_path = "mlTrees"
print(os.path.exists(""))
trees_object = generate_multiple_tree_object_from_newick(trees_path)
branch_lengths_optimized_per_node = {}
for tree in trees_object:
        for i, node in enumerate(tree.iter_descendants()):
            if node.name not in branch_lengths_optimized_per_node:
                branch_lengths_optimized_per_node[node.name]=[]
            # Do some analysis on node
            else:
                branch_lengths_optimized_per_node[node.name].append(node.dist)
