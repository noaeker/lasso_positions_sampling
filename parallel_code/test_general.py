from spr_prune_and_regraft import *
from config import *
import numpy as np
import matplotlib.pyplot as plt

trees_path = "mlTrees"
trees_object = generate_multiple_tree_object_from_newick(trees_path)
branch_lengths_optimized_per_node = {}
for tree in trees_object:
        for i, node in enumerate(tree.iter_descendants()):
            if node.name not in branch_lengths_optimized_per_node:
                branch_lengths_optimized_per_node[node.name]=[]
            # Do some analysis on node
            else:
                branch_lengths_optimized_per_node[node.name].append(node.dist)

print(branch_lengths_optimized)
print(np.var(branch_lengths_optimized))
#plt.hist(branch_lengths_optimized)
#plt.show()
branch_lengths_exponential = sample_exp(len(branch_lengths_optimized),1)
print(branch_lengths_exponential)
print(np.var(branch_lengths_exponential))
#plt.hist(branch_lengths_exponential)
#plt.show()


fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
axs[0].hist(branch_lengths_optimized)
axs[1].hist(branch_lengths_exponential)
plt.show()