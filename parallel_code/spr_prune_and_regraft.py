from ete3 import *
import logging


def prune_at_internal_node(tree_root_pointer, pruned_internal_node_name):
    main_tree_root_pointer_cp = tree_root_pointer.copy()  # working with a copy of the original tree
    pruned_node_pointer = main_tree_root_pointer_cp.search_nodes(name=pruned_internal_node_name)[0]
    pruned_node_parent_name = pruned_node_pointer.up.name
    pruned_subtree_pointer = pruned_node_pointer.detach()
    main_tree_root_pointer_cp.search_nodes(name=pruned_node_parent_name)[0].delete(
        preserve_branch_length=True)
    return pruned_node_parent_name, pruned_subtree_pointer, main_tree_root_pointer_cp

def print_subtree(tree, log_file, text):
    if log_file:
        log_file.write(text + " visualization: " + "\n" + tree.get_ascii(attributes=['name'],
                                                                                       show_internal=True) + "\n")
        log_file.write(text+ " newick " +str(tree.write(format=1))+"\n")
    else:
        logging.info(text + " visualization: " + "\n"+tree.get_ascii(attributes=['name'], show_internal=True))
        logging.info(str(text+ " newick " +tree.write(format=1)))

def regraft_as_sister_of_given_internal_node(future_sister_internal_node_name,pruned_subtree_pointer, tree_root_pointer_cp):

    pruned_subtree_pointer_cp = pruned_subtree_pointer.copy()
    tree_root_pointer_cp_cp = tree_root_pointer_cp.copy()
    future_sister_internal_node_pointer = tree_root_pointer_cp_cp.search_nodes(name=future_sister_internal_node_name)[0]
    future_sister_internal_node_dist_to_parent = future_sister_internal_node_pointer.dist
    future_sister_internal_node_parent_pointer = future_sister_internal_node_pointer.up
    pruned_sister_internal_node_pointer = future_sister_internal_node_pointer.detach()
    #print("pruned_sister_internal_node_pointer")
    #print_subtree(pruned_sister_internal_node_pointer)
    #print("Remaining tree")
    #print_subtree( tree_root_pointer_cp_cp)
    new_tree_adding_pruned_and_future_sister = Tree()
    new_tree_adding_pruned_and_future_sister.add_child(pruned_sister_internal_node_pointer,
                                                       dist=future_sister_internal_node_dist_to_parent / 2)
    new_tree_adding_pruned_and_future_sister.add_child(
        pruned_subtree_pointer_cp)
    new_tree_adding_pruned_and_future_sister.name = "Temp"
    #print("new_tree_adding_pruned_and_future_sister")
    #print_subtree(new_tree_adding_pruned_and_future_sister)
    #print("future_sister_internal_node_parent_pointer:")
    #print_subtree(future_sister_internal_node_parent_pointer)
    #print("tree_root_pointer_cp_cp:")
    #print_subtree(tree_root_pointer_cp_cp)
    future_sister_internal_node_parent_pointer.add_child(new_tree_adding_pruned_and_future_sister,
                                                         dist=future_sister_internal_node_dist_to_parent / 2)
    fixed_future_sister_internal_node_parent_pointer = fix_structure(tree_root_pointer_cp_cp)
    return fixed_future_sister_internal_node_parent_pointer

def add_internal_names(original_tree):
    for i, node in enumerate(original_tree.traverse()):
        if not node.is_leaf():
            node.name = "N{}".format(i)
    return original_tree


def fix_structure(tree):
    if tree.search_nodes(name="Temp"):
        print("tree before outgroup")
        temp_node_pointer = tree.search_nodes(name="Temp")[0]
        tree.set_outgroup(temp_node_pointer)
        print("tree after outgroup")
        new_tree_with_3_children = Tree()
        new_tree_with_3_children.add_child(temp_node_pointer.children[0])
        new_tree_with_3_children.add_child(temp_node_pointer.children[1])
        for child in tree.children:
            if child.name!="Temp":
                new_tree_with_3_children.add_child(child)
        add_internal_names(new_tree_with_3_children)
        new_tree_with_3_children.get_tree_root().name = "ROOT"
        assert (len(new_tree_with_3_children.children)==3)
        return new_tree_with_3_children
    return tree
    # if tree.search_nodes(name="Temp"):
    #     temp_node_pointer = tree.search_nodes(name="Temp")[0]
    #     tree.set_outgroup(temp_node_pointer)
    #     #print("tree after outgroup:")
    #     #print_subtree(tree)
    #     temp_node_pointer.delete(
    #         preserve_branch_length=True)
    #     #add_internal_names(tree)
    #     return tree
    # return tree

def generate_tree_object(tree_path):
   starting_tree_object = Tree(newick=tree_path, format=1)
   add_internal_names(starting_tree_object)
   starting_tree_object.get_tree_root().name = "ROOT"
   return starting_tree_object


def compute_tree_divergence(tree_path):
    total_dist=0
    tree = generate_tree_object(tree_path)
    for node in tree.iter_descendants():
        # Do some analysis on node
        total_dist = total_dist+node.dist
    return total_dist

def assign_brlen(tree_path,brlen_list):
    tree = generate_tree_object(tree_path)
    for i,node in enumerate(tree.iter_descendants()):
        # Do some analysis on node
        node.dist=brlen_list[i]
    tree.write(format=1, outfile=tree_path)




#t = Tree('(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1,E:0.1);', format=1)

#t.write(format=1, outfile="tree_test")

#assign_brlen("tree_test",[0.2,0.3,0.4,0.5,0.6,0.7,0.8])







#print(E.up==B.up)
#add_internal_names(t)
#print(t.get_ascii(attributes=['name'],  show_internal=True))

#print(compute_tree_divergence(t))


# def get_SPR_neighbours(starting_tree):
#     spr_neighbours = []
#     for i, prune_node in enumerate(starting_tree.iter_descendants("levelorder")):
#         pruned_node_parent_name, pruned_subtree, remaining_tree = prune_at_internal_node(starting_tree,
#                                                                                          prune_node.name)  # subtree1 is the pruned subtree. subtree2 is the remaining subtree
#
#         for j, rgft_node in enumerate(remaining_tree.iter_descendants("levelorder")):
#                 print("starting_tree:")
#                 print_subtree(starting_tree,reg,
#                 print("prune_node:")
#                 print_subtree(prune_node,None,"")
#                 regrafted_tree = regraft_as_sister_of_given_internal_node(rgft_node.name,
#                                                                           pruned_subtree, remaining_tree)
#                 print("regraft node:")
#                 print_subtree(rgft_node,,
#                 print("regrafted tree:")
#                 print_subtree(regrafted_tree,,
#                 spr_neighbours.append(regrafted_tree)



#t = Tree('(((A:1.828378,B:0.750678):0.0150676,C:0.82589):0.232357,D:0.464714,E:0.707315);', format=1)
#A=t&"A"
# B=t&"B"
# E=t&"E"
# add_internal_names(t)
#print(t.get_ascii(attributes=['name'],  show_internal=True))
# pruned_node_parent_name, pruned_subtree_pointer, tree_root_pointer_cp = prune_at_internal_node(t, "N4")
# print("pruned subtree:")
# print(pruned_subtree_pointer.get_ascii(attributes=['name'],  show_internal=True))
# print("Remaining tree:")
# print(tree_root_pointer_cp.get_ascii(attributes=['name'],  show_internal=True))
# print("regrafting at E")
# regrafted_tree = regraft_as_sister_of_given_internal_node("E",pruned_subtree_pointer, tree_root_pointer_cp)
#
# print("Rrgrafted tree:")
# print(regrafted_tree.get_ascii(attributes=['name'],  show_internal=True))
# print(regrafted_tree.write(format=1))
# fix_structure(regrafted_tree)
# print(fix_structure(regrafted_tree).write(format=1))
