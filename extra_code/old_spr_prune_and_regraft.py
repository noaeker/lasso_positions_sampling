from ete3 import *
import logging



def add_internal_names(original_tree):
    for i, node in enumerate(original_tree.traverse()):
        if not node.is_leaf():
            node.name = "N{}".format(i)
    return original_tree


def prune_at_internal_node(tree_root_pointer, pruned_internal_node_name):
    main_tree_root_pointer_cp = tree_root_pointer.copy()  # working with a copy of the original tree
    pruned_node_pointer = main_tree_root_pointer_cp.search_nodes(name=pruned_internal_node_name)[0]
    pruned_node_parent_name = pruned_node_pointer.up.name
    pruned_subtree_pointer = pruned_node_pointer.detach()
    main_tree_root_pointer_cp.search_nodes(name=pruned_node_parent_name)[0].delete(
        preserve_branch_length=True)
    return pruned_node_parent_name, pruned_subtree_pointer, main_tree_root_pointer_cp


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



#t1 = Tree('((0011:0.1,0012:0.3)N1:0.0625,0027:0.1,(0017:0.1,(0018:0.2,(0029:0.1,((0008:0.1,0006:0.1)N12:0.025,(0002:0.05,0031:0.1)N13:0.025)N11:0.0125)N9:0.15)N7:0.0125)N3:0.03125);', format=1)
#print(t1.get_ascii(attributes=['name'],show_internal=True))

#pruned_node_parent_name, pruned_subtree, remaining_tree = prune_at_internal_node(t1,"N9")


#print("pruned subtree:")
#print(pruned_subtree.get_ascii(attributes=['name'],show_internal=True))
#print("Remaining tree:")
#print(remaining_tree .get_ascii(attributes=['name'],show_internal=True))



#regrafted_tree = regraft_as_sister_of_given_internal_node("N1",
#                                                                              pruned_subtree, remaining_tree)
#print("regrafted tree:")
#print(regrafted_tree.get_ascii(attributes=['name'],show_internal=True))