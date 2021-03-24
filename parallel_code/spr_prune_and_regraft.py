from ete3 import *
import logging


class Edge:
    def __init__(self, node_a, node_b):
        self.node_a= node_a
        self.node_b = node_b

    def __str__(self):
        return ("[a={a} b={b}]".format(a=self.node_a, b=self.node_b))

    def __eq__(self, other):
        """Overrides the default implementation"""
        if ((self.node_a==other.node_a) and (self.node_b==other.node_b)) or ((self.node_b==other.node_a) and (self.node_a==other.node_b)):
            return True
        else:
            return False




def print_subtree(tree, log_file, text):
    if log_file:
        log_file.write(text + " visualization: " + "\n" + tree.get_ascii(attributes=['name'],
                                                                                       show_internal=True) + "\n")
        log_file.write(text+ " newick " +str(tree.write(format=1))+"\n")
    else:
        logging.info(text + " visualization: " + "\n"+tree.get_ascii(attributes=['name'], show_internal=True))
        logging.info(str(text+ " newick " +tree.write(format=1)))



def get_list_of_edges(starting_tree):
    edges_list=[]
    main_tree_root_pointer_cp = starting_tree.copy()
    for i, prune_node in enumerate(main_tree_root_pointer_cp.iter_descendants("levelorder")):
        if prune_node.up:
            edge = Edge(node_a=prune_node.name,node_b=prune_node.up.name)
            edges_list.append(edge)
    return edges_list

def get_possible_spr_moves(edges_set):
    possible_moves = []
    for edge1 in edges_set:
        for edge2 in edges_set:
            if not ((edge1.node_a==edge2.node_a) or (edge1.node_b==edge2.node_b) or (edge1.node_b==edge2.node_a) or (edge1.node_a==edge2.node_b)):
                possible_moves.append((edge1,edge2))
    return possible_moves



def add_subtree_to_basetree(subtree_root,basetree_root,regraft_edge,length_regraft_edge,length_pruned_edge):
    future_sister_tree_to_pruned_tree = (basetree_root & regraft_edge.node_a).detach()
    new_tree_adding_pruned_and_future_sister = Tree()
    new_tree_adding_pruned_and_future_sister.add_child( subtree_root.copy(),
                                                       dist=length_pruned_edge)
    new_tree_adding_pruned_and_future_sister.add_child( future_sister_tree_to_pruned_tree, dist=length_regraft_edge/2)
    (basetree_root & regraft_edge.node_b).add_child(new_tree_adding_pruned_and_future_sister, dist=length_regraft_edge / 2)
    basetree_root.unroot()
    return basetree_root




def generate_neighbour(base_tree, possible_move):
    base_tree= base_tree.copy() #not working on original tree
    pruned_edge, regraft_edge= possible_move
    length_regraft_edge = (base_tree & regraft_edge.node_a).dist
    length_pruned_edge = (base_tree & pruned_edge.node_a).dist
    if base_tree.get_common_ancestor(regraft_edge.node_a, pruned_edge.node_a).name==pruned_edge.node_a:
        new_base_tree=(base_tree & pruned_edge.node_a).detach()
        new_subtree_to_be_regrafted = base_tree
        if not (new_subtree_to_be_regrafted & pruned_edge.node_b).name==new_subtree_to_be_regrafted.get_tree_root().name:
            new_subtree_to_be_regrafted.set_outgroup(new_subtree_to_be_regrafted & pruned_edge.node_b)
        (new_subtree_to_be_regrafted & pruned_edge.node_b).delete(preserve_branch_length=True)
        output_tree = add_subtree_to_basetree(new_subtree_to_be_regrafted,new_base_tree,regraft_edge,length_regraft_edge,length_pruned_edge)
    else:
        pruned_subtree = (base_tree & pruned_edge.node_a).detach()
        (base_tree & pruned_edge.node_b).delete(preserve_branch_length=True)
        output_tree = add_subtree_to_basetree(pruned_subtree,base_tree,regraft_edge,length_regraft_edge,length_pruned_edge)
    return output_tree







def add_internal_names(original_tree):
    for i, node in enumerate(original_tree.traverse()):
        if not node.is_leaf():
            node.name = "N{}".format(i)
    return original_tree




def generate_tree_object_from_newick(tree_path):
   starting_tree_object = Tree(newick=tree_path, format=1)
   add_internal_names(starting_tree_object)
   starting_tree_object.get_tree_root().name = "ROOT"
   return starting_tree_object


def generate_multiple_tree_object_from_newick(trees_path):
   with open(trees_path) as trees_path:
       newicks = trees_path.read().split("\n")
       newicks = [t for t in newicks if len(t)>0]
       tree_objects = [generate_tree_object_from_newick(newick) for newick in newicks]
       return tree_objects


def get_tree_string(tree_path):
   tree_object =  Tree(newick=tree_path, format=1)
   return (tree_object.write(format=1))


def compute_tree_divergence(tree_path):
    total_dist=0
    tree = generate_tree_object_from_newick(tree_path)
    for node in tree.iter_descendants():
        # Do some analysis on node
        total_dist = total_dist+node.dist
    return total_dist

def assign_brlen_to_tree_object(tree_object, brlen_list):
    for i,node in enumerate(tree_object.iter_descendants()):
        # Do some analysis on node
        node.dist=brlen_list[i]
    return tree_object



#
#
# with open('test_spr/trees.test','w') as trees:
#     edges_list = get_list_of_edges(t1)
#     moves = get_possible_spr_moves(edges_list)
#     for move in moves:
#         print("original tree")
#         print(t1.get_ascii(attributes=['name'], show_internal=True))
#         print(move[0])
#         print(move[1])
#         tree = generate_neighbour(t1, move)
#         print("Final tree")
#         print(tree.get_ascii(attributes=['name']
#                              , show_internal=True))
#         trees.write(tree.write(format=1)+"\n")
#     edges_list = get_list_of_edges(t2)
#     moves = get_possible_spr_moves(edges_list)
#     for move in moves:
#         print("original tree")
#         print(t2.get_ascii(attributes=['name'], show_internal=True))
#         print(move[0])
#         print(move[1])
#         tree = generate_neighbour(t2, move)
#         print("Final tree")
#         print(tree.get_ascii(attributes=['name']
#                              , show_internal=True))
#         trees.write(tree.write(format=1) + "\n")
#     print(len(moves))
#
#


# print("Alternative tree")
#t3= Tree('(0008:0.1,0006:0.2,((0002:0.05,0031:0.1)N4:0.0625,(0029:0.075,(0018:0.1,((0027:0.05,(0011:0.1,0012:0.2)N15:0.1)N12:0.0125,0017:0.05)N11:0.00625)N9:0.225)N5:0.03125)N3:0.05);',format=1)

#'(0024:0.100000,(0023:0.100000,0027:0.100000):0.100000,(0011:0.100000,0012:0.100000):0.100000);'
# add_internal_names(t3)
# t3.get_tree_root().name="ROOT"
# print(t3.get_ascii(attributes=['name'],show_internal=True))
# edges_list = get_list_of_edges(t3)
# moves = get_possible_spr_moves(edges_list)
# move=(Edge(node_a="N3",node_b="ROOT"),Edge(node_a="N4",node_b="N3") )
# tree = generate_neighbour(t3, move)
# print(tree.get_ascii(attributes=['name'], show_internal=True))
#
#
# with open('trees.test','w') as trees:
#     1==1
#     for move in moves:
#         print("original tree")
#         print(t3.get_ascii(attributes=['name'], show_internal=True))
#         print(move[0])
#         print(move[1])
#         tree = generate_neighbour(t3, move)
#         print("Final tree")
#         print(tree.get_ascii(attributes=['name']
#                              , show_internal=True))
#         trees.write(tree.write(format=1)+"\n")
#     print(len(moves))


#print(t3.get_ascii(attributes=['name'],show_internal=True))
#t3.set_outgroup(t3&"N12")
#print(t3.get_ascii(attributes=['name'],show_internal=True))
#t3.unroot()
#print(t3.get_ascii(attributes=['name'],show_internal=True))

#print((t3&"0008").detach())
#print(t3)
#print((t3&"0006").detach())
#print(t3)




#pruned_node_parent_name, pruned_subtree, remaining_tree = prune_at_internal_node(t3,"N9")
#print("pruned subtree:")
#print(pruned_subtree.get_ascii(attributes=['name'],show_internal=True))
#print("Remaining tree:")
#print(remaining_tree .get_ascii(attributes=['name'],show_internal=True))
#regrafted_tree = regraft_as_sister_of_given_internal_node("N1", pruned_subtree, remaining_tree)

#print("regrafted tree:")
#print(regrafted_tree.get_ascii(attributes=['name'],show_internal=True))





#print("Output tree:")
#t2 = Tree('((0011:0.1,0012:0.3)N1:0.03125,(0029:0.1,((0008:0.1,0006:0.1)N10:0.025,(0002:0.05,0031:0.1)N11:0.025)N7:0.0125)N2:0.15,(0027:0.1,(0017:0.1,0018:0.2125)N9:0.03125)N3:0.015625);', format=1)
#print(t2.get_ascii(attributes=['name'],show_internal=True))


#pruned_node_parent_name, pruned_subtree, remaining_tree = prune_at_internal_node(t2,"N9")





















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
