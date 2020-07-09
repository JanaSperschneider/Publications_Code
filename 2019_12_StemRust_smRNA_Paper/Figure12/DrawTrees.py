from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace

def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=20)
        faces.add_face_to_node(N, node, 0, position="aligned")

# Set dashed blue lines in all leaves
nst1 = NodeStyle()
nst1["bgcolor"] = "LightSteelBlue"
nst2 = NodeStyle()
nst2["bgcolor"] = "Moccasin"
nst3 = NodeStyle()
nst3["bgcolor"] = "DarkSeaGreen"
nst4 = NodeStyle()
nst4["bgcolor"] = "Khaki"
nst5 = NodeStyle()
nst5["bgcolor"] = "Lavender"
nst6 = NodeStyle()
nst6["bgcolor"] = "Gold"
nst7 = NodeStyle()
nst7["bgcolor"] = "PowderBlue"
nst8 = NodeStyle()
nst8["bgcolor"] = "Linen"

f = open('AGOs_UniProt_ReferenceSet_Muscle_3.8.31.fasttree', 'r')
content = f.readlines()
f.close()
tree_string = "".join(content)
t = Tree(tree_string)
#t = Tree( "(((((Ml_AGO2A:0.00187,Ml_AGO2B:0.00055)1.000:0.19006,(Mlp_AGO1:0.04247,(MlAG3A:0.00740,MlAG3B:0.01238)1.000:0.07012)1.000:0.18319)1.000:0.53996,(((Mlp_AGO2:0.02370,(Ml_AGO1A:0.0,Ml_AGO1B:0.0):0.03382)1.000:0.28776,(PGT21_AGO2A:0.00055,(PGT_AGO1:0.01647,PGT21_AGO2B:0.00054)0.971:0.00640)1.000:0.40103)1.000:0.73105,((At_AGO6:0.41636,(At_AGO4:0.19578,(At_AGO8:0.26817,At_AGO9:0.16101)0.973:0.09249)0.998:0.21283)1.000:0.57476,((At_AGO7:0.72975,(At_AGO2:0.16484,At_AGO3:0.19311)1.000:0.82454)0.996:0.22989,(At_AGO5:0.37673,(At_AGO1:0.16478,At_AGO10:0.26880)1.000:0.19807)1.000:0.26777)0.931:0.12518)1.000:0.58386)1.000:0.66830)1.000:0.39945,(PGT21_AGO1A:0.00053,(PGT21_AGO1B:0.00358,PGT_AGO2:0.00054)0.797:0.00200)0.989:0.08187)0.988:0.09125,(Pst_AGO2A:0.08160,Pst_AGO2B:0.00113)0.923:0.00905,(Pst_AGO1B:0.05708,Pst_AGO1A:0.00055)0.956:0.01219);" )

# Calculate the midpoint node
R = t.get_midpoint_outgroup()
# and set it as tree outgroup
t.set_outgroup(R)
print t

for n in t.traverse():
    n.dist = 0

D_leaf_color = {"PGT21_AGO2A":"red", "PGT21_AGO2B":"red", "PGT21_AGO1A":"red", "PGT21_AGO1B":"red"}

for node in t.traverse():
    # Hide node circles
    #node.img_style['size'] = 0
    
    if node.is_leaf():
        color = D_leaf_color.get(node.name, None)
        if color:
            name_face = TextFace(node.name, fgcolor=color, fsize=18)
            node.add_face(name_face, column=0, position='branch-right')
        else:
            name_face = TextFace(node.name, fgcolor="black", fsize=18)
            node.add_face(name_face, column=0, position='branch-right')

n1 = t.get_common_ancestor("PGT21_AGO2A", "Ml_AGO2A")
n1.set_style(nst1)

n2 = t.get_common_ancestor("PGT21_AGO1A", "Ml_AGO1A")
n2.set_style(nst2)

n3 = t.get_common_ancestor("At_AGO1", "At_AGO5")
n3.set_style(nst3)

n4 = t.get_common_ancestor("At_AGO8", "At_AGO6")
n4.set_style(nst4)

n5 = t.get_common_ancestor("At_AGO2", "At_AGO7")
n5.set_style(nst5)

#----------------------------------------------------------------------------------------------------------------
ts = TreeStyle()
ts.mode = "c" # draw tree in circular mode

ts.scale =  200 # 120 pixels per branch length unit
ts.show_branch_support = True
#ts.branch_vertical_margin = 0 # 10 pixels between adjacent branches
#ts.layout_fn = layout

# Disable the default tip names config
ts.show_leaf_name = False
ts.show_scale = False

t.render("AGOs.png", dpi=600, tree_style=ts)
#----------------------------------------------------------------------------------------------------------------
ts = TreeStyle()

ts.show_branch_support = True

# Disable the default tip names config
ts.show_leaf_name = False

ts.scale =  300 # x pixels per branch length unit

t.render("AGOs_ClassicTree.png", dpi=600, tree_style=ts)
#----------------------------------------------------------------------------------------------------------------
