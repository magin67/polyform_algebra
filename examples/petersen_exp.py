import networkx as nx

from polyform_algebra.combinations.polyform import Polyform

G_petersen = nx.petersen_graph()
G = G_petersen
poly = Polyform.from_networkx(G)
W = poly.exp()

print(W.sizes())
