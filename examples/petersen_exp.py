import networkx as nx

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.combinations.polyform import Polyform

G_petersen = nx.petersen_graph()
G = G_petersen
poly = Polyform.from_networkx(G)
W = poly.exp()

print(W.sizes())
