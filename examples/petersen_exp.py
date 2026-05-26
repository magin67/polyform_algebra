import networkx as nx

from polyform_algebra.combinations.polyform import Polyform

# from polyform_algebra.objects.element import Point, Vector
# from polyform_algebra.objects.quform import quForm

# from polyform_algebra.combinations.polysimplex import Polysimplex
# from polyform_algebra.combinations.polyform import Polyform
# from polyform_algebra._core.basis import Basis


G_petersen = nx.petersen_graph()
L_petersen = Polyform.from_networkx(G_petersen)

W = L_petersen.exp()

print(W.sizes())
