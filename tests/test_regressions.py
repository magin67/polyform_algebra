"""Регрессии вводной серии PMG, заметки 8–10.

Источник: https://github.com/magin67/pmg-notes/tree/v5/pmg-library/ru/Введение
Граница и правило Лейбница; область поляризации; фиксированный носитель метрики.
"""
import networkx as nx
import numpy as np
import pytest
from polyform_algebra import Element, Point, Vector, Simplex, Polysimplex, Polyform, quForm


def test_scalar_unit_and_zero():
    """Скалярная единица не равна нулю: 1·a=a, 0·a=0 и ∂a=1.

    Регрессия: канонизация превращала пустое произведение (единицу) в ноль."""
    a = Point("a")
    one = Simplex.one()
    zero = Simplex.zero()
    assert one != zero
    assert Polysimplex.one().is_one()
    assert Simplex([]).canonical().is_one()
    assert (one * Simplex([a])) == Simplex([a])
    assert (zero * Simplex([a])).is_zero()
    assert a.boundary() == Polysimplex.one()


@pytest.mark.parametrize("sign", [1, -1])
def test_boundary_expansion_and_orientation(sign):
    """Проверяем (ab)=b−a, (abc)=bc−ac+ab и ∂²=0 при обеих ориентациях."""
    a, b, c = Point.create_list("abc")
    assert Simplex([(a,b)], sign=sign).to_polysimplex() == sign * (b-a)
    assert Simplex([(a,b,c)], sign=sign).to_polysimplex() == sign * (b*c-a*c+a*b)
    assert Simplex([a,b], sign=sign).boundary() == sign * (b-a)
    assert Simplex([a,b,c], sign=sign).boundary().boundary().is_zero()


def test_expansion_preserves_vectors_and_products():
    """Раскрытие границ сохраняет произведение и не обнуляет атомарный вектор.

    Нулевая граница вектора не означает, что сам вектор равен нулю."""
    a,b,c,d = Point.create_list("abcd")
    v = Vector("v")
    assert Simplex([v]).to_polysimplex() == v.as_polysimplex()
    assert Simplex([(a,b),(c,d)]).to_polysimplex() == (b-a)*(d-c)
    assert Simplex([a,(b,c)]).to_polysimplex() == a*(c-b)


def test_boundary_leibniz_and_affine_weight():
    """Проверяем ∂(XY)=(∂X)Y+(−1)^grade(X)X∂Y и ∂e=ε(e).

    Наличие границы среди множителей не обнуляет границу всего произведения."""
    a,b,c,d = Point.create_list("abcd")
    edge = Simplex([(c,d)])
    assert Simplex([a,(c,d)]).boundary() == edge.as_polysimplex()
    assert Simplex([a,b,(c,d)]).boundary() == (b-a)*edge
    assert Simplex([(a,b),(c,d)]).boundary().is_zero()
    e = Element("e", multiplicity=2)
    assert Simplex([e,a]).boundary() == 2*a-e


def test_compact_first_grade_polarization():
    """Поляризация не зависит от записи вектора: (ab) и b−a дают один результат."""
    a,b,c = Point.create_list("abc")
    assert Simplex([(a,b)]) @ Simplex([(a,c)]) == (b-a) @ (c-a)


@pytest.mark.parametrize("left,right", [(2,2),(1,2),(0,0)])
def test_unsupported_scalar_product_grades_raise(left,right):
    """Фиксируем ограничение текущего API: @ поддерживает только первый грейд.

    Это проверка явного отказа, а не математического запрета других грейдов."""
    points = Point.create_list("abc")
    x = Simplex(points[:left])
    y = Simplex(points[:right])
    with pytest.raises(ValueError, match="grade 1"):
        x @ y


def test_isolated_vertices_survive_matrix_and_algebra():
    """Вершина 2 изолирована, но должна сохранять строку и столбец лапласиана.

    Даже обнуление коэффициентов не должно стирать исходный набор вершин."""
    graph = nx.Graph()
    graph.add_nodes_from([0,1,2])
    graph.add_edge(0,1)
    lap, points = Polyform.from_networkx(graph, return_points="list")
    expected = np.array([[1,-1,0],[-1,1,0],[0,0,0]])
    for form in [lap, Polyform(lap), lap+Polyform.zero(), 2*lap/2, -(-lap), lap*Polyform.one(), lap.extract_grade(1)]:
        matrix,index = form.to_matrix()
        assert list(index) == points
        np.testing.assert_allclose(matrix,expected)
    zero = lap-lap
    assert zero.to_matrix()[0].shape == (3,3)
    assert (0*lap).to_matrix()[0].shape == (3,3)
    assert len(lap.exp().vertices) == 3


def test_edgeless_graph_matrix():
    """Для трёх изолированных вершин L=0 и exp(L)=1, но остовов нет.

    Поэтому матрица остаётся 3×3, а нормировка на остовной коэффициент запрещена."""
    lap = Polyform.from_networkx(nx.empty_graph(3))
    np.testing.assert_array_equal(lap.to_matrix()[0],np.zeros((3,3)))
    assert lap.exp() == Polyform.one()
    with pytest.raises(ValueError, match="spanning"):
        lap.exp().norm(quForm.one())


def test_laplacian_round_trip_with_isolated_vertex():
    """Переход матрица → полиформа → матрица сохраняет нулевую строку изолированной вершины."""
    points=Point.create_list("abc")
    matrix=np.array([[1,-1,0],[-1,1,0],[0,0,0]])
    lap=Polyform.from_laplacian(matrix,points)
    np.testing.assert_allclose(lap.to_matrix()[0],matrix)


def test_disconnected_metric_rejects_normalization_but_has_potential():
    """Два отдельных ребра не имеют остова: τ(M)=0, но τ(M·[(12)]²)=1.

    Связь 1–2 соединяет компоненты в дерево. Потенциал определён, деление на τ(M) — нет."""
    graph=nx.Graph([(0,1),(2,3)])
    lap,p=Polyform.from_networkx(graph,return_points="list")
    metric=lap.exp()
    assert metric.potential(quForm.one()) == 0
    assert metric.potential(quForm([(p[1],p[2])])) == 1
    for query in [quForm.one(),quForm([(p[0],p[2])]),Polyform.zero()]:
        with pytest.raises(ValueError,match="spanning"):
            metric.norm(query)


def test_connected_metrics_and_edge_deletion():
    """У пути 0–1–2 сопротивление концов равно 2; удаление ребра разрушает связность.

    Дополнительно проверяем линейность потенциала и запрет измерений вне набора вершин."""
    lap,p=Polyform.from_networkx(nx.path_graph(3),return_points="list")
    metric=lap.exp()
    assert metric.norm(quForm([(p[0],p[2])])) == 2
    # Линейность: потенциал единицы равен 1, формы концов пути — 2.
    assert metric.potential(Polyform.one()+quForm([(p[0],p[2])])) == 3
    with pytest.raises(ValueError,match="outside"):
        metric.norm(quForm([(p[0],Point("outside"))]))
    removed=lap-quForm([(p[1],p[2])])
    with pytest.raises(ValueError,match="spanning"):
        removed.exp().norm(quForm([(p[0],p[2])]))


def test_triangle_metric_and_higher_quadratic_form():
    """Эталон K3: сопротивление 2/3, скалярное произведение 1/3, норма границы 1/3.

    Последняя равна det([[2/3,1/3],[1/3,2/3]]); это квадрат удвоенной площади."""
    lap,p=Polyform.from_networkx(nx.complete_graph(3),return_points="list")
    metric=lap.exp()
    assert metric.norm(quForm([(p[0],p[1])])) == pytest.approx(2/3)
    assert metric.norm(quForm([tuple(p)])) == pytest.approx(1/3)
    assert metric.norm((p[1]-p[0])@(p[2]-p[0])) == pytest.approx(1/3)


def test_quadratic_zero_is_not_unit():
    """Нулевая форма и единичная форма различаются и как объекты, и как ключи множества."""
    assert quForm.zero() != quForm.one()
    assert len({quForm.zero(),quForm.one()}) == 2
    assert Polyform.zero() != Polyform.one()


def test_context_survives_powers_and_connected_variation():
    """Степени нулевой полиформы сохраняют вершины; добавление двух рёбер создаёт путь.

    Проверяем также, что построение нового графа не изменяет исходную полиформу."""
    lap,p=Polyform.from_networkx(nx.empty_graph(3),return_points="list")
    assert (lap**0).vertices == lap.vertices
    assert (lap**2).vertices == lap.vertices
    connected = lap + quForm([(p[0],p[1])]) + quForm([(p[1],p[2])])
    assert connected.exp().norm(quForm([(p[0],p[2])])) == 2
    assert lap.is_zero()


def test_compact_higher_grade_polarization_is_not_silently_zero():
    """Неподдерживаемая поляризация границы треугольника должна сообщить об ограничении.

    Раньше X@X молча давало ноль, хотя метрическая оценка в K3 равна 1/3."""
    a,b,c=Point.create_list("abc")
    x=Simplex([(a,b,c)])
    with pytest.raises(ValueError,match="grade 1"):
        x @ x


def test_oriented_products_and_leibniz_against_expansion():
    """Сверяем компактные границы с раскрытием через вершины для 1600 пар.

    Перестановки вершин проверяют ориентацию; раскрытие проверяет слияние границ."""
    from itertools import permutations, product
    points=Point.create_list("abcd")
    # 4 точки + 12 ориентированных пар + 24 ориентированные тройки = 40 объектов.
    factors=[Simplex([x]) for x in points]
    factors += [Simplex([x]) for k in (2,3) for x in permutations(points,k)]

    def expand(obj):
        return sum((c*t.to_polysimplex() for t,c in obj.terms.items()),Polysimplex.zero())

    for x,y in product(factors,repeat=2):
        assert (x*y).to_polysimplex() == x.to_polysimplex()*y.to_polysimplex(), (x,y)
        expected=x.boundary()*y+((-1)**x.grade())*(x*y.boundary())
        assert expand((x*y).boundary()) == expand(expected), (x,y)


@pytest.mark.parametrize("graph", [nx.path_graph(5),nx.cycle_graph(4),nx.complete_graph(4)])
def test_graph_metric_matches_independent_matrix_computation(graph):
    """Сверяем алгебру с независимым матричным расчётом на взвешенных P5, C4, K4.

    По Кирхгофу τ — минор лапласиана; сопротивление R_ij=Gr_ii+Gr_jj−2Gr_ij."""
    for i,(u,v) in enumerate(graph.edges()):
        graph[u][v]["weight"]=i+1
    lap,points=Polyform.from_networkx(graph,return_points="list")
    matrix,_=lap.to_matrix()
    expected=np.zeros_like(matrix)
    for u,v,data in graph.edges(data=True):
        w=data["weight"]
        expected[u,u]+=w
        expected[v,v]+=w
        expected[u,v]-=w
        expected[v,u]-=w
    np.testing.assert_allclose(matrix,expected)
    # Псевдообратная даёт матрицу Грина независимо от операций с полиформами.
    green=np.linalg.pinv(expected)
    metric=lap.exp()
    assert metric.potential(quForm.one()) == pytest.approx(np.linalg.det(expected[1:,1:]))
    for i in range(len(points)):
        for j in range(i+1,len(points)):
            resistance=green[i,i]+green[j,j]-2*green[i,j]
            assert metric.norm(quForm([(points[i],points[j])])) == pytest.approx(resistance)
