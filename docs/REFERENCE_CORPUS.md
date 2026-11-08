# Независимый эталонный корпус PMG

Этап 6. Корпус дополняет проверки API этапов 1–5; математическое ядро и прежний API не изменены.

## Источники и воспроизводимость

Зафиксирована ветка заметок `v5`, ревизия `4bea8e6e52c1299510dab05e5d96dbe6ea79d796`.
`tests/fixtures/reference_corpus_v1.json` содержит явные сети, рациональные веса строками,
ожидаемые числа, пути восьми исходных заметок, постоянные ссылки и Git blob SHA.
Двенадцать страниц PythonExamples имеют те же blob SHA, что и ранее проверенный снимок
`pythonexamples_v5.json`; их исполняемые регрессии сохранены.

Запуск: `python -m pytest tests/reference -q`; весь проект: `python -m pytest -q`.
Источники не загружаются при запуске тестов. Обновление заметок требует отдельной сверки
определений и ожидаемых ответов; автоматической перезаписи эталонов нет.

## Независимость и границы

`tests/reference/oracles.py` не импортирует PMG. Эталоны используют полный перебор
остовных лесов и разделяющих двухкомпонентных лесов, заземлённое решение SymPy,
определитель Кэли—Менгера, миноры матрицы инцидентности и формулу Якоби.
Все 20 сетей первоначальной ревизии повторно проверяются перебором лесов;
все 76 простых графов на 0–4 вершинах — через независимо построенные компоненты.
T15 проверяет все 24 порядка вершин тетраэдра, ориентацию, зависимые векторы и масштабы
в точном и численном режимах. Малый ненулевой мост допускает диагностированный отказ
численного алгоритма, но не ложное разбиение сети.

`test_note_corpus.py` отдельно проверяет экспоненту и дополнительные миноры,
граничные леса и склейку путей, изменение одного ребра, высшие грейды, структурное
тождество, симметрические формы K5 и два изомера из трёх шестиугольников.
Для каждого изомера перебираются все 65536 подмножеств рёбер; сравниваются обычные
и укоренённые леса, характеристический полином, сопротивления и спектральные моменты.
Это выбранные конечные примеры из заметок, а не запуск всего внешнего каталога `checks`
и не доказательство общих теорем. Тематические операции в публичный API не добавлены.

SymPy используется и ядром, и частью эталонов, поэтому полностью независимого
вычислительного стека нет. Перебор лесов на Fraction обеспечивает дополнительный
алгоритмически независимый путь для сопротивлений и остовных коэффициентов.
В проверке формальных коэффициентов `_pairs` читает только фактический результат;
ожидаемые коэффициенты получаются независимо из миноров инцидентности.
Это единственная намеренная зависимость корпуса от частного представления форм.

## Карта приёмки T01–T42

Машиночитаемый индекс: `tests/fixtures/reference_coverage.json`. Тест проверяет наличие
всех 42 идентификаторов и существование указанных тестов. Индекс обозначает точки входа,
а не автоматически доказывает полноту математических условий; дополнительные проверки
расположены в тех же модулях. Старые тематические тесты сохраняют совместимость.

| Требование | Проверки |
| --- | --- |
| T01 | `tests/test_core_network_math.py::test_electrical_scale_and_single_edge` |
| T02 | `tests/test_core_network_math.py::test_weighted_chain_matrices` |
| T03 | `tests/test_core_network_math.py::test_triangle_and_gram_geometry` |
| T04 | `tests/test_core_network_math.py::test_triangle_and_gram_geometry` |
| T05 | `tests/test_core_reductions.py::test_weighted_path` |
| T06 | `tests/test_core_reductions.py::test_three_reductions_differ` |
| T07 | `tests/test_core_reductions.py::test_star_and_sequential_elimination` |
| T08 | `tests/test_core_network_math.py::test_disconnected_empty_singleton_and_infinity` |
| T09 | `tests/test_core_network_math.py::test_disconnected_empty_singleton_and_infinity` |
| T10 | `tests/test_core_network_math.py::test_all_simple_graphs_up_to_four_vertices` |
| T11 | `tests/test_core_network_math.py::test_invalid_resistance`, `tests/test_core_network_math.py::test_all_simple_graphs_up_to_four_vertices` |
| T12 | `tests/test_core_network_math.py::test_exact_scale_and_permutation_laws`, `tests/test_core_data.py::test_all_permutations_preserve_labeled_network_and_coordinates` |
| T13 | `tests/test_core_network_math.py::test_exact_scale_and_permutation_laws` |
| T14 | `tests/test_core_network_math.py::test_triangle_and_gram_geometry` |
| T15 | `tests/reference/test_contract_corpus.py::test_t15_volume_orientation_basepoint_and_scale` |
| T16 | `tests/test_core_reductions.py::test_general_network_identities` |
| T17 | `tests/test_core_reductions.py::test_star_and_sequential_elimination` |
| T18 | `tests/test_core_reductions.py::test_contraction_and_induced_degrees` |
| T19 | `tests/test_core_algebra.py::test_wedge_algebra_laws` |
| T20 | `tests/test_core_algebra.py::test_boundary_squared_and_leibniz` |
| T21 | `tests/test_core_algebra.py::test_weighted_boundary_regression` |
| T22 | `tests/test_core_algebra.py::test_simple_sign_square_and_commuting_forms`, `tests/test_core_algebra.py::test_potential_metric_identity_and_nonmultiplicative_evaluation` |
| T23 | `tests/test_core_algebra.py::test_tetrahedron_identity_and_product` |
| T24 | `tests/test_core_algebra.py::test_polarization_normalization_and_cancellation` |
| T25 | `tests/test_core_algebra.py::test_exponential_factorials_and_mixed_grades` |
| T26 | `tests/test_core_algebra.py::test_formal_tree_coefficients_all_graphs_up_to_four_vertices` |
| T27 | `tests/test_core_algebra.py::test_potential_metric_identity_and_nonmultiplicative_evaluation` |
| T28 | `tests/test_core_network_math.py::test_electrical_scale_and_single_edge` |
| T29 | `tests/test_core_algebra.py::test_exterior_normalization_and_context` |
| T30 | `tests/test_core_algebra.py::test_weighted_boundary_regression` |
| T31 | `tests/test_core_reductions.py::test_numeric_huge_internal_contraction` |
| T32 | `tests/reference/test_contract_corpus.py::test_weak_bridge_is_resolved_or_reported` |
| T33 | `tests/test_core_data.py::test_network_invalid_edges`, `tests/test_core_data.py::test_context_identity_and_order_are_not_names`, `tests/test_core_data.py::test_exact_rejects_nonrational_inputs` |
| T34 | `tests/test_core_network_math.py::test_invalid_laplacians` |
| T35 | `tests/test_core_reductions.py::test_invalid_boundary`, `tests/test_core_reductions.py::test_invalid_partition` |
| T36 | `tests/test_core_algebra.py::test_unsupported_square_and_polarization`, `tests/test_core_algebra.py::test_numeric_formal_objects_are_rejected` |
| T37 | `tests/test_core_algebra.py::test_limits_cover_intermediates_and_keep_inputs` |
| T38 | `tests/test_core_data.py::test_conversion_rejects_overflow_and_underflow` |
| T39 | `tests/reference/test_contract_corpus.py::test_audit_answers_against_two_forest_enumeration` |
| T40 | `tests/reference/test_contract_corpus.py::test_all_small_graph_components_and_projectors` |
| T41 | `tests/test_core_algebra.py::test_1600_oriented_simple_products_against_minors` |
| T42 | `tests/reference/test_note_corpus.py` |
