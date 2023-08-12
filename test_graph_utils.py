import unittest
import graph_utils
import numpy as np

class GraphTestCase(unittest.TestCase):

    def test_givenBasicGraph_whenDijkstra_thenShortestPathsReturned(self):
        # GIVEN
        g = graph_utils._initialize_graph(9)
        graph_utils._add_edge(g, 0, 1, 4)
        graph_utils._add_edge(g, 0, 6, 7)
        graph_utils._add_edge(g, 1, 6, 11)
        graph_utils._add_edge(g, 1, 7, 20)
        graph_utils._add_edge(g, 1, 2, 9)
        graph_utils._add_edge(g, 2, 3, 6)
        graph_utils._add_edge(g, 2, 4, 2)
        graph_utils._add_edge(g, 3, 4, 10)
        graph_utils._add_edge(g, 3, 5, 5)
        graph_utils._add_edge(g, 4, 5, 15)
        graph_utils._add_edge(g, 4, 7, 1)
        graph_utils._add_edge(g, 4, 8, 5)
        graph_utils._add_edge(g, 5, 8, 12)
        graph_utils._add_edge(g, 6, 7, 1)
        graph_utils._add_edge(g, 7, 8, 3) 
        # WHEN
        D, shortest_paths = graph_utils.dijkstra(g, 0)
        # THEN
        D_expected = {0: 0, 1: 4.0, 2: 11.0, 3: 17.0, 4: 9.0, 5: 22.0, 6: 7.0, 7: 8.0, 8: 11.0}
        shortest_paths_expected = {0: [0], 
                                   1: [0, 1], 
                                   2: [0, 6, 7, 4, 2], 
                                   3: [0, 6, 7, 4, 2, 3], 
                                   4: [0, 6, 7, 4], 
                                   5: [0, 6, 7, 4, 2, 3, 5],
                                   6: [0, 6], 
                                   7: [0, 6, 7], 
                                   8: [0, 6, 7, 8]}
        self.assertEqual(D, D_expected)
        self.assertEqual(shortest_paths, shortest_paths_expected)

    # test example from https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.floyd_warshall.html
    def test_givenBasicGraph_whenFloydWarshal_thenShortestPathsReturned(self):
        # GIVEN
        g = graph_utils._initialize_graph(4)
        graph_utils._add_edge(g, 0, 1, 1)
        graph_utils._add_edge(g, 0, 2, 2)
        graph_utils._add_edge(g, 1, 3, 1)
        graph_utils._add_edge(g, 2, 0, 2)
        graph_utils._add_edge(g, 2, 3, 3)
        # WHEN
        predecessors = graph_utils.floyd_warshall(g).tolist()
        # THEN
        predecessors_expected = [[-9999,    0,     0,     1],
                                 [   1, -9999,     0,     1],
                                 [   2,     0, -9999,     2],
                                 [   1,     3,     3, -9999]]
        self.assertEqual(predecessors, predecessors_expected)


if __name__ == '__main__':
    unittest.main()