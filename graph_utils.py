from queue import PriorityQueue
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import csgraph
    
def dijkstra(graph, start_node):
    n = len(graph)
    D = {v : float('inf') for v in range(n)}
    shortest_path = {node : [start_node] for node in range(n)}
    D[start_node] = 0

    visited = []
    pq = PriorityQueue()
    pq.put((0, start_node))

    while not pq.empty():
        (dist, current_node) = pq.get()
        visited.append(current_node)

        for neighbor in range(n):
            if graph[current_node, neighbor] != 0:
                distance = graph[current_node, neighbor]
                if neighbor not in visited:
                    old_cost = D[neighbor]
                    new_cost = D[current_node] + distance
                    if new_cost < old_cost:
                        pq.put((new_cost, neighbor))
                        D[neighbor] = new_cost
                        shortest_path[neighbor] = shortest_path[current_node] + [neighbor]
    return D, shortest_path

# todo adapt output
def scipy_dijkstra(graph, start_node):
    sparse_edges = csr_matrix(graph)
    return csgraph.dijkstra(csgraph=sparse_edges, 
                            directed=True, 
                            return_predecessors=True)


def _initialize_graph(num_of_vertices):
    return np.zeros((num_of_vertices,num_of_vertices))

def _add_edge(g, u, v, weight):
    g[u][v] = weight
    g[v][u] = weight


