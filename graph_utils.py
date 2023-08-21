from queue import PriorityQueue
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import csgraph

NO_PATH_INDICATOR = -9999

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
            if not graph[current_node, neighbor] == NO_EDGE_INDICATOR:
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
def dijkstra_scipy(graph, start_node):
    sparse_edges = csr_matrix(graph)
    return csgraph.dijkstra(csgraph=sparse_edges, 
                            directed=True, 
                            return_predecessors=True)

def floyd_warshall(graph):
    sparse_graph = csr_matrix(graph)
    _, predecessors = csgraph.floyd_warshall(csgraph=sparse_graph, directed=True, return_predecessors=True)
    return predecessors

NO_EDGE_INDICATOR = 0

def floyd_warshall_py(graph, hubs):
    n = len(graph)
    predecesors = _initialize_graph_with_num(n, NO_PATH_INDICATOR)
    for i in range(n):
        for j in range(n): 
            if graph[i][j] != NO_EDGE_INDICATOR or i==j:
                predecesors[i,j] = i

    # print("predecesor initialized")
    # print(predecesors)
    for k in hubs:
         for i in range(n):
            for j in range(n):
                if graph[i][k] == NO_EDGE_INDICATOR or graph[k][j] == NO_EDGE_INDICATOR:
                    continue

                if graph[i][j] == NO_EDGE_INDICATOR or graph[i][k] + graph[k][j] < graph[i][j]:
                    graph[i][j] = graph[i][k] + graph[k][j]
                    predecesors[i,j] = predecesors[k,j]
    for i in range(n):
        predecesors[i,i] = NO_PATH_INDICATOR
    # print("costs")
    # print(graph)
    return predecesors

def _initialize_graph(num_of_vertices):
    return _initialize_graph_with_num(num_of_vertices, NO_EDGE_INDICATOR)

def _initialize_graph_with_num(num_of_vertices, num):
    graph = np.zeros((num_of_vertices,num_of_vertices), dtype=int)
    for i in range(num_of_vertices):
        for j in range(num_of_vertices):
            graph[i,j] = num
    return graph

def _add_edge(g, u, v, weight):
    g[u][v] = weight
    g[v][u] = weight


