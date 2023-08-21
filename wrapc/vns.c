#include "vns.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#define NO_EDGE_INDICATOR 0
#define NO_PATH_INDICATOR -9999

/*
python implementation:
NO_PATH_INDICATOR = -9999
def _normal_paths_calulation(n, demand, cost_graph, predecessors):
    tmp_total_cost = 0
    for orig in range(n):
        for dest in range(n):
            left = predecessors[orig, dest]
            right = dest
            while True:
                if left == NO_PATH_INDICATOR:
                    break
                tmp_total_cost += demand[orig, dest] * cost_graph[left, right]
                if left == orig:
                    break
                right = left
                left = predecessors[orig, left]
    return tmp_total_cost
*/
// #define NO_PATH_INDICATOR -9999
double normal_paths_calculation(int n, double** demand, double** cost_graph, int** predecessors) {
    double tmp_total_cost = 0;
    for(int orig = 0; orig < n; orig++) {
        for(int dest = 0; dest < n; dest++) {
            int left = predecessors[orig][dest];
            int right = dest;
            while (1) {
                if (left == NO_PATH_INDICATOR) {
                    break;
                }
                tmp_total_cost += demand[orig][dest] * cost_graph[left][right];
                if (left == orig) {
                    break;
                }
                right = left;
                left = predecessors[orig][left];
            }
        }
    }
    return tmp_total_cost;
}

/*
--------python implementation:----------
NO_EDGE_INDICATOR = 0
NO_PATH_INDICATOR = -9999

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
*/


int** floyd_warshall_c_impl(double** graph, int* hubs, int n, int p) {
    int** predecessors;
    // allocate predecessors
    predecessors = (int **)malloc(n * sizeof(int*));
    for(int i = 0; i < n; i++) predecessors[i] = (int *)malloc(n * sizeof(int));
    // initialize predecessors
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            predecessors[i][j] = NO_PATH_INDICATOR;
        }
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (graph[i][j] != NO_EDGE_INDICATOR || i==j) {
                predecessors[i][j] = i;
            }
        }
    }

    for(int k = 0; k < p; k++){
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                int hub = hubs[k];
                // this is not needed, is it? there will always be edge, considering how we built the graph
                if (graph[i][hub] == NO_EDGE_INDICATOR || graph[hub][j] == NO_EDGE_INDICATOR) {
                    continue;
                }
                if (graph[i][j] == NO_EDGE_INDICATOR || graph[i][hub] + graph[hub][j] < graph[i][j]) {
                    graph[i][j] = graph[i][hub] + graph[hub][j];
                    predecessors[i][j] = predecessors[hub][j];   
                }
            }
        }
    }

    for(int i = 0; i < n; i++) {
        predecessors[i][i] = NO_PATH_INDICATOR;
    }

    return predecessors;
}


int isHub(int node, int* hub_bitmap) {
    return hub_bitmap[node];
}

/*
def bitmap(n, list):
    return [1 if i in list else 0 for i in range(n)]
*/

int* bitmap(int n, int* arr, int p) {
    int* bitmap = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        bitmap[i] = 0;
    }
    for (int k = 0; k < p; k++) {
        bitmap[arr[k]] = 1;
    }
    return bitmap;
}

/*
def _closest(node, hubs, cost_graph):
    hub_dist = {hub: cost_graph[node, hub] for hub in hubs}
    return min(hub_dist, key=hub_dist.get)
*/

int closest_hub_cost(int node, int* hubs, int p, double** cost_graph, double** demand, int n) {
    double min_cost = DBL_MAX;
    for (int k = 0; k < p; k++) {
        double cost = (cost_graph[node][hubs[k]] + cost_graph[hubs[k]][node]);
        if(cost < min_cost) {
            min_cost = cost;
        }
    }
    return  demand[node][node] * min_cost;
}

/*
def _peculiar_paths_calulations(n, problem, hubs, cost_graph, discounts):
    tmp_total_cost = 0
    for node in range(n):
        if node in hubs:
            continue
        closest_hub = _closest(node, hubs, cost_graph)
        # adding cost of [node -> closest_hub -> node] path
        tmp_total_cost += problem.demand[node, node] * problem.distances[node, closest_hub] * discounts[node, closest_hub] + \
                      problem.demand[node, node] * problem.distances[closest_hub, node] * discounts[closest_hub, node]
    return tmp_total_cost
*/

double peculiar_paths_calulation(int n, int p, double** demand, int* hubs, int* hub_bitmap, double** cost_graph) {
    double total_cost = 0;
    for (int node = 0; node < n; node++) {
        // if node is hub-node, continue
        if (isHub(node, hub_bitmap)) {
            continue;
        }
        total_cost +=  closest_hub_cost(node, hubs, p, cost_graph, demand, n);
    }
    return total_cost;
}


/*
def _get_valid_cost_graph(n, hubs, distances, discounts):
    nodes = get_nodes(n)
    cost_graph = distances*discounts
    for A in nodes:
        for B in nodes:
            # connection between non-hub nodes is not allowed
            if A != B and (A not in hubs) and (B not in hubs):
                cost_graph[A, B] = NO_EDGE_INDICATOR
    return cost_graph
*/

double** get_cost_graph(int* hub_bitmap, int n, int p, double** distances, double** discounts) {
    // initialize cost_graph
    double** cost_graph = (double **)malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++) cost_graph[i] = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cost_graph[i][j] = distances[i][j] * discounts[i][j];
        }
    }

    for (int node1 = 0; node1 < n; node1++) {
        for (int node2 = 0; node2 < n; node2++) {
            // if i and j are both non-hub nodes, connection is not allowed
            if (node1 != node2 && !isHub(node1, hub_bitmap) && !isHub(node2, hub_bitmap)) {
                cost_graph[node1][node2] = NO_EDGE_INDICATOR;
            }
        }
    }
    return cost_graph;
}

/*
def get_discount_matrix(n, hubs, alpha, delta, ksi):
    discounts = np.ones((n,n))
    nodes = get_nodes(n)
    for A in nodes:
        for B in nodes:
            if (A in hubs) and (B in hubs):
                discounts[A, B] = alpha
            if (A not in hubs) and (B in hubs):
                discounts[A, B] = ksi
            if (A in hubs) and (B not in hubs):
                discounts[A, B] = delta
    return discounts
*/

double** get_discount_matrix(int n, int* hub_bitmap, double alpha, double delta, double ksi) {
    double** discounts = (double **)malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++) discounts[i] = (double *)malloc(n * sizeof(double));
    // initialize
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            discounts[i][j] = 1;
        }
    }

    for (int node1 = 0; node1 < n; node1++) {
        for (int node2 = 0; node2 < n; node2++) {
            if (isHub(node1, hub_bitmap) && isHub(node2, hub_bitmap)) {
                discounts[node1][node2] = alpha;
            }
            if (!isHub(node1, hub_bitmap) && isHub(node2, hub_bitmap)) {
                discounts[node1][node2] = ksi;
            }
            if (isHub(node1, hub_bitmap) && !isHub(node2, hub_bitmap)) {
                discounts[node1][node2] = delta;
            }
        }
    }
    return discounts;
}


/*
def _prepare(hubs, problem):
    nodes = get_nodes(problem.n)
    discounts = get_discount_matrix(problem.n, hubs, problem.alpha, problem.delta, problem.ksi)
    # prepare graph
    cost_graph = _get_valid_cost_graph(problem.n, hubs, problem.distances, discounts)
    # calculate cost
    return nodes, discounts, cost_graph
*/

/*
def get_solution_cost_fw(hubs, problem, use_c=False):
    total_cost = 0
    nodes, discounts, cost_graph = _prepare(hubs, problem)
    # calculate cost
    predecessors = floyd_warshall_c(cost_graph, hubs, problem.n, problem.p)

    if use_c:
        total_cost += normal_paths_calulation_c(problem.n, problem.demand, cost_graph, predecessors)
    else:
        total_cost += _normal_paths_calulation(problem.n, problem.demand, cost_graph, predecessors)

    # origin equals destination paths
    total_cost += _peculiar_paths_calulations(problem.n, problem, hubs, cost_graph, discounts)

    return total_cost
*/

void free_matrix_double(double** matrix, int n) {
    for(int i = 0; i < n; i++)
        free(matrix[i]);
    free(matrix);
}

void free_matrix_int(int** matrix, int n) {
    for(int i = 0; i < n; i++)
        free(matrix[i]);
    free(matrix);
}

double get_solution_cost_c_impl(int* hubs, int n, int p, double** distances, double** demand, double alpha, double delta, double ksi) {
    double total_cost = 0;
    // prepare
    int* hub_bitmap = bitmap(n, hubs, p);
    double** discounts = get_discount_matrix(n, hub_bitmap, alpha, delta, ksi);
    double** cost_graph = get_cost_graph(hub_bitmap, n, p, distances, discounts);

    int** predecessors =  floyd_warshall_c_impl(cost_graph, hubs, n, p);

    total_cost += normal_paths_calculation(n, demand, cost_graph, predecessors);
    total_cost += peculiar_paths_calulation(n, p, demand, hubs, hub_bitmap, cost_graph);

    free(hub_bitmap);
    free_matrix_double(discounts, n);
    free_matrix_double(cost_graph, n);
    free_matrix_int(predecessors, n);
    // free(distances);
    // free(demand);

    return total_cost;
}
