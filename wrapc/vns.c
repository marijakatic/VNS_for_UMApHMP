#include "vns.h"
#include <stdlib.h>
#include <stdio.h>

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
    // move this into macro 
    int NO_PATH_INDICATOR = -9999;
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
def _initialize_graph(num_of_vertices):
    return _initialize_graph_with_num(num_of_vertices, NO_EDGE_INDICATOR)

def _initialize_graph_with_num(num_of_vertices, num):
    graph = np.zeros((num_of_vertices,num_of_vertices), dtype=int)
    for i in range(num_of_vertices):
        for j in range(num_of_vertices):
            graph[i,j] = num
    return graph

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
    // move these into macro
    int NO_EDGE_INDICATOR = 0;
    int NO_PATH_INDICATOR = -9999;
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
    // printf("%d %d\n", n, p);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (graph[i][j] != NO_EDGE_INDICATOR || i==j) {
                predecessors[i][j] = i;
            }
        }
    }

    // for(int i = 0; i < n; i++) {
    //     for(int j = 0; j < n; j++) {
    //        printf("%d ", predecessors[i][j]);
    //     }
    //     printf("\n");
    // }

    // printf("hubs\n");
    // for(int k = 0; k < p; k++){
    //     printf("%d ", hubs[k]);
    // }
    // printf("\n");

    for(int k = 0; k < p; k++){
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                int hub = hubs[k];
                // if(i == 0 && j == 3) {
                //     printf("slucaj 0,3: hub=%d, graph[i][j]=%f, graph[i][hub]=%f, graph[hub][j]=%f\n", hub, graph[i][j], graph[i][hub], graph[hub][j]);
                // }
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