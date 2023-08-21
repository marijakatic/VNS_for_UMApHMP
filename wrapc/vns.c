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
