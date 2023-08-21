from cffi import FFI
import numpy as np

ffibuilder = FFI()
ffibuilder.cdef("double normal_paths_calculation(int n, double** demand, double** cost_graph, int** predecessors);")
ffibuilder.cdef("int** floyd_warshall_c_impl(double** graph, int* hubs, int n, int p);")
ffibuilder.cdef("double get_solution_cost_c_impl(int* hubs, int n, int p, double** distances, double** demand, double alpha, double delta, double ksi);")

ffibuilder.set_source("py_vns",'#include "wrapc/vns.h"',
                        sources=["wrapc/vns.c"])
ffibuilder.compile()

def _cast_matrix_double(matrix, ffi):
    ap = ffi.new("double* [%d]" % (matrix.shape[0]))
    ptr = ffi.cast("double *", matrix.ctypes.data)
    for i in range(matrix.shape[0]):
        ap[i] = ptr + i*matrix.shape[1]
    return ap

def _cast_matrix_int(matrix, ffi):
    ap = ffi.new("int* [%d]" % (matrix.shape[0]))
    ptr = ffi.cast("int *", matrix.ctypes.data)
    for i in range(matrix.shape[0]):
        ap[i] = ptr + i*matrix.shape[1]
    return ap

# ------------------------------------------
#        wrap c implementations:
# ------------------------------------------


from py_vns.lib import normal_paths_calculation
from py_vns.lib import floyd_warshall_c_impl
from py_vns.lib import get_solution_cost_c_impl


def normal_paths_calulation_c(n, demand, cost_graph, predecessors):
    return normal_paths_calculation(n,
                               _cast_matrix_double(demand, ffibuilder),
                               _cast_matrix_double(cost_graph, ffibuilder),
                               _cast_matrix_int(predecessors, ffibuilder))

def floyd_warshall_c(graph, hubs, n, p):
    ret = floyd_warshall_c_impl(_cast_matrix_double(graph, ffibuilder),
                                hubs,
                                n,
                                p)
    return np.array([[ret[i][j] for j in range(n)] for i in range(n)])

def get_solution_cost_c(hubs, problem):
    return get_solution_cost_c_impl(hubs,
                                    problem.n,
                                    problem.p,
                                    _cast_matrix_double(problem.distances, ffibuilder),
                                    _cast_matrix_double(problem.demand, ffibuilder),
                                    problem.alpha,
                                    problem.delta,
                                    problem.ksi)
