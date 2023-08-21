from cffi import FFI

ffibuilder = FFI()
ffibuilder.cdef("double normal_paths_calculation(int n, double** demand, double** cost_graph, int** predecessors);")
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

from py_vns.lib import normal_paths_calculation

def normal_paths_calulation_c(n, demand, cost_graph, predecessors):
    return normal_paths_calculation(n,
                               _cast_matrix_double(demand, ffibuilder),
                               _cast_matrix_double(cost_graph, ffibuilder),
                               _cast_matrix_int(predecessors, ffibuilder))
