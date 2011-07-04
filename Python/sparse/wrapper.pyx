"""
  Sparse Matrix/Tensor Library
"""

import math

# Import declarations from .h file 
cdef extern from "Test.h" namespace "test":
    cdef cppclass SomeClass:
        SomeClass()
        int method1()
        

class SparseTensor:
    """ Actual python class visible to user """
    cdef SparseMatrix *cpp_sp_mtx 

    def __init__(self):
        self.cpp_sp_mtx = new SomeClass()

    def max(self):
       self.cpp_sp_mtx.max()

    def dotproduct(self)


        #self.some = new SomeClass()
        #cdef SomeClass *some = new SomeClass()
        #del some

#
#some.method1()
#del test
