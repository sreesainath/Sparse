"""
  Sparse Matrix/Tensor Library
"""

import math

# Import declarations from .h file 
cdef extern from "Sparse_matrix.h" namespace "sparse":
    cdef cppclass Sparse_matrix:
        Sparse_matrix()
        void getdata()
        void printdata()
        void getshape()
        double max()
        double min()
        double mean()
        void sort()
        void square()
        void sortindices()
        double* cumprod()
        void cumsum()
        void __pos__()
        void __neg__()
        void __abs__()
        void sq_root()
        Sparse_matrix addition(Sparse_matrix)
        double* dotprod(Sparse_matrix &)
        Sparse_matrix slices(int **)
        void clip(double, double)
        void swapaxes(int,int)
        void flatten()
     
cdef class SparseTensor:
    """ Actual python class visible to user """
    cdef Sparse_matrix *cpp_sp_mtx 

    def __init__(self):
        self.cpp_sp_mtx = new Sparse_matrix()

    def getdata(self):
        self.cpp_sp_mtx.getdata()

    def printdata(self):
        self.cpp_sp_mtx.printdata()

    def getshape(self):
        self.cpp_sp_mtx.getshape()

    def max(self):
        cdef double maxi
        maxi = self.cpp_sp_mtx.max()
        return maxi
    
    def min(self):
        cdef double mini
        mini = self.cpp_sp_mtx.min()
        return mini

    def mean(self):
        cdef double meanie
        meanie = self.cpp_sp_mtx.mean()
        return meanie
    
    def sort(self):
        self.cpp_sp_mtx.sort()

    def square(self):
        self.cpp_sp_mtx.square()
        
    def sortindices(self):
        self.cpp_sp_mtx.sortindices()

    def double* cumprod(self):
        self.cpp_sp_mtx.cumprod()
            
    def cumsum(self):
        self.cpp_sp_mtx.cumsum()
       
    def __pos__(self):
        self.cpp_sp_mtx.__pos__()

    def __neg__(self):
        self.cpp_sp_mtx.__neg__()

    def __abs__(self):
        self.cpp_sp_mtx.__abs__()

    def sq_root(self):
        self.cpp_sp_mtx.sq_root()

    def clip(self, double minval , double maxval):
        self.cpp_sp_mtx.clip(minval,maxval)

    def swapaxes(self, int axis1, int axis2):
        self.cpp_sp_mtx.swapaxes(axis1,axis2)

    def flatten(self):
        self.cpp_sp_mtx.flatten()

#    def addition(self, s):
#        self.cpp_sp_mtx.addition(s)

#    def dotprod(self, Sparse_matrix s):
#        self.cpp_sp_mtx.dotprod(s)

