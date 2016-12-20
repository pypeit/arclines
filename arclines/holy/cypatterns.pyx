import numpy as np
cimport numpy as np
cimport cython
DTYPE = np.float64
ctypedef np.float_t DTYPE_t
ITYPE = np.int64
ctypedef np.int_t ITYPE_t

cdef extern from "math.h":
    double csqrt "sqrt" (double)
    double cexp "exp" (double)
    double clog "log" (double)
    double cpow "pow" (double, double)

def triangles(np.ndarray[DTYPE_t, ndim=1] detlines not None,
              np.ndarray[DTYPE_t, ndim=1] linelist not None,
              int numsrch):

    cdef int nptn = 3  # Number of lines used to create a pattern
    cdef int l, ll, sz_l, sz_d, x, xx
    cdef int cnt, nup

    sz_d = detlines.shape[0]
    sz_l = linelist.shape[0]

    # Count the number of linelist patterns that will be created
    cnt=0
    for l in range(0, sz_l-nptn+1):
        nup = (l+nptn-1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l+nptn-1, nup):
            for x in range(l+1,ll-1):
                cnt += 1
    print cnt
    return cnt, cnt

    cdef np.ndarray[ITYPE_t, ndim=2] index = np.zeros((cnt,nptn), dtype=ITYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] pattern = np.zeros((cnt,nptn-2), dtype=DTYPE)

    # Generate the patterns
    cnt=0
    for l in range(0, sz_l-nptn+1):
        nup = (l+nptn-1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l+nptn-1, nup):
            if (linelist[ll]-linelist[l]) > maxlin: continue
            # Create a pattern with these two endpoints
            for x in range(l+1,ll-2):
                for xx in range(x+1,ll-1):
                    index[cnt,0] = l
                    index[cnt,1] = x
                    index[cnt,2] = xx
                    index[cnt,3] = ll
                    pattern[cnt,0] = (linelist[x] -linelist[l])/(linelist[ll]-linelist[l])
                    pattern[cnt,1] = (linelist[xx]-linelist[l])/(linelist[ll]-linelist[l])
                    cnt += 1
    return pattern, index
