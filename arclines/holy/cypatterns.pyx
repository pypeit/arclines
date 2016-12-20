import numpy as np
cimport numpy as np
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
              int detsrch, int lstsrch, double tol):
    """
    detlines - list of detected lines in pixels (sorted, increasing)
    linelist - list of lines that should be detected (sorted, increasing)
    detsrch  - Number of consecutive elements in detlines to use to create a pattern (-1 means all lines in detlines)
    lstsrch  - Number of consecutive elements in linelist to use to create a pattern (-1 means all lines in detlines)
    tol      - A tolerance that is used to determine if a match is successful
    """

    cdef int nptn = 3  # Number of lines used to create a pattern
    cdef int d, dd, sz_d, xd, cntdet, dup
    cdef int l, ll, sz_l, xl, cntlst, lup
    cdef double dval, lval, tst

    sz_d = detlines.shape[0]
    sz_l = linelist.shape[0]

    # Count the number of detlines patterns that will be created
    cntdet = 0
    for d in range(0, sz_d-nptn+1):
        dup = (d+nptn-1) + detsrch
        if dup > sz_d:
            dup = sz_d
        if detsrch == -1:
            dup = sz_d
        for dd in range(d+nptn-1, dup):
            for xd in range(d+1, dd-1):
                cntdet += 1

    # Count the number of linelist patterns that will be created
    cntlst = 0
    for l in range(0, sz_l-nptn+1):
        lup = (l+nptn-1) + lstsrch
        if lup > sz_l:
            lup = sz_l
        if lstsrch == -1:
            lup = sz_l
        for ll in range(l+nptn-1, lup):
            for xl in range(l+1, ll-1):
                cntlst += 1
    print 1 + 2 + 6*(sz_d-3), 1 + 2 + 6*(sz_l-3)
    print cntdet, cntlst, sz_d, sz_l

    cdef np.ndarray[ITYPE_t, ndim=2] lindex = np.zeros((cntdet*cntlst, nptn), dtype=ITYPE)
    cdef np.ndarray[ITYPE_t, ndim=2] dindex = np.zeros((cntdet*cntlst, nptn), dtype=ITYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] disps = np.zeros((cntdet*cntlst), dtype=DTYPE)

    # Test each detlines combination
    cntdet = 0
    for d in range(0, sz_d-nptn+1):
        dup = (l+nptn-1) + detsrch
        if dup > sz_d:
            dup = sz_d
        if detsrch == -1:
            dup = sz_d
        for dd in range(d+nptn-1, dup):
            for xd in range(d+1, dd-1):
                # Create the test pattern
                dval = (detlines[xd]-detlines[d])/(detlines[dd]-detlines[d])
                # Search through all possible patterns in the linelist
                for l in range(0, sz_l-nptn+1):
                    lup = (l+nptn-1) + lstsrch
                    if lup > sz_l:
                        lup = sz_l
                    if lstsrch == -1:
                        lup = sz_l
                    for ll in range(l+nptn-1, lup):
                        for xl in range(l+1, ll-1):
                            lval = (linelist[xl]-linelist[l])/(linelist[ll]-linelist[l])
                            tst = lval-dval
                            if tst < 0.0:
                                tst *= -1.0
                            if tst <= tol:
                                lindex[cntdet,0] = l
                                lindex[cntdet,1] = xl
                                lindex[cntdet,2] = ll
                                dindex[cntdet,0] = d
                                dindex[cntdet,1] = xd
                                dindex[cntdet,2] = dd
                                disps[cntdet] = (linelist[ll]-linelist[l]) / (detlines[dd]-detlines[d])
                            cntdet += 1
    return dindex, lindex, disps
