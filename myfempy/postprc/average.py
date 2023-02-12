#!/usr/bin/env python
"""
Average Nodes Calculator version 2
"""
import Cython
import numpy as np
import scipy.sparse as sp


def elem2nodes_filter(nnode: cython.int, nelem: cython.int, dofe: cython.int, inci):
    """_summary_

    Returns:
        _description_
    """
    # ith: cython.int[nelem * (dofe * dofe)]
    # jth: cython.int[nelem * (dofe * dofe)]
    # val: cython.double[nelem * (dofe * dofe)]

    ith = np.zeros((nelem * (dofe * dofe)), dtype=int)
    jth = np.zeros((nelem * (dofe * dofe)), dtype=int)
    val = np.zeros((nelem * (dofe * dofe)), dtype=int)
    q0 = 0
    for i in range(nnode):
        elmlist = inci[(np.asarray(np.where(inci[:, 4:] == i + 1)))[0][:], 0]
        q1 = elmlist.size
        ith[q0 : q1 + q0] = i
        jth[q0 : q1 + q0] = elmlist - 1
        val[q0 : q1 + q0] = elmlist
        q0 = q1 + q0
    S = sp.csc_matrix((val, (ith, jth)), shape=(nnode, nelem))
    return S

def results_average(results_elm, nnode: cython.int, nelem: cython.int, dofe: cython.int, inci):
    """_summary_

    Arguments:
        results_elm -- _description_

    Returns:
        _description_
    """
    # results_avr: cython.double[nnode]

    S = elem2nodes_filter(nnode, nelem, dofe, inci)
    results_avr = np.zeros((nnode), dtype=float)

    results_avr = [np.mean(results_elm[(S[mm, :].nonzero())[1]]) for mm in range(nnode)]

    # for mm in range(nnode):
    #     results_avr[mm] = np.mean(results_elm[(S[mm, :].nonzero())[1]])
    return results_avr