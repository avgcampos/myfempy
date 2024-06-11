# import cython
import numpy as np

INT32 = np.int32
FLT64 = np.float64

from myfempy.core.solver.assembler import getLoc, getMatrix


# @profile
def getMatrixAssemblerSymm(
    Model, inci, coord, tabmat, tabgeo, elemdof, intgauss, type_assembler
):
    # ith_band = []
    # jth_band = []
    # val_band = []
    # ith_diag = []
    # jth_diag = []
    # val_diag = []

    dim_band = int(0.5 * (elemdof * elemdof - elemdof) * inci.shape[0])
    dim_diag = int(elemdof * inci.shape[0])

    ith_band = np.zeros((dim_band,), dtype=INT32)
    jth_band = np.zeros((dim_band,), dtype=INT32)
    val_band = np.zeros((dim_band,), dtype=FLT64)
    ith_diag = np.zeros((dim_diag,), dtype=INT32)
    val_diag = np.zeros((dim_diag,), dtype=FLT64)

    KI = 0
    KJ = 0
    VAL = 0.0

    nb = 0
    nd = 0
    for ee in range(inci.shape[0]):
        mat = getMatrix(
            Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler
        )
        loc = getLoc(Model, inci, ee)
        for ii in range(elemdof):
            KI = loc[ii]
            for jj in range(ii, elemdof):
                KJ = loc[jj]
                VAL = mat[ii, jj]
                if KI == KJ:
                    ith_diag[nd] = KI
                    val_diag[nd] = VAL
                    nd += 1
                else:
                    ith_band[nb] = KI
                    jth_band[nb] = KJ
                    val_band[nb] = VAL
                    nb += 1
    return ith_diag, val_diag, ith_band, jth_band, val_band
