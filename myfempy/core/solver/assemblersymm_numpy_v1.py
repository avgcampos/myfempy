
# @profile
def getVectorizationSymm(ith_band, jth_band, val_band, ith_diag, val_diag, nb, nd, loc, matrix, ee, elemdof):

    KI = 0
    KJ = 0
    VAL = 0.0

    for ii in range(elemdof):
        KI = loc[ii]
        for jj in range(ii, elemdof):
            KJ = loc[jj]
            VAL = matrix[ii, jj]
            if KI == KJ:
                ith_diag[nd] = KI
                val_diag[nd] = VAL
                nd += 1
            else:
                ith_band[nb] = KI
                jth_band[nb] = KJ
                val_band[nb] = VAL
                nb += 1
    return ith_diag, val_diag, ith_band, jth_band, val_band, nb, nd
