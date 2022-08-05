#!/usr/bin/env python
__doc__ = """
LEGACY MESH GEN
"""
import numpy as np


def get_legacy_line2(GEOMETRY):
    nel = GEOMETRY["nx"]
    lx = GEOMETRY["lx"]
    step = lx / nel
    conec = np.zeros((nel, 3))
    coord = np.zeros((nel + 1, 4))
    coord[:, 0] = np.arange(1, nel + 2, 1).tolist()
    coord[:, 1] = np.arange(0, lx + step, step).tolist()
    conec[:, 0] = np.arange(1, nel + 1, 1).tolist()
    conec[:, 1] = np.arange(1, nel + 1, 1).tolist()
    conec[:, 2] = np.arange(2, nel + 2, 1).tolist()
    return conec, coord


def get_legacy_tria3(GEOMETRY):
    nelx = GEOMETRY["nx"]
    nely = GEOMETRY["ny"]
    nel = nelx * nely * 2
    nnx = nelx + 1
    nny = nely + 1
    nos = nnx * nny
    lx = GEOMETRY["lx"]
    ly = GEOMETRY["ly"]
    conec = np.zeros((nel, 4))
    for i in range(1, nel, 2):
        linha = int(np.ceil(i / (2 * nelx)))
        y = 2 * linha - 1
        n1 = (i + y) / 2
        n2 = n1 + 1
        n3 = n2 + nelx
        conec[i - 1, 0] = i
        conec[i - 1, 1] = n1
        conec[i - 1, 2] = n2
        conec[i - 1, 3] = n3
    for i in range(2, nel + 1, 2):
        linha = int(np.ceil(i / (2 * nelx)))
        y = 2 * linha
        n1 = (i + y) / 2
        n2 = n1 + nelx
        n3 = n2 + 1
        conec[i - 1, 0] = i
        conec[i - 1, 1] = n1
        conec[i - 1, 2] = n2
        conec[i - 1, 3] = n3
    nos = int(conec[nel - 1, 3])
    coord = np.zeros((nos, 4))
    coord[0, 0] = 1
    for i in range(2, nos + 1):
        linha = int(np.ceil(i / (nelx + 1))) - 1
        coord[i - 1, 0] = i
        coord[i - 1, 1] = ((i - 1) - linha * (nelx + 1)) * (lx / nelx)
        coord[i - 1, 2] = linha * (ly / nely)
    return conec, coord


def get_legacy_quad4(GEOMETRY):
    nelx = GEOMETRY["nx"]
    nely = GEOMETRY["ny"]
    nel = nelx * nely
    nnx = nelx + 1
    nny = nely + 1
    nos = nnx * nny
    lx = GEOMETRY["lx"]
    ly = GEOMETRY["ly"]
    dx = np.arange(0, lx + lx / nelx, lx / nelx)
    dy = np.arange(0, ly + ly / nely, ly / nely)
    xv, yv = np.meshgrid(dx, dy)
    conec = np.zeros((nel, 5))
    coord = np.zeros((nos, 4))
    for i in range(1, nny + 1, 1):
        coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 0] = np.arange(
            (i - 1) * nnx + 1, i * nnx + 1, 1
        ).tolist()
        coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 1] = xv[i - 1, :]
        coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 2] = yv[i - 1, :]
    for i in range(1, nely + 1, 1):
        conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 0] = np.arange(
            (i - 1) * nelx + 1, i * nelx + 1, 1
        ).tolist()
        conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 1] = np.arange(
            (i - 1) * nnx + 1, i * nnx, 1
        ).tolist()
        conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 2] = np.arange(
            (i - 1) * nnx + 2, i * nnx + 1, 1
        ).tolist()
        conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 3] = np.arange(
            i * nnx + 2, (i + 1) * nnx + 1, 1
        ).tolist()
        conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 4] = np.arange(
            i * nnx + 1, (i + 1) * nnx, 1
        ).tolist()
    return conec, coord
