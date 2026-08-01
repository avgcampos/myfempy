# distutils: language=c
# cython: language_level=3
# distutils: extra_compile_args=-fopenmp
# distutils: extra_link_args=-fopenmp
from cython cimport boundscheck, wraparound
import numpy as np
cimport numpy as np    # CORREÇÃO 1: Permite o uso de tipos C como np.int32_t e np.float64_t
cimport numpy as cnp
from cpython.mem cimport PyMem_Malloc, PyMem_Free

# OBRIGATÓRIO: Inicializa a API C do NumPy
cnp.import_array()

ctypedef np.int32_t INT32_t
ctypedef np.float64_t FLT64_t


@boundscheck(False)
@wraparound(False)         
def getVectorization(INT32_t [::1] ith, INT32_t [::1] jth, FLT64_t [::1] val, INT32_t [::1] loc, FLT64_t [:, ::1] matrix, INT32_t ee, INT32_t elemdof):
    cdef Py_ssize_t LOOP_MAX = elemdof  
    cdef FLT64_t VAL
    cdef INT32_t KI, KJ
    cdef Py_ssize_t ii, jj, idx  # CORREÇÃO 3: Calcula o índice antecipadamente com Py_ssize_t
    cdef FLT64_t [:, ::1] mat_view = matrix
    cdef INT32_t [::1] loc_view = loc
    cdef INT32_t [::1] ith_view = ith
    cdef INT32_t [::1] jth_view = jth 
    cdef FLT64_t [::1] val_view = val
    
    for ii in range(LOOP_MAX):
        KI = loc_view[ii]
        for jj in range(LOOP_MAX):
            KJ = loc_view[jj]
            VAL = mat_view[ii, jj]
            idx = (LOOP_MAX * LOOP_MAX) * ee + LOOP_MAX * ii + jj
            ith_view[idx] = KI
            jth_view[idx] = KJ
            val_view[idx] = VAL
    return ith, jth, val


@boundscheck(False)
@wraparound(False) 
def getLoadAssembler(FLT64_t [:, ::1] loadaply_obj, INT32_t nodetot, INT32_t nodedof):
    cdef Py_ssize_t n_rows = loadaply_obj.shape[0] if loadaply_obj.size > 0 else 0
    cdef INT32_t total_dof = nodedof * nodetot

    if n_rows == 0:
        return np.zeros((total_dof, 1), dtype=np.float64)

    cdef FLT64_t [:, ::1] loadaply = np.ascontiguousarray(loadaply_obj, dtype=np.float64)

    cdef Py_ssize_t i
    cdef INT32_t step_val, max_step = 0, steps = 1
    cdef INT32_t node, dof, fstep, gdlload
    cdef FLT64_t val

    for i in range(n_rows):
        step_val = <INT32_t>(loadaply[i, 3])
        if step_val > max_step:
            max_step = step_val

    if max_step > 0:
        steps = max_step

    cdef cnp.ndarray[FLT64_t, ndim=2] forcevec_obj = np.zeros((total_dof, steps), dtype=np.float64)
    cdef FLT64_t [:, ::1] forcevec = forcevec_obj

    for i in range(n_rows):
        node = <INT32_t>(loadaply[i, 0]) - 1
        dof = <INT32_t>(loadaply[i, 1]) - 1
        val = loadaply[i, 2]
        fstep = <INT32_t>(loadaply[i, 3]) - 1

        gdlload = nodedof * node + dof
        forcevec[gdlload, fstep] += val

    return forcevec_obj


@boundscheck(False)
@wraparound(False) 
def getConstrains(FLT64_t [:, ::1] constrains_obj, INT32_t nodetot, INT32_t nodedof):
    cdef Py_ssize_t n_rows = constrains_obj.shape[0] if constrains_obj.size > 0 else 0
    cdef Py_ssize_t total_dof = <Py_ssize_t>(nodedof * nodetot)

    cdef char* is_fixed = <char*> PyMem_Malloc(total_dof * sizeof(char))
    cdef char* is_const = <char*> PyMem_Malloc(total_dof * sizeof(char))

    if not is_fixed or not is_const:
        if is_fixed: PyMem_Free(is_fixed)
        if is_const: PyMem_Free(is_const)
        raise MemoryError("Falha ao alocar memória nos seletores de graus de liberdade.")

    # CORREÇÃO 2: Usa Py_ssize_t para indexar ponteiros C nativos (char*)
    cdef Py_ssize_t j, gdl, d
    for j in range(total_dof):
        is_fixed[j] = 0
        is_const[j] = 0

    cdef FLT64_t [:, ::1] constrains = np.ascontiguousarray(constrains_obj, dtype=np.float64)
    cdef Py_ssize_t i
    cdef INT32_t node, dof_spec
    cdef FLT64_t val

    for i in range(n_rows):
        node = <INT32_t>(constrains[i, 0]) - 1
        dof_spec = <INT32_t>(constrains[i, 1])
        val = constrains[i, 2]

        if dof_spec == 0:
            for d in range(nodedof):
                gdl = <Py_ssize_t>(node * nodedof + d)
                if val == 0.0:
                    is_fixed[gdl] = 1
                else:
                    is_const[gdl] = 1
        else:
            gdl = <Py_ssize_t>(node * nodedof + (dof_spec - 1))
            if val == 0.0:
                is_fixed[gdl] = 1
            else:
                is_const[gdl] = 1

    cdef INT32_t n_fixed = 0, n_const = 0, n_free = 0
    for j in range(total_dof):
        if is_fixed[j] == 1:
            n_fixed += 1
        elif is_const[j] == 1:
            n_const += 1
        else:
            n_free += 1

    cdef cnp.ndarray[INT32_t, ndim=1] fixedof_obj = np.empty(n_fixed, dtype=np.int32)
    cdef cnp.ndarray[INT32_t, ndim=1] constdof_obj = np.empty(n_const, dtype=np.int32)
    cdef cnp.ndarray[INT32_t, ndim=1] freedof_obj = np.empty(n_free, dtype=np.int32)

    cdef INT32_t [::1] fixedof = fixedof_obj
    cdef INT32_t [::1] constdof = constdof_obj
    cdef INT32_t [::1] freedof = freedof_obj

    cdef INT32_t idx_fixed = 0, idx_const = 0, idx_free = 0

    for j in range(total_dof):
        if is_fixed[j] == 1:
            fixedof[idx_fixed] = <INT32_t>j
            idx_fixed += 1
        elif is_const[j] == 1:
            constdof[idx_const] = <INT32_t>j
            idx_const += 1
        else:
            freedof[idx_free] = <INT32_t>j
            idx_free += 1

    PyMem_Free(is_fixed)
    PyMem_Free(is_const)

    return freedof_obj, fixedof_obj, constdof_obj


@boundscheck(False)
@wraparound(False) 
def getDirichletNH(FLT64_t [:, ::1] constrains_obj, INT32_t nodetot, INT32_t nodedof):
    cdef Py_ssize_t n_rows = constrains_obj.shape[0] if constrains_obj.size > 0 else 0
    cdef INT32_t total_dof = nodedof * nodetot

    if n_rows == 0:
        return np.zeros((total_dof, 1), dtype=np.float64)

    cdef FLT64_t [:, ::1] constrains = np.ascontiguousarray(constrains_obj, dtype=np.float64)

    cdef Py_ssize_t i
    cdef INT32_t cstep_raw, max_step = 0, steps = 1
    cdef INT32_t node, dof, cstep, gdlload
    cdef FLT64_t val

    for i in range(n_rows):
        cstep_raw = <INT32_t>(constrains[i, 3])
        if cstep_raw != 0 and cstep_raw > max_step:
            max_step = cstep_raw

    if max_step > 0:
        steps = max_step

    cdef cnp.ndarray[FLT64_t, ndim=2] Uc_obj = np.zeros((total_dof, steps), dtype=np.float64)
    cdef FLT64_t [:, ::1] Uc = Uc_obj

    if max_step == 0:
        return Uc_obj

    for i in range(n_rows):
        cstep_raw = <INT32_t>(constrains[i, 3])
        if cstep_raw == 0:
            continue

        node = <INT32_t>(constrains[i, 0]) - 1
        dof = <INT32_t>(constrains[i, 1]) - 1
        val = constrains[i, 2]
        cstep = cstep_raw - 1

        gdlload = nodedof * node + dof
        Uc[gdlload, cstep] = val

    return Uc_obj