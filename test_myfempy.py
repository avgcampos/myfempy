import numpy as np
import pytest
from myfempy import SteadyStateLinear, newAnalysis


@pytest.fixture
def fem_setup():
    """Fixture para configurar a geometria, propriedades e condições de contorno."""
    fea = newAnalysis(SteadyStateLinear)

    mat1 = {"NAME": "material1", "EXX": 1000, "VXY": 0.3}
    mat2 = {"NAME": "material2", "EXX": 2000, "VXY": 0.3}
    geo = {"NAME": "geo", "THICKN": 0.1}

    nodes = [
        [1, 0.00, 0.00, 0.00],
        [2, 0.50, 0.00, 0.00],
        [3, 1.00, 0.00, 0.00],
        [4, 0.00, 0.50, 0.00],
        [5, 0.50, 0.50, 0.00],
        [6, 1.00, 0.50, 0.00],
        [7, 0.00, 1.00, 0.00],
        [8, 0.50, 1.00, 0.00],
        [9, 1.00, 1.00, 0.00],
    ]

    conec = [
        [1, 1, 1, 1, 2, 5, 4],
        [2, 1, 1, 2, 3, 6, 5],
        [3, 1, 1, 4, 5, 8, 7],
        [4, 2, 1, 5, 6, 9, 8],
    ]

    modeldata = {
        "MESH": {"TYPE": "manual", "COORD": nodes, "INCI": conec},
        "ELEMENT": {"TYPE": "structplane", "SHAPE": "quad4"},
        "MATERIAL": {
            "MAT": "planestress",
            "TYPE": "isotropic",
            "PROPMAT": [mat1, mat2],
        },
        "GEOMETRY": {"GEO": "thickness", "PROPGEO": [geo]},
    }

    bcfix = {
        "TYPE": "fixed",
        "DOF": "full",
        "DIR": "node",
        "MESHNODE": [1, 4, 7],
    }

    force = {
        "TYPE": "forcenode",
        "DOF": "fy",
        "DIR": "node",
        "MESHNODE": [6],
        "VAL": [-100],
    }

    physicdata = {
        "PHYSIC": {
            "DOMAIN": "structural",
            "LOAD": [force],
            "BOUNDCOND": [bcfix],
        }
    }

    fea.Model(modeldata)
    fea.Physic(physicdata)

    return fea


def test_fea_solver_execution(fem_setup):
    """Testa se a etapa de solução roda e retorna uma estrutura de deslocamentos válida."""
    fea = fem_setup

    solverset = {"STEPSET": {"type": "table", "start": 0, "end": 1, "step": 1}}
    solverdata = fea.Solve(solverset)

    # Verifica se os dados do solver não são nulos e contêm a chave da solução
    assert solverdata is not None
    assert "solution" in solverdata
    assert "U" in solverdata["solution"]

    # Verifica se os deslocamentos foram calculados
    u_vec = np.array(solverdata["solution"]["U"])
    assert u_vec.size > 0
    assert not np.isnan(u_vec).any(), "Existem valores NaN na solução de deslocamentos"


def test_fixed_boundary_conditions(fem_setup):
    """Garante que os nós engastados (1, 4, 7) permaneçam com deslocamento zero."""
    fea = fem_setup

    solverset = {"STEPSET": {"type": "table", "start": 0, "end": 1, "step": 1}}
    solverdata = fea.Solve(solverset)
    u_vec = np.array(solverdata["solution"]["U"])

    # Assumindo que cada nó possui 2 DOFs (ux, uy) organizados sequencialmente:
    fixed_nodes = [1, 4, 7]
    for node_id in fixed_nodes:
        idx_x = (node_id - 1) * 2
        idx_y = idx_x + 1

        assert u_vec[idx_x] == pytest.approx(0.0, abs=1e-8)
        assert u_vec[idx_y] == pytest.approx(0.0, abs=1e-8)


def test_postprocessing_outputs(fem_setup):
    """Valida se o pós-processamento gera o dicionário com as tensões esperadas."""
    fea = fem_setup

    solverset = {"STEPSET": {"type": "table", "start": 0, "end": 1, "step": 1}}
    solverdata = fea.Solve(solverset)

    postprocset = {
        "SOLVERDATA": solverdata,
        "COMPUTER": {"structural": {"displ": True, "stress": True}},
        "PLOTSET": {"show": False, "filename": "test_output", "savepng": False},
        "REPORT": {
            "log": False,
            "get": {"nelem": True, "nnode": True, "tabmat": True},
        },
    }

    postprocdata = fea.PostProcess(postprocset)

    # Asserções do pós-processamento
    assert postprocdata is not None
    assert "STRESS_XX" in postprocdata

    stress_xx = np.array(postprocdata["STRESS_XX"])
    assert stress_xx.size > 0
    assert not np.isnan(stress_xx).any(), "Tensões obtidas contêm valores NaN"