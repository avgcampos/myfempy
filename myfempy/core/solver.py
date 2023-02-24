#!/usr/bin/env python
import sys
import time

from myfempy.core.assembler import Assembler
from myfempy.core.solverset import get_constrains_dofs, get_solve, step_setting
from myfempy.utils.utils import print_console

__doc__ = """
Solver Manager
"""


class Solver:
    """class solver
    
    SLD     -- scipy sparse linear solver 
    SLI     -- scipy sparse biconjugate gradient stabilized iteration solver 
    SLIPRE  -- scipy generalized minimal residual iteration solver
    EIG     -- scipy eigenvalues and eigenvectors solver
    FRF     -- scipy sparse linear steps(frequency) solver 
    
    """

    @staticmethod
    def get_static_solve(solverset: dict, modelinfo: dict):
        """get a static solution

        Arguments:
            solverset:dict  -- solver setting
            modelinfo:dict  -- F.E. model dict with full information needed

        Returns:
            solution:dict   -- solution
        """
        
        # print_console("solver")
        solve = get_solve(solverset["SOLVER"])
        start = time.time()
        KG = Assembler.assembler(modelinfo, key="stiffness")
        end = time.time()
        kg_mem_size = sys.getsizeof(KG.toarray()) / 1e6
        assembly_time = end - start
        print(" ")
        # print("STIFFNESS SIZE: ", kg_mem_size, " MB")
        # print("ASSEMBLY FULL TIME SPEND ", assembly_time, " SEC")
        F, KG = Assembler.loads(modelinfo, KG)
        freedof, fixedof = get_constrains_dofs(modelinfo)
        solverset["nsteps"] = step_setting(solverset["STEPSET"])
        solverset["coord"] = modelinfo["coord"]
        solverset["nodedof"] = modelinfo["nodedof"][0]
        fulldofs = (modelinfo["nodedof"][0]) * len(modelinfo["coord"])
        start = time.time()
        U = solve(fulldofs, KG, F, freedof, solverset)
        end = time.time()
        time_spend = end - start
        # print("\nSOLVE FULL TIME SPEND ", time_spend, " SEC")
        solvestatus = {
            "timeasb": assembly_time,
            "timesim": time_spend,
            "kgsize": kg_mem_size,
        }
        solution = {"U": U, "solvestatus": solvestatus}
        print_console("thank")
        return solution

    @staticmethod
    def get_modal_solve(solverset: dict, modelinfo: dict):
        """get a modal solution

        Arguments:
            solverset:dict  -- solver setting
            modelinfo:dict  -- F.E. model dict with full information needed

        Returns:
            solution:dict   -- solution
        """
        
        # print_console("solver")
        solve = get_solve(solverset["SOLVER"])
        start = time.time()
        KG = Assembler.assembler(modelinfo, key="stiffness")
        MG = Assembler.assembler(modelinfo, key="mass")
        end = time.time()
        kg_mem_size = sys.getsizeof(KG.toarray()) / 1e6
        assembly_time = end - start
        # print(" ")
        # print("STIFFNESS SIZE: ", kg_mem_size, " MB")
        # print("ASSEMBLY FULL TIME SPEND ", assembly_time, " SEC")
        F, KG = Assembler.loads(modelinfo, KG)
        freedof, fixedof = get_constrains_dofs(modelinfo)
        solverset["start"] = solverset["STEPSET"]["start"]
        solverset["end"] = solverset["STEPSET"]["end"]
        solverset["nsteps"] = step_setting(solverset["STEPSET"])
        solverset["coord"] = modelinfo["coord"]
        fulldofs = (modelinfo["nodedof"][0]) * len(modelinfo["coord"])
        start = time.time()
        U, w_range = solve(fulldofs, KG, MG, F, freedof, solverset)
        end = time.time()
        time_spend = end - start
        # print("\nSOLVE FULL TIME SPEND ", time_spend, " SEC")
        solvestatus = {
            "timeasb": assembly_time,
            "timesim": time_spend,
            "kgsize": kg_mem_size,
        }
        solution = {
            "U": U,
            "FREQ": w_range,
            "solvestatus": solvestatus,
        }
        # print_console("thank")
        return solution


if __name__ == "__main__":
    import doctest

    doctest.testmod()
