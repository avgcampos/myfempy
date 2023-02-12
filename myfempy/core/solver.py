#!/usr/bin/env python
from myfempy.utils.utils import print_console
from myfempy.core.solverset import get_constrains_dofs, step_setting, get_solve
from myfempy.core.assembler import Assembler
import time
import sys

__doc__ = """
Solver Manager
"""


class Solver:
    """_summary_

    Returns:
        _description_
    """

   
    @staticmethod
    def get_static_solve(solverset, modelinfo):
        """_summary_

        Arguments:
            solverset -- _description_
            modelinfo -- _description_

        Returns:
            _description_
        """
        # print_console("solver")
        solve = get_solve(solverset["SOLVER"])
        start = time.time()
        KG = Assembler.assembler(modelinfo, key="stiffness")
        end = time.time()
        kg_mem_size = sys.getsizeof(KG.toarray())/1e6
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
    def get_modal_solve(solverset, modelinfo):
        """_summary_

        Arguments:
            solverset -- _description_
            modelinfo -- _description_

        Returns:
            _description_
        """
        # print_console("solver")
        solve = get_solve(solverset["SOLVER"])
        start = time.time()
        KG = Assembler.assembler(modelinfo, key="stiffness")
        MG = Assembler.assembler(modelinfo, key="mass")
        end = time.time()
        kg_mem_size = sys.getsizeof(KG.toarray())/1e6
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
