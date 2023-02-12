#!/usr/bin/env python
"""
Physics Setting
"""
import numpy as np


def gen_force(forcelist):
    """_summary_

    Arguments:
        forcelist -- _description_

    Returns:
        _description_
    """
    nforc = len(forcelist)
    flist = np.zeros((1, 9))
    for fl in range(nforc):
        fap = forcelist[fl]
        for fs in range(len(fap["VAL"])):
            if "LOC" in fap.keys():
                if "TAG" in fap.keys():
                    linearray = np.array(
                        [
                            fap["DEF"],
                            fap["DOF"],
                            fap["VAL"][fs],
                            fap["DIR"],
                            fap["LOC"]["x"],
                            fap["LOC"]["y"],
                            fap["LOC"]["z"],
                            int(fs + 1),
                            fap["TAG"],
                        ]
                    )
                else:
                    linearray = np.array(
                        [
                            fap["DEF"],
                            fap["DOF"],
                            fap["VAL"][fs],
                            fap["DIR"],
                            fap["LOC"]["x"],
                            fap["LOC"]["y"],
                            fap["LOC"]["z"],
                            int(fs + 1),
                            0,
                        ]
                    )
            else:
                if "TAG" in fap.keys():
                    linearray = np.array(
                        [
                            fap["DEF"],
                            fap["DOF"],
                            fap["VAL"][fs],
                            fap["DIR"],
                            0.0,
                            0.0,
                            0.0,
                            int(fs + 1),
                            fap["TAG"],
                        ]
                    )
                else:
                    linearray = np.array(
                        [
                            fap["DEF"],
                            fap["DOF"],
                            fap["VAL"][fs],
                            fap["DIR"],
                            0.0,
                            0.0,
                            0.0,
                            int(fs + 1),
                            0,
                        ]
                    )
            flist = np.append(flist, [linearray], axis=0)
    flist = flist[1::][::]
    return flist


def gen_bound(boundcondlist):
    """_summary_

    Arguments:
        boundcondlist -- _description_

    Returns:
        _description_
    """
    nbound = len(boundcondlist)
    blist = np.zeros((1, 7))
    for bl in range(nbound):
        bap = boundcondlist[bl]
        if "LOC" in bap.keys():
            if "TAG" in bap.keys():
                linearray = np.array(
                    [
                        bap["DEF"],
                        bap["DOF"],
                        bap["DIR"],
                        bap["LOC"]["x"],
                        bap["LOC"]["y"],
                        bap["LOC"]["z"],
                        bap["TAG"],
                    ]
                )
            else:
                linearray = np.array(
                    [
                        bap["DEF"],
                        bap["DOF"],
                        bap["DIR"],
                        bap["LOC"]["x"],
                        bap["LOC"]["y"],
                        bap["LOC"]["z"],
                        0,
                    ]
                )
        else:
            if "TAG" in bap.keys():
                linearray = np.array(
                    [bap["DEF"], bap["DOF"], bap["DIR"], 0.0, 0.0, 0.0, bap["TAG"]]
                )
            else:
                linearray = np.array(
                    [bap["DEF"], bap["DOF"], bap["DIR"], 0.0, 0.0, 0.0, 0]
                )
        blist = np.append(blist, [linearray], axis=0)
    blist = blist[1::][::]
    return blist
