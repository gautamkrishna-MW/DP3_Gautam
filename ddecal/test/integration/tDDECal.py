# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import pytest
import os
import shutil
import uuid
from subprocess import check_call, check_output, CalledProcessError, STDOUT
import numpy as np

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf
from testconfig import TAQLEXE
from utils import assert_taql, untar_ms, get_taql_result

"""
Script can be invoked in two ways:
- as standalone from the build/ddecal/test/integration directory,
  using `pytest source/tDDECal.py` (extended with pytest options of your choice)
- using ctest, see DP3/ddecal/test/integration/CMakeLists.txt
"""

MSINTGZ = "tDDECal.in_MS.tgz"
MSIN = "tDDECal.MS"
CWD = os.getcwd()


@pytest.fixture(autouse=True)
def source_env():
    os.chdir(CWD)
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    os.chdir(tmpdir)

    untar_ms(f"{tcf.RESOURCEDIR}/{MSINTGZ}")
    check_call([tcf.MAKESOURCEDBEXE, f"in={MSIN}/sky.txt", f"out={MSIN}/sky"])

    # Tests are executed here
    yield

    # Post-test: clean up
    os.chdir(CWD)
    shutil.rmtree(tmpdir)


@pytest.fixture()
def copy_data_to_model_data():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=MODEL_DATA",
            "steps=[]",
        ]
    )


@pytest.fixture()
def create_corrupted_visibilities():
    taqlcommand = f"update {MSIN} set WEIGHT_SPECTRUM=1, FLAG=False"
    check_output([TAQLEXE, "-noph", taqlcommand])

    # Use ddecal to create template h5parm
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrumentcorrupted.h5",
            "ddecal.mode=complexgain",
        ]
    )

    # Modify h5 file for multiple solution intervals
    import h5py  # Don't import h5py when pytest is only collecting tests.

    h5file = h5py.File("instrumentcorrupted.h5", "r+")
    sol = h5file["sol000/amplitude000/val"]
    sol[:4, ..., 0, :] = np.sqrt(5)
    sol[4:, ..., 0, :] = np.sqrt(5 + 2)
    sol[:4, ..., 1, :] = np.sqrt(9)
    sol[4:, ..., 1, :] = np.sqrt(9 + 2)
    sol[:4, ..., 2, :] = np.sqrt(13)
    sol[4:, ..., 2, :] = np.sqrt(13 + 2)
    h5file.close()

    # Predict corrupted visibilities into DATA column
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=DATA",
            "steps=[h5parmpredict]",
            f"h5parmpredict.sourcedb={MSIN}/sky",
            "h5parmpredict.applycal.parmdb=instrumentcorrupted.h5",
            "h5parmpredict.applycal.correction=amplitude000",
        ]
    )


@pytest.mark.parametrize(
    "caltype",
    ["complexgain", "scalarcomplexgain", "amplitudeonly", "scalaramplitude"],
)
@pytest.mark.parametrize("solint", [0, 1, 2, 4])
@pytest.mark.parametrize("nchan", [1, 2, 5])
def test(
    create_corrupted_visibilities,
    copy_data_to_model_data,
    caltype,
    solint,
    nchan,
):
    # Subtract corrupted visibilities using multiple predict steps
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            f"ddecal.solint={solint}",
            f"ddecal.nchan={nchan}",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
        ]
    )

    # Calibrate on the original sources, caltype=$caltype
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3]",
            f"predict1.sourcedb={MSIN}/sky",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center,dec_off]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={MSIN}/sky",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[radec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={MSIN}/sky",
            "predict3.applycal.parmdb=instrument.h5",
            "predict3.sources=[ra_off]",
            "predict3.operation=subtract",
            "predict3.applycal.correction=amplitude000",
        ]
    )

    if solint == 0:
        tolerance = 0.15
    else:
        tolerance = 0.015
    taqlcommand_run = f"select norm_residual/norm_data FROM (select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*DATA))) as norm_data, sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*SUBTRACTED_DATA))) as norm_residual from {MSIN})"
    check_output([TAQLEXE, "-noph", taqlcommand_run])
    taql_command = f"select FROM (select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*DATA))) as norm_data, sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*SUBTRACTED_DATA))) as norm_residual from {MSIN}) where norm_residual/norm_data > {tolerance} or isinf(norm_residual/norm_data) or isnan(norm_residual/norm_data)"
    assert_taql(taql_command)


@pytest.mark.parametrize(
    "solutions_per_direction", [[1, 1, 1], [1, 3, 6], [3, 3, 3], [6, 6, 6]]
)
@pytest.mark.parametrize(
    "caltype",
    ["scalaramplitude", "scalar"],
)
def test_calibration_with_dd_intervals(
    create_corrupted_visibilities,
    copy_data_to_model_data,
    solutions_per_direction,
    caltype,
):
    # This test checks that the calibration with different solution intervals per direction gives the same result as the corruption applied in the fixture "create_corrupted_visibilities".

    # Run calibration on corrupted dataset and store solutions in instrument.h5
    sol_int = 6
    n_timeslots_in_ms = 6

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            f"ddecal.solint={sol_int}",
            "ddecal.nchan=1",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
            "ddecal.solveralgorithm=directioniterative",
            f"ddecal.solutions_per_direction={solutions_per_direction}",
        ]
    )

    # Verify that he values in instrument.h5 are close to the values in instrumentcorrupted.h5
    import h5py  # Don't import h5py when pytest is only collecting tests.

    ddecal_solutions = h5py.File("instrument.h5", "r")
    reference_solutions = h5py.File("instrumentcorrupted.h5", "r")

    assert (
        ddecal_solutions["sol000/amplitude000/val"].attrs["AXES"]
        == b"time,freq,ant,dir"
    )

    for direction_index in range(
        0, len(solutions_per_direction)
    ):  # loop over number of directions
        for solint_index in range(
            0, int(np.ceil(n_timeslots_in_ms / sol_int))
        ):  # loop over number of solution intervals
            for interval_in_direction_index in range(
                0, solutions_per_direction[direction_index]
            ):  # loop over dd-intervals within one solution interval
                values_ddecal = ddecal_solutions["sol000/amplitude000/val"][
                    solint_index * np.max(solutions_per_direction)
                    + interval_in_direction_index,
                    :,
                    :,
                    direction_index,
                ]
                corresponding_index = solint_index * np.max(
                    solutions_per_direction
                ) + int(sol_int / solutions_per_direction[direction_index])

                if (
                    corresponding_index
                    >= reference_solutions["sol000/amplitude000/val"].shape[0]
                ):
                    corresponding_index = (
                        reference_solutions["sol000/amplitude000/val"].shape[0]
                        - 1
                    )

                values_reference = reference_solutions[
                    "sol000/amplitude000/val"
                ][corresponding_index, :, :, direction_index, 0]

                assert (abs(values_ddecal - values_reference) < 1.2).all()


@pytest.mark.xfail(reason="Will be fixed in AST-924")
@pytest.mark.parametrize(
    "solutions_per_direction",
    [[1, 1, 1], [1, 2, 4], [2, 2, 2], [4, 4, 4]],
)
@pytest.mark.parametrize(
    "caltype",
    ["scalaramplitude", "scalar", "diagonal", "diagonalamplitude"],
)
def test_bug_ast_924(
    create_corrupted_visibilities,
    copy_data_to_model_data,
    solutions_per_direction,
    caltype,
):
    # This is the same test as "test_calibration_with_dd_intervals" , but with a different solution interval:
    # In this case sol_int is not an integer divisior of n_timeslots_in_ms

    sol_int = 4
    n_timeslots_in_ms = 6
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            f"ddecal.solint={sol_int}",
            "ddecal.nchan=1",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
            "ddecal.solveralgorithm=directioniterative",
            f"ddecal.solutions_per_direction={solutions_per_direction}",
        ]
    )

    # Verify that he values in instrument.h5 are close to the values in instrumentcorrupted.h5
    import h5py  # Don't import h5py when pytest is only collecting tests.

    ddecal_solutions = h5py.File("instrument.h5", "r")
    reference_solutions = h5py.File("instrumentcorrupted.h5", "r")

    for direction_index in range(
        0, len(solutions_per_direction)
    ):  # loop over number of directions
        for solint_index in range(
            0, int(np.ceil(n_timeslots_in_ms / sol_int))
        ):  # loop over number of solution intervals
            for interval_in_direction_index in range(
                0, solutions_per_direction[direction_index]
            ):  # loop over dd-intervals within one solution interval
                values_ddecal = ddecal_solutions["sol000/amplitude000/val"][
                    solint_index * np.max(solutions_per_direction)
                    + interval_in_direction_index,
                    :,
                    :,
                    direction_index,
                ]
                corresponding_index = solint_index * np.max(
                    solutions_per_direction
                ) + int(sol_int / solutions_per_direction[direction_index])

                if (
                    corresponding_index
                    >= reference_solutions["sol000/amplitude000/val"].shape[0]
                ):
                    corresponding_index = (
                        reference_solutions["sol000/amplitude000/val"].shape[0]
                        - 1
                    )

                values_reference = reference_solutions[
                    "sol000/amplitude000/val"
                ][corresponding_index, :, :, direction_index, 0]

                assert (abs(values_ddecal - values_reference) < 1.2).all()


@pytest.mark.xfail(reason="Will be implemented in AST-920")
@pytest.mark.parametrize(
    "solutions_per_direction", [[1, 1, 1], [1, 3, 6], [3, 3, 3], [6, 6, 6]]
)
@pytest.mark.parametrize(
    "caltype",
    [
        "amplitudeonly",
        "scalaramplitude",
        "scalar",
        "diagonal",
        "diagonalamplitude",
    ],
)
def test_subtract_with_dd_intervals(
    create_corrupted_visibilities,
    copy_data_to_model_data,
    solutions_per_direction,
    caltype,
):
    # This test checks that the subtraction operation works correctly
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DURING_DDECAL",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            f"ddecal.solint=6",
            "ddecal.nchan=2",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode={caltype}",
            "ddecal.solveralgorithm=directioniterative",
            f"ddecal.solutions_per_direction={solutions_per_direction}",
            "ddecal.subtract=true",
        ]
    )

    common_predict_command = [
        tcf.DP3EXE,
        "checkparset=1",
        "numthreads=1",
        f"msin={MSIN}",
        "msout=.",
        "msout.datacolumn=SUBTRACTED_DURING_PREDICT",
        "steps=[predict1,predict2,predict3]",
        f"predict1.sourcedb={MSIN}/sky",
        "predict1.applycal.parmdb=instrument.h5",
        "predict1.sources=[center,dec_off]",
        "predict1.operation=subtract",
        f"predict2.sourcedb={MSIN}/sky",
        "predict2.applycal.parmdb=instrument.h5",
        "predict2.sources=[radec_off]",
        "predict2.operation=subtract",
        f"predict3.sourcedb={MSIN}/sky",
        "predict3.applycal.parmdb=instrument.h5",
        "predict3.sources=[ra_off]",
        "predict3.operation=subtract",
    ]
    corrections_amplitude_only = [
        "predict1.applycal.correction=amplitude000",
        "predict2.applycal.correction=amplitude000",
        "predict3.applycal.correction=amplitude000",
    ]

    corrections_amplitude_and_phase = [
        "predict1.applycal.steps=[amplitude,phase]",
        "predict1.applycal.amplitude.correction=amplitude000",
        "predict1.applycal.phase.correction=phase000",
        "predict2.applycal.steps=[amplitude,phase]",
        "predict2.applycal.amplitude.correction=amplitude000",
        "predict2.applycal.phase.correction=phase000",
        "predict3.applycal.steps=[amplitude,phase]",
        "predict3.applycal.amplitude.correction=amplitude000",
        "predict3.applycal.phase.correction=phase000",
    ]

    if caltype == "scalar" or caltype == "diagonal":
        check_call(common_predict_command + corrections_amplitude_and_phase)

    else:
        check_call(common_predict_command + corrections_amplitude_only)

    # Quantify the difference between the subtraction in the DDECal step and in the Predict step
    predict_residual = float(
        check_output(
            [
                tcf.TAQLEXE,
                "-nopr",
                "-noph",
                f"select sqrt(abs(gsumsqr(WEIGHT_SPECTRUM*(SUBTRACTED_DURING_DDECAL- SUBTRACTED_DURING_PREDICT)))) from {MSIN}",
            ]
        )
    )

    tolerance = 2
    assert predict_residual < tolerance


@pytest.mark.parametrize("skymodel", ["sky", "sky.txt"])
def test_h5parm_predict(skymodel):
    # make calibration solutions
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/{skymodel}",
            # By omitting the direction, the directions are determined from the
            # sourcedb. This gives the same results as:
            # "ddecal.directions=[[center],[dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode=diagonal",
        ]
    )

    # subtract using multiple predict steps
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3,predict4]",
            f"predict1.sourcedb={MSIN}/{skymodel}",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={MSIN}/{skymodel}",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[dec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={MSIN}/{skymodel}",
            "predict3.applycal.parmdb=instrument.h5",
            "predict3.sources=[radec_off]",
            "predict3.operation=subtract",
            "predict3.applycal.correction=amplitude000",
            f"predict4.sourcedb={MSIN}/{skymodel}",
            "predict4.applycal.parmdb=instrument.h5",
            "predict4.sources=[ra_off]",
            "predict4.operation=subtract",
            "predict4.applycal.correction=amplitude000",
        ]
    )

    # subtract using h5parmpredict
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA_H5PARM",
            "steps=[h5parmpredict]",
            f"h5parmpredict.sourcedb={MSIN}/{skymodel}",
            "h5parmpredict.applycal.parmdb=instrument.h5",
            "h5parmpredict.operation=subtract",
            "h5parmpredict.applycal.correction=amplitude000",
        ]
    )

    # echo "Check that h5parmpredict creates the same output as multiple predict steps"
    taql_command = f"select from (select abs(gsumsqr(SUBTRACTED_DATA-SUBTRACTED_DATA_H5PARM)) as diff from {MSIN}) where diff>1.e-6"
    assert_taql(taql_command)


def test_oneapplycal_from_buffer():
    skymodel = "sky.txt"

    # Apply ddecal and oneapplycal sequentially
    check_call(
        [
            tcf.DP3EXE,
            "numthreads=1",
            f"msin={MSIN}",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/{skymodel}",
            "ddecal.h5parm=cal.h5",
            "ddecal.directions=[[center,dec_off,ra_off,radec_off]]",
            "msout=.",
            "msout.datacolumn=DATA_REF_BUF",
            "checkparset=1",
        ]
    )
    check_call(
        [
            tcf.DP3EXE,
            "numthreads=1",
            f"msin={MSIN}",
            "msin.datacolumn=DATA_REF_BUF",
            "steps=[applycal]",
            "applycal.parmdb=cal.h5",
            "applycal.steps=[phase,ampl]",
            "applycal.phase.correction=phase000",
            "applycal.ampl.correction=amplitude000",
            "msout=.",
            "msout.datacolumn=DATA_REF_BUF",
            "checkparset=1",
        ]
    )

    # DDECal and immediate application
    check_call(
        [
            tcf.DP3EXE,
            "numthreads=1",
            f"msin={MSIN}",
            "steps=[ddecal,applycal]",
            f"ddecal.sourcedb={MSIN}/{skymodel}",
            "ddecal.directions=[[center,dec_off,ra_off,radec_off]]",
            "ddecal.storebuffer=True",
            "msout=.",
            "msout.datacolumn=DATA_NEW_BUF",
            "applycal.parmdb=",
            "checkparset=1",
        ]
    )

    # echo "Check that h5parmpredict creates the same output as multiple predict steps"
    taql_command = f"select from (select abs(gsumsqr(DATA_REF_BUF-DATA_NEW_BUF)) as diff from {MSIN}) where diff>1.e-6"
    assert_taql(taql_command)


def test_pre_apply():
    # make calibration solutions
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument.h5",
            f"ddecal.mode=diagonal",
        ]
    )

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[predict1,predict2,predict3]",
            f"predict1.sourcedb={MSIN}/sky",
            "predict1.applycal.parmdb=instrument.h5",
            "predict1.sources=[center,dec_off]",
            "predict1.operation=subtract",
            "predict1.applycal.correction=amplitude000",
            f"predict2.sourcedb={MSIN}/sky",
            "predict2.applycal.parmdb=instrument.h5",
            "predict2.sources=[radec_off]",
            "predict2.operation=subtract",
            "predict2.applycal.correction=amplitude000",
            f"predict3.sourcedb={MSIN}/sky",
            "predict3.applycal.parmdb=instrument.h5",
            "predict3.sources=[ra_off]",
            "predict3.operation=subtract",
            "predict3.applycal.correction=amplitude000",
        ]
    )

    # Check that preapply runs (output is not tested)
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.applycal.parmdb=instrument.h5",
            "ddecal.applycal.steps=applyampl",
            "ddecal.applycal.applyampl.correction=amplitude000",
            "ddecal.h5parm=instrument2.h5",
            "ddecal.mode=scalarcomplexgain",
        ]
    )


def test_check_tec():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument-tec.h5",
            "ddecal.mode=tec",
        ]
    )


def test_check_tec_and_phase():
    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msin.baseline='!CS001HBA0'",
            "steps=[ddecal]",
            f"ddecal.sourcedb={MSIN}/sky",
            "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
            "ddecal.h5parm=instrument-tecandphase.h5 ddecal.mode=tecandphase",
        ]
    )


@pytest.mark.parametrize(
    "solutions_per_direction", [None, [1], [2, 3, 1], [5, 5, 5], [2, 0]]
)
def test_dd_solution_intervals(solutions_per_direction):
    base_command = [
        tcf.DP3EXE,
        "checkparset=1",
        "numthreads=1",
        f"msin={MSIN}",
        "msout=.",
        "steps=[ddecal]",
        f"ddecal.sourcedb={MSIN}/sky",
        "ddecal.solint=6",
        "ddecal.directions=[[center,dec_off],[ra_off],[radec_off]]",
        "ddecal.h5parm=instrument-tec.h5",
        "ddecal.mode=scalaramplitude",
        "ddecal.solveralgorithm=directioniterative",
    ]
    try:
        check_output(
            base_command
            if solutions_per_direction is None
            else base_command
            + [f"ddecal.solutions_per_direction={solutions_per_direction}"],
            stderr=STDOUT,
        )
    except CalledProcessError as e:
        if solutions_per_direction == [5, 5, 5]:
            if not e.output.decode().startswith(
                "\nstd exception detected: Values in ddecal"
            ):
                raise e
        elif solutions_per_direction == [2, 0]:
            if not e.output.decode().startswith(
                "\nstd exception detected: All entries in ddecal"
            ):
                raise e
        else:
            raise (e)


def test_modelnextsteps(copy_data_to_model_data):
    # Multiply MODEL_DATA by 42
    taqlcommand_run = f"update {MSIN} set MODEL_DATA=DATA*42"
    check_output([TAQLEXE, "-noph", taqlcommand_run])

    check_call(
        [
            tcf.DP3EXE,
            "checkparset=1",
            "numthreads=1",
            f"msin={MSIN}",
            "msout=.",
            "msout.datacolumn=SUBTRACTED_DATA",
            "steps=[ddecal]",
            "ddecal.modeldatacolumns=[MODEL_DATA]",
            "ddecal.modelnextsteps.MODEL_DATA=[scaledata]",
            "ddecal.h5parm=instrument-modeldata",
            "ddecal.solint=2",
            "ddecal.nchan=3",
            "scaledata.stations='*'",
            "scaledata.scalesize=False",
            "scaledata.coeffs=1",
            "ddecal.subtract=True",
        ]
    )
    taql_command = f"select from (select abs(sumsqr(SUBTRACTED_DATA)/sumsqr(DATA)) as diff from {MSIN}) where diff>1.e-6"
    assert_taql(taql_command)


def test_bda_constaints():
    import h5py  # Don't import h5py when pytest is only collecting tests.

    common = [
        f"msin={MSIN}",
        "msout.overwrite=true",
        "msin.datacolumn=DATA",
        "solve.type=ddecal",
        "solve.mode=complexgain",
        "solve.usebeammodel=True",
        "solve.beammode=array_factor",
        "solve.maxiter=150",
        "solve.directions=[[center,dec_off],[ra_off],[radec_off]]",
        "numthreads=6",
        "solve.onebeamperpatch=False",
        "solve.propagatesolutions=True",
        "solve.smoothnessconstraint=4000000.0",
        "solve.solint=37",
        "solve.nchan=10",
        "solve.solveralgorithm=hybrid",
        f"solve.sourcedb={MSIN}/sky",
        "solve.stepsize=0.2",
        "solve.tolerance=0.005",
    ]

    check_call(
        [
            tcf.DP3EXE,
            "msout=tmp_bda.ms",
            "steps=[avg,solve]",
            "avg.type=bdaaverager",
            "avg.timebase=20000",
            "avg.frequencybase=10000",
            "avg.maxinterval=37",
            "solve.h5parm=test_bda.h5parm",
        ]
        + common
    )

    check_call(
        [
            tcf.DP3EXE,
            "msout=tmp_no_bda.ms",
            "steps=[solve]",
            "solve.h5parm=test_no_bda.h5parm",
        ]
        + common
    )

    f_bda = h5py.File("test_bda.h5parm", "r")
    f_no_bda = h5py.File("test_no_bda.h5parm", "r")

    ampl_bda = f_bda["sol000/amplitude000/val"]
    ampl_no_bda = f_no_bda["sol000/amplitude000/val"]
    phase_bda = f_bda["sol000/phase000/val"]
    phase_no_bda = f_no_bda["sol000/phase000/val"]

    np.testing.assert_allclose(
        ampl_bda, ampl_no_bda, rtol=0.05, atol=0, equal_nan=True
    )
    np.testing.assert_allclose(
        phase_bda, phase_no_bda, rtol=0.3, atol=0, equal_nan=True
    )
