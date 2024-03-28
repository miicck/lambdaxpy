import os
import numpy as np


def compare_x_xpy(x_file: str, xpy_file: str):
    data = None
    deguass = []
    nef = []
    with open(x_file) as f:
        for line in f:
            split = line.split()
            if split[1] == "=":
                deguass.append(float(split[-1]))
                nef.append(float(split[-4]))

            if split == ["lambda", "omega_log", "T_c"]:
                data = []
                continue

            if data is None:
                continue

            data.append([float(x) for x in split])

    x_data = np.array(data).T

    data = []
    with open(xpy_file) as f:
        for i, line in enumerate(f):
            split = line.split()
            if i > 0:
                data.append([float(x) for x in split])

    xpy_data = np.array(data).T

    # Degauss should be exactly the same
    assert np.allclose(deguass, xpy_data[0])

    # N(Ef) should be exactly the same
    assert np.allclose(nef, xpy_data[1])

    # Lambda should be exactly the same, both codes use the analytic expression
    assert np.allclose(x_data[0], xpy_data[2])

    # Omega log should be close, lambda.x.py uses analytic, lambda.x uses smearing
    assert np.allclose(x_data[1], xpy_data[3], rtol=0.01, atol=1e-5)

    # Same goes for Tc
    assert np.allclose(x_data[2], xpy_data[4], rtol=0.01, atol=1e-5)


def run_directory(directory: str, manual_smearing: float = None):
    assert os.path.isdir(directory)

    args = ""
    if manual_smearing is not None:
        args += f"--a2f_smearing {manual_smearing}"

    try:
        # Run lambda.x.py
        os.system(f"cd {directory} && lambdaxpy {args} lambda.in > lambda.out")
        os.system(f"cd {directory} && lambda.x < lambda.in > lambda.x.out")

        # Check output
        compare_x_xpy(os.path.join(directory, "lambda.x.out"),
                      os.path.join(directory, "lambda.out"))

    finally:

        # Cleanup
        for f in ["lambda.out", "lambda.x.out", "alpha2F.dat", "lambda.dat"]:
            f = os.path.join(directory, f)
            if os.path.isfile(f):
                os.remove(f)


def test_mg2irh6():
    run_directory(os.path.join(os.path.dirname(__file__), "data", "Mg2IrH6"))


def test_all_dirs():
    root = os.path.join(os.path.dirname(__file__), "data")
    for d in os.listdir(root):
        d = os.path.join(root, d)
        if os.path.isdir(d):
            run_directory(d)
            run_directory(d, manual_smearing=1.0)
