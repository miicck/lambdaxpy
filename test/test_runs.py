import os

LAM_X_PY = os.path.join(os.path.dirname(os.path.dirname(__file__)), "lambda.x.py")
assert os.path.isfile(LAM_X_PY)


def run_directory(directory: str):
    assert os.path.isdir(directory)

    # Check input is present, output is not
    inp_file = os.path.join(directory, "lambda.in")
    out_file = os.path.join(directory, "lambda.out")
    assert os.path.isfile(inp_file)
    assert not os.path.isfile(out_file)

    try:
        # Run lambda.x.py
        os.system(f"cd {directory} && python {LAM_X_PY} lambda.in > lambda.out")

        # Check output
        with open(out_file) as f:
            for line in f:
                print(line.strip())

    finally:

        # Cleanup
        if os.path.isfile(out_file):
            os.remove(out_file)


def test_mg2irh6():
    run_directory("data/Mg2IrH6")
