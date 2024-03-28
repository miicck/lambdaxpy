from typing import Iterable
import os
import numpy as np
import argparse

# Constants
RY_TO_CMM = 109736.75775046606
RY_TO_K = 157887.6633481157
THZ_TO_RY = 0.0003039657635040861


class LambdaElphInput:
    """
    Class representing the elph.inp_lambda.n files read by lambda.x.
    """

    def __init__(self, filename: str):
        assert os.path.isfile(filename), f"File does not exist: {filename}"

        self.filename = filename
        self.sigmas = []
        self.efs = []
        self.dos_at_efs = []
        self.mode_lambdas = []
        self.freqs_squared = []

        i_next_modes = None

        with open(filename) as f:
            for i, line in enumerate(f):
                line = line.strip()
                split = line.split()

                if i == 0:
                    # Read q-point, number of smearing values and number of modes
                    self.q = [float(x) for x in split[:3]]
                    self.nsig = int(split[3])
                    self.nmodes = int(split[4])
                    continue

                if line.startswith("Gaussian Broadening:"):
                    self.sigmas.append(float(split[2]))
                    i_next_modes = None
                    self.mode_lambdas.append([])
                    continue

                if line.startswith("DOS ="):
                    self.dos_at_efs.append(float(split[2]))
                    self.efs.append(float(split[-2]))
                    i_next_modes = i + 1
                    continue

                if len(self.mode_lambdas) == 0:
                    # Haven't hit lambdas yet, read frequencies
                    self.freqs_squared.extend([float(x) for x in split])
                    continue

                if i_next_modes is not None and i >= i_next_modes:
                    self.mode_lambdas[-1].append(float(line.split("=")[1].split()[0]))
                    continue

                print(i, line)

        self.mode_lambdas = np.array(self.mode_lambdas).T

        assert len(self.sigmas) == self.nsig
        assert self.mode_lambdas.shape == (self.nmodes, self.nsig)
        assert len(self.efs) == self.nsig
        assert len(self.dos_at_efs) == self.nsig
        assert len(self.freqs_squared) == self.nmodes

    @property
    def frequencies(self) -> np.ndarray:
        return np.array([f ** 0.5 if f >= 0.0 else 0.0 for f in self.freqs_squared])


class LambdaInput:
    """
    Class representing the lambda.x input file
    """

    def __init__(self, filename: str):
        assert os.path.isfile(filename), f"File does not exist: {filename}"

        n_qpoints = None
        self.q_points = []
        self.q_files = []
        self.mu_star = None

        with open(filename) as f:
            for i, line in enumerate(f):
                line = line.strip()
                split = line.split()

                if i == 0:
                    # Ignore settings, we will work these out
                    continue

                if i == 1:
                    n_qpoints = int(line)
                    continue

                if 1 < i <= n_qpoints + 1:
                    self.q_points.append([float(x) for x in split])
                    continue

                if n_qpoints + 1 < i <= 2 * n_qpoints + 1:
                    self.q_files.append(line)
                    continue

                self.mu_star = float(split[0])
                break

        assert len(self.q_points) == n_qpoints
        assert len(self.q_files) == n_qpoints
        assert self.mu_star is not None


def main():
    # Setup arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", action="store",
                        help="The input file for lambda.x. "
                             "The first line is ignored - these parameters are worked out automatically.")
    parser.add_argument("--plot_a2f", action="store_true", help="Plot a2F(omega)")
    parser.add_argument("--a2f_smearing", action="store", required=False,
                        help="Manually specify the smearing width (in THz) used to evaluate a2F")

    # Parse arguments
    args = parser.parse_args()

    # Parse input file
    input_file = LambdaInput(args.input_file)

    # Parse elph input files
    lambda_input_files = [LambdaElphInput(x) for x in input_file.q_files]
    assert len(lambda_input_files) == len(input_file.q_points), \
        "Number of elph input files != number of q-points"

    # Grab number of double-delta smearing values
    nsigma = lambda_input_files[0].nsig
    for l in lambda_input_files:
        assert l.nsig == nsigma, \
            f"Inconsistent number of smearing values between {l.filename} and {lambda_input_files[0].filename}"

    # Get maximum phonon frequency
    max_freq = max(max(l.frequencies) for l in lambda_input_files)

    # Work out an appropriate smearing
    # (half of the averge distance between a mode and the next mode)
    all_freqs = []
    for l in lambda_input_files:
        all_freqs.extend(l.frequencies)
    smearing = np.mean(np.diff(sorted(all_freqs))) / 2.0

    # Manual smearing width override
    if args.a2f_smearing is not None:
        try:
            smearing = float(args.a2f_smearing) * THZ_TO_RY
        except ValueError:
            raise Exception(f"Could not parse smearing value from: '{args.a2f_smearing}'")

    # Generate frequency grid, initialize a2F and lambda arrays
    omega = np.linspace(0, max_freq + 10 * smearing, 10000)
    a2f = np.zeros((nsigma, len(omega)))
    lambdas = np.zeros((nsigma))

    # Grab normalized q-point weights
    weights = np.array([q[-1] for q in input_file.q_points])
    weights /= sum(weights)

    # Grab density of states @ Ef and degauss
    dosef = lambda_input_files[0].dos_at_efs
    deguass = lambda_input_files[0].sigmas
    for l in lambda_input_files:
        assert np.allclose(dosef, l.dos_at_efs), \
            f"Inconsistent DOS @ Ef from {l.filename} and {lambda_input_files[0].filename}"
        assert np.allclose(deguass, l.sigmas), \
            f"Inconsistent degauss from {l.filename} and {lambda_input_files[0].filename}"

    # Sum smeared Eliashberg function
    for w, l in zip(weights, lambda_input_files):

        # Frequency, lambda values for each smearing at that frequency
        for f, lams in zip(l.frequencies, l.mode_lambdas):

            assert len(lams) == nsigma, \
                f"Number of lambda values differ from number of deguass values in {l.filename}"

            # Evaluate guassian (see QE w0gauss method)
            x = (f - omega) / smearing
            sqrtpm1 = np.pi ** (-0.5)
            guassian = np.exp(-x ** 2) * sqrtpm1 / smearing

            # Sum contribution to a2F from this mode
            for n in range(nsigma):
                a2f[n] += w * lams[n] * f * 0.5 * guassian
                lambdas[n] += w * lams[n]

    # Evaluate omega log from frequencies and lambda values for each mode
    omega_logs = np.ones((nsigma))  # Initialized to 1, as is evaluated using logarithmic mean
    for w, l in zip(weights, lambda_input_files):  # Loop over q-points
        for f, lams in zip(l.frequencies, l.mode_lambdas):  # Loop over frequencies
            for n in range(nsigma):  # Loop over double-delta smearing values
                p = w * lams[n] / lambdas[n]  # \lambda_{qn} / \lambda
                omega_logs[n] *= f ** p

    # Convert frequencies to cm^{-1}
    omega *= RY_TO_CMM

    # Convert omega log to K
    omega_logs *= RY_TO_K

    # Print output
    print(f"{'degauss':>15} {'N(Ef)':>15} {'lambda':>15} {'omega_log':>15} {'T_c':>15}")
    for i in range(nsigma):
        ad_exp = np.exp(-1.04 * (1 + lambdas[i]) / (lambdas[i] - input_file.mu_star * (1 + 0.62 * lambdas[i])))
        tc = omega_logs[i] * ad_exp / 1.2
        print(f"{deguass[i]:>15.8f} {dosef[i]:>15.8f} {lambdas[i]:>15.8f} {omega_logs[i]:>15.8f} {tc:>15.8f}")

    # Plot a2F if requested
    if args.plot_a2f:
        import matplotlib.pyplot as plt

        for i in np.argsort(deguass):
            c = i / (len(deguass) - 1)
            c = (1 - c, c, 0)
            plt.plot(omega, a2f[i, :], color=c, label=f"deguass = {deguass[i]}")

        plt.xlim([0, max(omega)])
        plt.xlabel(r"$\omega$ (cm$^{-1}$)")
        plt.ylabel(r"$\alpha^2F(\omega)$")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    main()
