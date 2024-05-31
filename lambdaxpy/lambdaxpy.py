from typing import Iterable, Tuple
import os
import numpy as np
import argparse

# Constants
RY_TO_CMM = 109736.75775046606
RY_TO_K = 157887.6633481157
RY_TO_THZ = 3289.8441866350436
THZ_TO_RY = 1 / RY_TO_THZ
CMM_TO_K = 1.4387855681606843


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

            if f < 1e-10:
                # Ignore -ve, or zero, frequency contributions
                # to lambda or eliashberg function
                continue

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

            if f < 1e-10:
                # Ignore -ve, or zero, frequency contributions
                # to omega log
                continue

            for n in range(nsigma):  # Loop over double-delta smearing values
                if lambdas[n] > 1e-10:
                    p = w * lams[n] / lambdas[n]  # \lambda_{qn} / \lambda
                    omega_logs[n] *= f ** p
                else:
                    # Total lambda is zero
                    omega_logs[n] = 0.0

    # Create alpha2F.dat file
    with open("alpha2F.dat", "w") as f:
        f.write(f"{'# E(THz)':<10} " + " ".join(f"{x:>10.5f}" for x in deguass) + "\n")
        for i, w in enumerate(omega):
            f.write(f"{w * RY_TO_THZ:>10.5f} " + " ".join(f"{x:>10.5f}" for x in a2f[:, i]) + "\n")

    # Convert frequencies to cm^{-1}
    omega *= RY_TO_CMM

    # Convert omega log to K
    omega_logs *= RY_TO_K

    def allen_dynes_exponent(lam: float) -> float:
        denom = lam - input_file.mu_star * (1 + 0.62 * lam)
        return np.exp(-1.04 * (1 + lam) / denom) if denom > 0 else 0.0

    # Print output
    print(f"{'degauss':>15} {'N(Ef)':>15} {'lambda':>15} {'omega_log':>15} {'T_c':>15}")
    for i in range(nsigma):
        tc = omega_logs[i] * allen_dynes_exponent(lambdas[i]) / 1.2
        print(f"{deguass[i]:>15.8f} {dosef[i]:>15.8f} {lambdas[i]:>15.8f} {omega_logs[i]:>15.8f} {tc:>15.8f}")

    # Plot a2F if requested
    if args.plot_a2f:
        import matplotlib.pyplot as plt

        deguass_colors = {}
        for n in np.argsort(deguass):
            c = n / (len(deguass) - 1)
            c = (1 - c, c, 0)
            deguass_colors[n] = c

        plt.subplot(338)
        ddg = deguass[1] - deguass[0]
        for n, s in enumerate(deguass):
            plt.axhline(s, color=deguass_colors[n])
            plt.annotate(f"degauss = {s}", (0, s + ddg / 4))
        plt.ylabel("degauss")
        plt.ylim(min(deguass) - ddg, max(deguass) + ddg)
        plt.axis("off")

        for n in np.argsort(deguass):

            plt.subplot(311)
            plt.plot(omega, a2f[n, :], color=deguass_colors[n])
            plt.xlim([0, max(omega)])
            plt.xlabel(r"$\omega$ (cm$^{-1}$)")
            plt.ylabel(r"$\alpha^2F(\omega)$")

            cum_lam = [0]
            cum_omega_log_exp = [0]
            cum_omega_log = [0]
            cum_tc = [0]

            for i in range(1, len(omega)):
                # d omega for this integration rectangle
                d_omega = omega[i] - omega[i - 1]

                # averge a2f and omega values for this rectangle
                a_av = (a2f[n, i] + a2f[n, i - 1]) / 2.0
                w_av = (omega[i] + omega[i - 1]) / 2.0

                # increment in lambda, expontential appearing in expression for omega log
                d_lambda = 2 * d_omega * a_av / w_av
                d_omega_log_exp = 2 * d_omega * np.log(w_av) * a_av / w_av

                # Update cumulative values
                cum_lam.append(cum_lam[-1] + d_lambda)
                cum_omega_log_exp.append(cum_omega_log_exp[-1] + d_omega_log_exp)
                cum_omega_log.append(CMM_TO_K * np.exp(cum_omega_log_exp[-1] / cum_lam[-1]))
                cum_tc.append(cum_omega_log[-1] * allen_dynes_exponent(cum_lam[-1]) / 1.2)

            plt.subplot(334)
            plt.plot(omega, cum_lam, color=deguass_colors[n])
            plt.xlabel(r"$\omega$ (cm$^{-1}$)")
            plt.ylabel(r"Cumulative $\lambda$")

            plt.subplot(335)
            plt.plot(omega, cum_omega_log, color=deguass_colors[n])
            plt.xlabel(r"$\omega$ (cm$^{-1}$)")
            plt.ylabel(r"Cumulative $\omega_{log}$")

            plt.subplot(336)
            plt.plot(omega, cum_tc, color=deguass_colors[n])
            plt.xlabel(r"$\omega$ (cm$^{-1}$)")
            plt.ylabel(r"Cumulative $T_c$")

        plt.show()


if __name__ == "__main__":
    main()
