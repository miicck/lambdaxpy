## lambdaxpy

A replacement for the lambda.x program shipped with quantum espresso. Improvements:

* An analytic expression for omega_log is used (see below). This means that the following results no longer depend on the smeared a2F function (and, therefore, the smearing parameters):
  * lambda (this was already the case in lambda.x)
  * omega_log
  * The resulting Allen-Dynes Tc

![image](https://github.com/miicck/lambdaxpy/assets/8690175/fa762527-8e6e-4060-8932-f6f0be6166ec)

* Sensible paramters for smearing a2F are chosen automatically from the input phonon frequencies:
  * No need to specify a maximum a2F frequency
  * No need to specify a smearing width (although you can if you like with --a2f_smearing)
* Deals with some numerical instabilities in lambda.x (see blue lines below)

![image](https://github.com/miicck/lambdaxpy/assets/8690175/3396470e-f99c-4a06-ab2b-f6564e9f4bd3)

* Optional a2F plotting with --plot_a2f, including cumulative values for the above quantitites:

![image](https://github.com/miicck/lambdaxpy/assets/8690175/fa6798a0-2e14-4572-8136-6f79d4c1b3cf)

# Installation
Clone this repository and cd into the resulting directory. The following will install the python package, and run the tests

    ./install.sh

Have a look inside if you don't trust me.

# Usage
Where you would have done

    lambda.x < lambda.in > lambda.out

Now do

    lambdaxpy lambda.in > lambda.out

As was the case for lambda.x, this will also produce an alpha2F.dat file containing the eliashberg function at each degauss value.

Note that the first line of lambda.in (containing parameters for lambda.x) is ignored - these parameters are worked out automatically where needed to smear a2F, and are unnecassary to evaluate lambda, omega_log and the Allen-Dynes Tc (see above).
