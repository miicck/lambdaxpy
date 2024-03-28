## lambdaxpy
A replacement for the lambda.x program shipped with quantum espresso. Improvements:

* An analytic expression for omega_log is used. This means that the following results no longer depend on the smeared a2F function (and, therefore, the smearing parameters):
  * lambda (this was already the case in lambda.x)
  * omega_log
  * The resulting Allen-Dynes Tc
* Sensible paramters for smearing a2F are chosen automatically from the input phonon frequencies:
  * No need to specify a maximum a2F frequency
  * No need to specify a smearing width (although you can if you like)
* Optional a2F plotting

# Installation
The following will install the python package, and run the tests

    ./install.sh

Have a look inside if you don't trust me.

# Usage
Where you would have done

    lambda.x < lambda.in > lambda.out

Now do

    lambdaxpy lambda.in > lambda.out

As was the case for lambda.x, this will also produce an alpha2F.dat file containing the eliashberg function at each degauss value.

Note that the first line of lambda.in (containing parameters for lambda.x) is ignored - these parameters are worked out automatically where needed to smear a2F, and are unnecassary to evaluate lambda, omega_log and the Allen-Dynes Tc (see above).
