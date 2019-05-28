# This script is to calculate the roots of the
# transcendental equation with the help of Newton-Raphson
# method. The transcendental equation is given by:
# z tanz = Bi;

# length of the symmetric square portion
L <- 0.05

# convection coefficient
h <- 100

# thermal conductivity
k <- 63.9

# Biot number
Bi <- (h*L)/k

# number of roots to be calculated
n <- 1:150

# initial guesses for those 100 roots
guess <- c(numeric(length(n))+((n-1)*pi))
solution <- guess

# tolerance for calculation
tolerance <- 1e-6

# the actual calculation
for (i in n) {
    # displace pi by a tenth to get an initial guess
    z <- (i-1)*pi + 0.1*pi
    # for each root, we need to iterate to get the closest approximation
    for (j in 1:1000) {
        # the function to be calculated
        func <- z*sin(z) - Bi*cos(z)
        # its derivative
        deriv <- (1 + Bi)*sin(z) + z*cos(z)
        # the slope
        dz <- - func/deriv
        # new guess based on the slope
        z <- z + dz
        # check if the solution is accurate enough
        if (abs(func) < tolerance) {
            break
        }
    }
    # diverging solution needs to be handled
    if (j > 1000) {
        print("Solution did not converge i = ", i)
    }
    # replace the initial guess with the actual solution
    solution[i] <- z
}

# make a data frame for simple comparison and plotting
transComparison <- data.frame(guess, solution)
