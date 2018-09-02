
function lsgd(A, b, mu, x0, nIters)
#
# Syntax: x = lsgd(A, b, mu, x0, nIters)
#
# Inputs: A is an m x n matrix
# b is a vector of length m
# mu is the step size to use, and must satisfy
# 0 < mu < 2 / norm(A)^2 to guarantee convergence
# x0 is the initial starting vector (of length n) to use
# nIters is the number of iterations to perform
#
# Outputs: x is a vector of length n containing the approximate solution
#
# Description: Performs gradient descent to solve the least squares problem
#
# \min x \|b - A x\|_2
#

    
    b  =  vec(b)
    x0  =  vec(x0)
    
    
    x  =  x0
    for k in 1:nIters
        x  -=  mu*(A'*(A*x-b))
    end
    
    return x 
end
