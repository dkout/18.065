function lsngd(A, b, mu, x0, nIters)
#
# Syntax: x = lsngd(A, b, mu, x0, nIters)
#
# Inputs: A is an m x n matrix
# b is a vector of length m
# mu is the step size to use, and must satisfy
# 0 < mu < 1 / norm(A)^2 to guarantee convergence
# x0 is the initial starting vector (of length n) to use
# nIters is the number of iterations to perform
#
# Outputs: x is a vector of length n containing the approximate solution
#
# Description: Performs Nesterov-accelerated gradient descent to solve
# the least squares problem
# \min x \|b - A x\|_2
#
    b  =  vec(b)
    x0  =  vec(x0)

    
    t  =  0
    xLast  =  x0
    x  =  x0
    
    for _ =  1:nIters
        
        tLast  =  t
        t  =  0.5 * (1  +  sqrt(1  +  4 * t^2))
        z  =  x  +  ((tLast  -  1)  /  t)*(x  -  xLast)
        
        xLast  =  x
        x  =  z  -  mu*(A'*(A*z-b))
        
    end
    
    return x
end