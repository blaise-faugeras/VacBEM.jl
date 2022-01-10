function gauleg(a,b,n)
#==
    1D Gauss-Legendre quadrature rule.

    quadrule = gauleg(a,b,n) computes the n-point Gauss-Legendre
    quadrature rule on the interval [a,b] up to mashine precision.

    Note that all quadrature rules obtained from gauleg are of order 2*n-1.

   The output named tuple quadrule contains the following fields:
    w n-by-1 vector specifying the weights of the quadrature rule.
    X n-by-1 vector specifying the abscissae of the quadrature rule.

   Example:
    quad = gauleg(0,1,10); 
    quad.x 
    quad.w
=#  
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n==1 
        x = 0; w = 2;
    else
        d = zeros(n)
        c = zeros(n-1)
        for i in 1:(n-1)
            c[i]=i/sqrt(4*i*i-1)
        end
        J=LinearAlgebra.SymTridiagonal(d,c)
        x=LinearAlgebra.eigvals(J)
        vec=LinearAlgebra.eigvecs(J)
        w=2*vec[1,:].*vec[1,:]
        
    end
    
    # The above rule computed a quadrature rule for the reference interval.
    # The reference interval is always [-1,1], and if [a,b] != [-1,1] then we
    # re-map the abscissas (x) and rescale the weights.

    xm = (b+a)/2; # midpoint
    xl = (b-a)/2; # area of the requested interval / area of the reference interval

    x = xm .+ xl*x;
    w = w * xl;

    return (x=x,w=w) #named tuple
end
