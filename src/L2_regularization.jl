function obj_f(f, p)
    return sum((p[1] * f - p[2]).^2)
end


function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::Type{linear_tikhonov}, order::Int=0)

    A = sparse([K; √(α) .* NMRInversions.Γ(size(K, 2), order)])
    b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    optf = Optimization.OptimizationFunction(obj_f, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, ones(size(K, 2)), (A,b), lcons=zeros(size(K, 2)), ucons=Inf * ones(size(K, 2)))
    f = OptimizationOptimJL.solve(prob, OptimizationOptimJL.LBFGS(), maxiters=5000, maxtime=100)
    r = K * f - g

    return f, r
end
