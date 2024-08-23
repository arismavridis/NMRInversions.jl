
function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::Type{optim_nnls}, order::Int=0)

    A = sparse([K; √(α) .* NMRInversions.Γ(size(K, 2), order)])
    b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    f = solve_nnls(A, b)

    r = K * f - g

    return f, r
end


function obj_f(x, p)
    return sum((p[1] * x - p[2]).^2)
end


function solve_nnls(A::AbstractMatrix, b::AbstractVector)

    optf = Optimization.OptimizationFunction(obj_f, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, ones(size(A, 2)), (A, b), lb=zeros(size(A, 2)), ub=Inf * ones(size(A, 2)))
    x = OptimizationOptimJL.solve(prob, OptimizationOptimJL.LBFGS(), maxiters=5000, maxtime=100)

    return x
end