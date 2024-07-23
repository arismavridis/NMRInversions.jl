function objective_f(c::AbstractVector, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    f = 0.5 * c' * ((p[1] * I + (G * p[3]')) * c) - c' * p[2]
end

function gradient_f(gradient, c, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    gradient .= (p[1] * I + G * p[3]') * c - p[2]
end

function hessian_f(H, c, p=(α, data, K))
    G = copy(p[3])
    G[:, (p[3]'*c.<=0)] .= 0
    H .= p[1] * I + G * p[3]'
end


function solve_regularization(K::AbstractMatrix, g::AbstractVector, α::Real, solver::Type{song}, order::Int=0)

    optf = Optimization.OptimizationFunction(objective_f, grad=gradient_f, hess=hessian_f)
    prob = Optimization.OptimizationProblem(optf, ones(size(K, 1)), (α, g, K))
    c = OptimizationOptimJL.solve(prob, OptimizationOptimJL.NewtonTrustRegion(), x_tol=1e-8, maxiters=5000, maxtime=100)
    f = vec(max.(0, (K' * c)))
    r = K * f - g

    return f, r

end



