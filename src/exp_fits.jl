using Optimization, OptimizationOptimJL, Plots

# x = collect(range(0, 5, 100))
# y = 2.0 * exp.(-(1/5) .* x) + 3.0 * exp.(-(1/0.1) .* x)
# y = 2.5 * exp.(-(1 / 0.5) .* x) + 0.04 .* randn(length(x))
# plot(x, y)


function single_exp_loss(u, p)
    x = p[1]
    y = p[2]
    fit = u[1] * exp.(-(1 / u[2]) * x)
    return sum((fit .- y) .^ 2)
end

function double_exp_loss(u, p)
    x = p[1]
    y = p[2]
    fit = u[1] * (u[2] * exp.(-(1 / u[3]) * x) 
        + (1 - u[2]) * exp.(-(1 / u[4]) * x))
    return sum((fit .- y) .^ 2)
end

function triple_exp_loss(u, p)
    x = p[1]
    y = p[2]
    fit = u[1] * exp.(-(1 / u[4]) * x) 
        + u[2] * exp.(-(1 / u[5]) * x) 
        + u[3] * exp.(-(1 / u[6]) * x)
    return sum((fit .- y) .^ 2)
end

function m_exp(seq,x,y,p0=[1.0, 1.0])
    optf = Optimization.OptimizationFunction(single_exp_loss, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(single_exp_loss, p0, [x, y])
    OptimizationOptimJL.solve(prob, OptimizationOptimJL.NelderMead(), maxiters=5000, maxtime=100)
end

function d_exp(seq,x,y,p0=[1.0, 0.9, 1.0, 1.0])
    optf = Optimization.OptimizationFunction(double_exp_loss, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(double_exp_loss, p0, [x, y])
    OptimizationOptimJL.solve(prob, OptimizationOptimJL.NelderMead(), maxiters=5000, maxtime=100)
end
