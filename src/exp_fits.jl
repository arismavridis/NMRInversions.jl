export expfit_struct
"""
Structure containing information about multiexponential fits.

The fields are as follows:
- `seq`: The pulse sequence
- `x` : The  

"""
struct expfit_struct
    seq::Type{<:NMRInversions.pulse_sequence1D}
    x::Vector
    y::Vector
    u::Vector
    u0::Vector
    r::Vector
    eq::String
    eqn::String
end

export mexp
"""
    mexp(u, x)
Construct a multi-exponential function from a sequence of parameters u, and evaluate it at the values x. \n
The `u` vector should be of the form [a1, b1, a2, b2, ...], 
where a's are the amplitudes and b's are the reciprocals of the decay constants.
The length of `u` must be even.

Examples: \n
mexp([a,b] , x) = a * exp.( (1/b) * x) \n
mexp([a,b,c,d] , x) = a * exp.( (1/b) * x) + c * exp.((1/d) * x) \n
mexp([a,b,c,d,e,f] , x) = a * exp.( (1/b) * x) + c * exp.((1/d) * x) + e * exp.((1/f) * x) \n
. . .

"""
mexp(u, x) = sum([u[i] * exp.(-(1 / u[i+1]) * x) for i in 1:2:length(u)])

"Loss function for multi-exponential fit, returns the sum of squared differences between the data and the model."
function mexp_loss(u, p)
    x = p[1]
    y = p[2]
    return sum((mexp(u, x) .- y) .^ 2)
end


export expfit
"""
    expfit(n, seq, x, y; u0, solver=BFGS())
Fit an exponential function to the data x, y using n exponential terms. \n

"""
function expfit(n::Int, seq, x::Vector, y::Vector;
    u0=rand(2n), solver=OptimizationOptimJL.BFGS())

    yfit = copy(y)
    # Transform recovery to decay if needed
    if seq == IR
        yfit .= (y[end] .- y) ./ (2 .* y[end])
    elseif seq == SR
        y .= y[end] .- y
    end

    # Solve the optimization
    optf = Optimization.OptimizationFunction(mexp_loss, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, u0, [x, y], lb=zeros(length(u0)), ub=Inf .* ones(length(u0)))
    u = OptimizationOptimJL.solve(prob, solver, maxiters=5000, maxtime=100)

    # Determine what's the x-axis of the seq (time or bfactor)
    seq == NMRInversions.PFG ? x_ax = "b" : x_ax = "t"

    # Get the fit's equation as a string
    eq = join([(i == 1 ? "" : " + ") * string(round(u[i], sigdigits=2)) * " * exp($x_ax/" * string(round(u[i+1], sigdigits=2)) * ")" for i in 1:2:length(u)])

    # Normalised version of the fit equation
    un = sum([u[i] for i in 1:2:length(u)])
    eqn = string(round(un, sigdigits=2)) * " * (" *
          join([(i == 1 ? "" : " + ") * string(round(u[i] / un, sigdigits=2)) *
                " * exp($x_ax/" * string(round(u[i+1], sigdigits=2)) * ")" for i in 1:2:length(u)]) * ")"

    # Calculate the residuals
    r = mexp(u, x) - y

    if n == 1
        return expfit_struct(seq, x, y, u, u0, r, eq, eq)
    else
        return expfit_struct(seq, x, y, u, u0, r, eq, eqn)
    end

end

function expfit(n::Int, res::NMRInversions.input1D; kwargs...)
    expfit(n, res.seq, res.x, res.y, kwargs...)
end

