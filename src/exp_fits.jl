export expfit_struct
"""
Structure containing information about multiexponential fits.

The fields are as follows:
- `seq`: The pulse sequence
- `x` : The x acquisition values (e.g. time for relaxation or b-factor for diffusion).
- `y` : The y acquisition values.
- `u` : The fitted parameters for the `mexp` function.
- `u0` : The initial parameters for the `mexp` function.
- `r` : The residuals.
- `eq` : The equation of the fitted function.
- `eqn` : The equation of the initial function.

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
    mexp(seq, u, x)
Construct a multi-exponential function from parameters u, and evaluate it at the values x. \n
The `u` vector should be of the form [a1, b1, a2, b2, ...], 
where a's are the amplitudes and b's are the reciprocals of the decay constants.
The length of `u` must be an even number.

Examples: \n
mexp(CPMG, [a,b] , x) = a * exp.( (1/b) * x) \n
mexp(CPMG, [a,b,c,d] , x) = a * exp.( (1/b) * x) + c * exp.((1/d) * x) \n
mexp(CPMG, [a,b,c,d,e,f] , x) = a * exp.( (1/b) * x) + c * exp.((1/d) * x) + e * exp.((1/f) * x) \n
(where a,b,c,d,e,f are all numbers)
. . .

"""
function mexp(seq, u, x)

    if seq == IR
        un = sum([u[i] for i in 1:2:length(u)])
        return un .- 2 .* mexp(u,x)

    elseif seq == SR
        un = sum([u[i] for i in 1:2:length(u)])
        return un .-  mexp(u,x)

    elseif seq in [CPMG, PFG]
        return mexp(u, x)

    end

end
mexp(u, x) = sum([u[i] * exp.(-(1 / u[i+1]) * x) for i in 1:2:length(u)])


"Loss function for multi-exponential fit, returns the sum of squared differences between the data and the model."
function mexp_loss(u, p)
    x = p[1]
    y = p[2]
    seq = p[3]
    return sum((mexp(seq, u, x) .- y) .^ 2)
end


export expfit
"""
    expfit(n, seq, x, y; u0, solver=BFGS())
Fit an exponential function to the data x, y using n exponential terms. \n

"""
function expfit(n::Int, seq::Type{<:NMRInversions.pulse_sequence1D}, x::Vector, y::Vector;
    u0::Vector=rand(2n), solver=OptimizationOptimJL.BFGS())

    u0 = Float64.(u0)

    # Solve the optimization
    optf = Optimization.OptimizationFunction(mexp_loss, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, u0, (x, y, seq), lb=zeros(length(u0)), ub=Inf .* ones(length(u0)))
    u = OptimizationOptimJL.solve(prob, solver, maxiters=5000, maxtime=100)

    # Determine what's the x-axis of the seq (time or bfactor)
    seq == NMRInversions.PFG ? x_ax = "b" : x_ax = "t"

    un = sum([u[i] for i in 1:2:length(u)]) 

    # Get the fit's equation as a string
    eq = join([(i == 1 ? "" : " + ") * string(round(u[i], sigdigits=2)) * " * exp($x_ax/" * string(round(u[i+1], sigdigits=2)) * ")" for i in 1:2:length(u)])

    # Normalised version of the fit equation
    eqn = string(round(un, sigdigits=2)) * " * (" * join([(i == 1 ? "" : " + ") * string(round(u[i] / un, sigdigits=2)) * " * exp($x_ax/" * string(round(u[i+1], sigdigits=2)) * ")" for i in 1:2:length(u)]) * ")"

    if seq == IR
        unstr = string(round(un,sigdigits = 2))
        eq = unstr * " - 2 * (" * eq * ")" 
        eqn = unstr * " - 2 * " *  eqn   

    elseif seq == SR
        unstr = string(round(un,sigdigits = 2))
        eq = unstr * " -" * "(" * eq * ")" 
        eqn = unstr * " -" *  eqn   
    end

    # Calculate the residuals
    r = mexp(seq, u, x) - y

    if n == 1
        return expfit_struct(seq, x, y, u, u0, r, eq, eq)
    else
        return expfit_struct(seq, x, y, u, u0, r, eq, eqn)
    end

end

function expfit(n::Int, res::NMRInversions.input1D; kwargs...)
    expfit(n, res.seq, res.x, res.y, kwargs...)
end

