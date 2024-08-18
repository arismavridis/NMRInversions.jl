struct input1D
    exptype::Type{<:inversion1D}
    x::AbstractVector{<:Real}
    y::AbstractVector
end

struct input2D
    exptype::Type{<:inversion2D}
    x_direct::AbstractVector{<:Real}
    x_indirect::AbstractVector{<:Real}
    data::AbstractMatrix
end

struct invres1D
    exptype::Type{<:inversion1D}
    X::Vector
    f::Vector
    r::Vector
    SNR::Real
    alpha::Real
    x_fit::Vector
    y_fit::Vector
end

struct invres2D
    exptype::Type{<:inversion2D}
    X_dir::Vector
    X_indir::Vector
    f::Vector
    r::Vector
    SNR::Real
    alpha::Real
end

function import_1D(file=pick_file(pwd()))
    data = readdlm(file, ',')
    x = vec(data[:, 1])
    y = vec(data[:, 2])
    return x, y
end

function create_decay(exptype="T1", ;
    coeffs=[1 1], noise=0.01,
    tlb=0.01, tub=5, points=32,
    logspace=false)
    coeffs = coeffs[:, :]

    if exptype == "T1"
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)

    elseif exptype == "T2"
        kernel_eq = (t, T) -> exp(-t / T)

    end

    #Create a time vector

    t = []
    if logspace == true
        t = exp10.(range(log10(tlb), log10(tub), points))
    elseif logspace == false
        t = collect(range(tlb, tub, points))
    end

    #Create a synthetic signal
    signal = vec(sum(coeffs[:, 1] .* kernel_eq.(t', coeffs[:, 2]), dims=1)')
    display(typeof(signal))
    signal = signal .+ (noise * maximum(signal)) .* randn(length(t))

    display(scatter(t, signal, legend=false))
    return t, signal
end

function read_acqu(filename, parameter)

    p = ""
    open(filename) do io
        readuntil(io, parameter*" = ")
        p = readline(io)
    end
    return replace(p, "\"" => "")
end


function import_spinsolve(directory::String=pick_folder(pwd()))

    directory = abspath(directory)
    cd(directory)

    # Read experiment parameters
    acqu = readdlm(joinpath(directory, "acqu.par.bak"))
    n_echoes = acqu[21, 3]
    t_echo = acqu[12, 3] * 1e-6
    τ_steps = acqu[36, 3]
    τ_min = acqu[20, 3] * 1e-3
    τ_max = acqu[19, 3] * 1e-3
    experiment = acqu[13, 3]

    Raw = readdlm(joinpath(directory, experiment * ".dat"), ' ')

    if size(Raw, 2) == 1
        Raw = readdlm(joinpath(directory, experiment * ".dat"), ',')
    end

    Data = collect(transpose(complex.(Raw[:, 1:2:end], Raw[:, 2:2:end])))

    ## Make time arrays
    # Time array in direct dimension
    t_direct = collect(1:n_echoes) * t_echo

    # Time array in direct dimension
    if acqu[18, 3] == "yes" # if log spacing is selected, do log array
        t_indirect = exp10.(range(log10(τ_min), log10(τ_max), τ_steps))
    else                   # otherwise, do a linear array
        t_indirect = collect(range(τ_min, τ_max, τ_steps))
    end

    if experiment == "T1IRT2"
        exptype = IRCPMG
    end

    return input2D(exptype, t_direct, t_indirect, Data)

end

function writeresults(name::String, res::invres2D)

    open(name, "w") do io
        write(io, "Pulse Sequence : " * string(res.exptype) * "\n")
        write(io, "SNR : " * string(res.SNR) * "\n")
        write(io, "alpha : " * string(res.alpha) * "\n")
        write(io, "Direct Dimension : " * join(res.X_dir, ", ") * "\n")
        write(io, "Indirect Dimension : " * join(res.X_indir, ", ") * "\n")
        write(io, "Inversion Results : " * join(res.f, ", ") * "\n")
        write(io, "Residuals : " * join(res.r, ", ") * "\n")
    end

end

function readresults(file::String=pick_file(pwd()))

    open(file) do io

        readuntil(io, "Pulse Sequence : ")
        PulseSequence = eval(Meta.parse(readline(io)))
        readuntil(io, "SNR : ")
        SNR = parse.(Float64, readline(io))
        readuntil(io, "alpha : ")
        α = parse.(Float64, readline(io))
        readuntil(io, "Direct Dimension : ")
        dir = parse.(Float64, split(readline(io), ','))
        readuntil(io, "Indirect Dimension : ")
        indir = parse.(Float64, split(readline(io), ','))
        readuntil(io, "Inversion Results : ")
        f = parse.(Float64, split(readline(io), ','))
        readuntil(io, "Residuals : ")
        r = parse.(Float64, split(readline(io), ','))

        return invres2D(PulseSequence, dir, indir, f, r, SNR, α)
    end
end


function entropy(u::AbstractArray, p::Tuple)

    phc0 = u[1]
    phc1 = u[2]
    R⁰ = p[1]
    I⁰ = p[2]
    γ = p[3]

    n = length(R⁰)
    ϕ = phc0 .+ phc1 .* collect(range(1, n) ./ n)

    R, _ = phase_shift(R⁰, I⁰, ϕ)

    h = abs.(diff(R)) ./ sum(abs.(diff(R)))

    return -sum(h .* log.(h)) + γ * sum((R .^ 2)[findall(x -> x < 0, R)])

end


function minimize_entropy(Re, Im, γ)

    optf = Optimization.OptimizationFunction(entropy)
    prob = Optimization.OptimizationProblem(optf, [0.1, 0.1], (Re, Im, γ))# lb=[0.0001, 0.0001], ub=[2π, 2π])

    uopt = OptimizationOptimJL.solve(prob, OptimizationOptimJL.SimulatedAnnealing())

    n = length(Re)
    ϕ = uopt[1] .+ uopt[2] .* collect(range(1, n) ./ n)

    Rₙ, Iₙ = phase_shift(Re, Im, ϕ)

    p = plot(a[:, 1], Rₙ, label="Real")
    p = plot!(a[:, 1], Iₙ, label="Imaginary")
    display(p)

    return Rₙ, Iₙ

end


function im_cost(u, p)  # sum of squares of the imaginary Part

    Re = p[1]
    Im = p[2]
    ϕ = u[1]

    _, Iₙ = phase_shift(Re, Im, ϕ)
    return sum(Iₙ .^ 2)

end


function phase_shift(Re, Im, ϕ)

    Re_new = Re .* cos(ϕ) - Im .* sin(ϕ)
    Im_new = Im .* cos(ϕ) + Re .* sin(ϕ)

    return Re_new, Im_new
end


function autophase(re, im, startingpoint::Real) # startingpoint is a value between -1 and 1

    ϕ_range = range(0, 2π, 500)
    Re1_vs_φ = re[1] .* cos.(ϕ_range) - im[1] .* sin.(ϕ_range)

    ϕ₀ = ϕ_range[argmin(abs.(Re1_vs_φ .- startingpoint * maximum(Re1_vs_φ)))]

    optf = Optimization.OptimizationFunction(im_cost, Optimization.AutoForwardDiff())
    prob = Optimization.OptimizationProblem(optf, [ϕ₀], (re, im), lb=[ϕ₀ - 1], ub=[ϕ₀ + 1], x_tol=1e-5)
    ϕ = OptimizationOptimJL.solve(prob, OptimizationOptimJL.BFGS())[1]

    Rₙ, Iₙ = phase_shift(re, im, ϕ)

    return Rₙ, Iₙ, ϕ
end

function import_geospec(filedir::String=pick_file(pwd()))

    cd(dirname(filedir))

    data = []
    pulse_sequence_number::Int16 = 0
    dimensions = [0, 0]

    open(filedir) do io

        readuntil(io, "TestType=")
        pulse_sequence_number = parse(Int16, readline(io))
        readuntil(io, "Dimensions=")
        dimensions .= parse.(Int16, split(readline(io), ','))
        readuntil(io, "[Data]")
        data = readdlm(io, '\t', Float64, skipstart=2)
    end

    y_re = data[:, 3]
    y_im = data[:, 4]


    typedict = Dict(
        3 => CPMG,
        7 => IR,
        105 => PFG,
        106 => IRCPMG,
        108 => PFGCPMG
        # ,110 => CPMGCPMG
    )

    exptype = typedict[pulse_sequence_number]

    if exptype in [IR, IRCPMG]
        y_re, y_im, ϕ = autophase(y_re, y_im, -1)
    else
        y_re, y_im, ϕ = autophase(y_re, y_im, 1)
    end

    display("Data phase corrected by $(round(ϕ,digits=3)) radians.")

    if exptype == IRCPMG

        return input2D(IRCPMG, data[1:dimensions[1], 1] .* (1 / 1000), data[1:dimensions[1]:end, 2] .* (1 / 1000), reshape(complex.(y_re, y_im), dimensions[1], dimensions[2]))

    elseif exptype == PFGCPMG

        return input2D(PFGCPMG, data[1:dimensions[1], 1], data[1:dimensions[1]:end, 2] .* (1 / 1000), reshape(complex.(y_re, y_im), dimensions[1], dimensions[2]))

    elseif exptype == PFG

        return input1D(exptype, data[:, 1], complex.(y_re, y_im))

    elseif exptype in [IR, CPMG]

        return input1D(exptype, data[:, 1] .* (1 / 1000), complex.(y_re, y_im)) # Converts time to seconds
    end
end
