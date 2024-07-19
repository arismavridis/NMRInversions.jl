
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


function import_spinsolve(directory::String=pick_folder(pwd()))

    Raw = readdlm(joinpath(directory, "T1IRT2.dat"), ' ')

    if size(Raw, 2) == 1
        Raw = readdlm(joinpath(directory, "T1IRT2.dat"), ',')
    end
    Raw = transpose(complex.(Raw[:, 1:2:end], Raw[:, 2:2:end]))

    # Read experiment parameters
    acqu = readdlm(joinpath(directory, "acqu.par.bak"))
    n_echoes = acqu[21, 3]
    t_echo = acqu[12, 3] * 1e-6
    τ_steps = acqu[36, 3]
    τ_min = acqu[20, 3] * 1e-3
    τ_max = acqu[19, 3] * 1e-3

    ## Make time arrays
    # Time array in direct dimension
    t_direct = collect(1:n_echoes) * t_echo

    # Time array in direct dimension
    if acqu[18, 3] == "yes" # if log spacing is selected, do log array
        t_indirect = exp10.(range(log10(τ_min), log10(τ_max), τ_steps))
    else                   # otherwise, do a linear array
        t_indirect = collect(range(τ_min, τ_max, τ_steps))
    end

    return t_direct, t_indirect, Raw

end

function readresults(file::String=pick_file(pwd()))

    open(file) do io

        readuntil(io, "Pulse Sequence : ")
        PulseSequence = readline(io)
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


function import_geospec(directory::String=pick_file(pwd()))

    open(directory) do io
        readuntil(io, "[Data]")

        data = readdlm(io, '\t', Float64, skipstart=2)

        return data
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

function autophase(Re, Im)

    ϕ_range = range(0, 2π, 20000)
    Re1_vs_φ = Re[1] .* cos.(ϕ_range) - Im[1] .* sin.(ϕ_range)
    ϕ = ϕ_range[argmax(Re1_vs_φ)]

    # optf = Optimization.OptimizationFunction(im_cost)
    # prob = Optimization.OptimizationProblem(optf, [ϕ₀], (Re, Im), lb=[ϕ₀ - 0.1], ub=[ϕ₀ + 0.1])
    # ϕ = OptimizationOptimJL.solve(prob, OptimizationOptimJL.ParticleSwarm([0], [2π], 10))[1]

    # cons(res, u, p) = (res .= p[1][1])
    # optf = Optimization.OptimizationFunction(im_cost, Optimization.AutoForwardDiff(); cons=cons)
    # prob = Optimization.OptimizationProblem(optf, [rand()], (Re, Im), lcons=[0.001], ucons=[(abs(Re[1]))])
    # ϕ = OptimizationOptimJL.solve(prob, OptimizationOptimJL.IPNewton())[1]

    Rₙ, Iₙ = phase_shift(Re, Im, ϕ)

    return Rₙ, Iₙ, ϕ
end


function testcorrection()
    # Create real and imaginary parts
    Re = exp.(-range(1, 20, 1000)) + randn(1000) .* 0.01
    Im = randn(1000) .* 0.01

    # Import data
    # data = import_geospec("/otherdata/9847zg/stratum_nmr/plug9_IR.txt")
    # Re = data[:, 3]
    # Im = data[:, 4]

    p1 = plot([Re, Im], label=["Original real" "Original Imaginary"])

    # Get them out of phase
    ϕd = rand() * 2π
    Re, Im = phase_shift(Re, Im, ϕd)
    
    p2 = plot([Re, Im], label=["Dephased real" "Dephased Imaginary"])
    
    Rₙ, Iₙ, ϕc = autophase(Re, Im)
    p3 = plot([Rₙ, Iₙ], label=["Corrected real" "Corrected Imaginary"])

    ϕ_range = range(0, 2π, 20000)
    Re1_vs_φ = Re[1] .* cos.(ϕ_range) - Im[1] .* sin.(ϕ_range)
    Im_sum_vs_φ = [im_cost([ϕ], (Re, Im)) for ϕ in ϕ_range]

    p4 = plot(ϕ_range, Re1_vs_φ, xlabel="ϕ", label="Re[1]")
    p4 = plot!(ϕ_range, (Im_sum_vs_φ ./maximum(Im_sum_vs_φ)).*maximum(Re1_vs_φ) , xlabel="ϕ", label="sum(Im.^2)", legend=:topleft)
    p4 = vline!(p4, [ϕc], label="corrected phase")
    display(plot(p1, p2, p3, p4))

    display("The correction error is $(2π - (ϕd + ϕc)) radians")

end

