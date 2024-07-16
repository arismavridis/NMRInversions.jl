
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
        readuntil(io,"[Data]")

        data = readdlm(io,'\t',Float64,skipstart=2)

        p = plot(data[:,1],data[:,3],label="real")
        p = plot!(data[:,1],data[:,4],label="imaginary")
        display(p)
        return data
    end

end



function autophase_entropy(R, I)

    model = Model(HiGHS.Optimizer)


    n = length(R)

    @variable(model, phc0)
    @variable(model, phc1)
    @variable(model, ϕ)

    @objective(model, Min, -sum(h .* ln.(h)))

    @constraint(model, h .== abs(R))
    @constrain(model, ϕ .== phc0 .+ phc1 .* collect(range(1,n) ./ n))

    ϕ = value.(ϕ)
    
    Rₙ .= R * cos(ϕ) - I * sin(ϕ) 
    Iₙ .= I * cos(ϕ) + R * sin(ϕ) 

    return Rₙ, Iₙ
end

function autophase_integral(R,I)

    model = Model(HiGHS.Optimizer)
    @variable(model, ϕ)

end

