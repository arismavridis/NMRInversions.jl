
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

    Raw = readdlm(joinpath(directory , "T1IRT2.dat"), ' ')

    if size(Raw, 2) == 1
        Raw = readdlm(joinpath(directory ,"T1IRT2.dat"), ',')
    end
    Raw = transpose(complex.(Raw[:, 1:2:end], Raw[:, 2:2:end]))

    # Read experiment parameters
    acqu = readdlm(joinpath(directory , "acqu.par.bak"))
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
