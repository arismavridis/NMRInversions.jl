struct svd_kernel_struct
    K::AbstractMatrix
    g::AbstractVector
    U::AbstractMatrix
    S::AbstractVector
    V::AbstractMatrix
end

"""
Create a kernel for the inversion of 1D data.
x is the experiment x axis (time or b factor etc.)
X is the range for the output x axis (T1, T2, D etc.)

If data vector is provided, SVD is performed on the kernel.
If data vector is complex, the SNR is calculated and the SVD is automatically truncated accordingly.
"""
function create_kernel(exptype::Type{<:inversion1D}, x::Vector, X::Vector)
    if exptype == IR
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)
    elseif exptype in [CPMG, PFG]
        kernel_eq = (t, T) -> exp(-t / T)
    end

    return kernel_eq.(x, X')
end

function create_kernel(exptype::Type{<:inversion1D}, x::Vector, X::Vector, g::Vector{<:Real})
    if exptype == IR
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)
    elseif exptype in [CPMG, PFG]
        kernel_eq = (t, T) -> exp(-t / T)
    end

    usv = svd(kernel_eq.(x, X'))

    K_new = Diagonal(usv.S) * usv.V'
    g_new = usv.U' * g
    
    return svd_kernel_struct(K_new, g_new, usv.U, usv.S, usv.V)

end

function create_kernel(exptype::Type{<:inversion1D}, x::Vector, X::Vector, g::Vector{<:Complex})
    if exptype == IR
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)
    elseif exptype in [CPMG, PFG]
        kernel_eq = (t, T) -> exp(-t / T)
    end
   
    SNR = calc_snr(g)
    usv = svd(kernel_eq.(x, X'))
    indices = findall(i -> i .> (1 / SNR), usv.S) # find elements in S12 above the noise threshold

    display("SVD truncated to $(length(indices)) singular values out of $(length(usv.S))")

    U = usv.U[:, indices]
    S = usv.S[indices]
    V = usv.V[:, indices]

    K_new = Diagonal(S) * V'
    g_new = U' * real.(g)
    
    return svd_kernel_struct(K_new, g_new, U, S, V)

end

"""
Calculate the Signal-to-Noise Ratio (SNR) from complex data,
where the real part is mostly signal and the imaginary part is mostly noise.
σ_n is the STD of the latter half of the imaginary signal 
(former half might contain signal residues as well) 
"""
function calc_snr(data::AbstractMatrix{<:Complex}) # For matrix data

    real_data = real.(data)
    imag_data = imag.(data)

    noise = imag_data[(floor(Int64, (size(imag_data, 1) / 2))):end, :]
    σ_n = sqrt(sum((noise .- sum(noise) / length(noise)) .^ 2) / (length(noise) - 1))
    SNR = maximum(abs.(real_data)) / σ_n

    return SNR
end
function calc_snr(data::AbstractVector{<:Complex}) # For vector data

    real_data = real.(data)
    imag_data = imag.(data)

    noise = imag_data[(floor(Int64, (size(imag_data, 1) / 2))):end]
    σ_n = sqrt(sum((noise .- sum(noise) / length(noise)) .^ 2) / (length(noise) - 1))
    SNR = maximum(abs.(real_data)) / σ_n

    return SNR
end


"""
Create a kernel for the inversion of 2D data.
t_direct is the direct dimension acquisition parameter
t_indirect is the indirect dimension acquisition parameter
Raw is the 2D data matrix of complex data
"""
function create_kernel(exptype::Type{<:inversion2D},
    x_direct::AbstractVector, x_indirect::AbstractVector,
    X_direct::AbstractVector, X_indirect::AbstractVector,
    Data::AbstractMatrix{<:Complex})

    G = real.(Data)
    SNR = calc_snr(Data)

    if SNR < 1000
        @warn("The SNR is $(round(SNR, digits=1)), which is below the recommended value of 1000. 
            Consider running experiment with more scans.")
    end

    # Generate Kernels
    if exptype == IRCPMG
        K_dir = create_kernel(CPMG, x_direct, X_direct)
        K_indir = create_kernel(IR, x_indirect, X_indirect)
    end

    ## Perform SVD truncation
    usv_dir = svd(K_dir) #paper (13)
    usv_indir = svd(K_indir) #paper (14)

    # finding which singular components are contributed from K1 and K2
    S21 = usv_indir.S * usv_dir.S' #Using outer product instead of Kronecker, to make indices more obvious
    indices = findall(i -> i .> (1 / SNR), S21) # find elements in S12 above the noise threshold
    s̃ = S21[indices]
    ñ = length(s̃)

    si = (first.(Tuple.(indices)))  # direct dimension singular vector indices
    sj = (last.(Tuple.(indices)))   # indirect dimension singular vector indices

    g̃ = diag(usv_indir.U[:, si]' * G' * usv_dir.U[:, sj])

    Ũ₀ = Array{Float64}(undef, 0, 0) # change that to the actual U at some point, this is just a placeholder

    V1t = repeat(usv_dir.V[:, sj], size(usv_indir.V, 1), 1)
    V2t = reshape(repeat(usv_indir.V[:, si]', size(usv_dir.V, 1), 1), ñ, size(V1t, 1))'
    Ṽ₀ = V1t .* V2t
    K̃₀ = Diagonal(s̃) * Ṽ₀'

    return svd_kernel_struct(K̃₀, g̃, Ũ₀, s̃, Ṽ₀)
end


## Data compression

function compress_data(t_direct::AbstractVector, G::AbstractMatrix, bins::Int=64)

    # Compress direct dimension to the length of bins

    # Use logarithmically spaced windows for a window average
    windows = zeros(bins)

    x = 1
    while windows[1] == 0
        windows = (exp10.(range(log10(x), log10(10), bins))) # make log array
        windows = (windows / sum(windows)) * size(G)[1] # normalize it so that elements add up to correct size
        windows = floor.(Int, LowerTriangular(ones(Int, bins, bins)) * windows) # make valid indices
        x += 1
        #sanity check, this sum below should be almost equal to the uncompressed array size
        #sum(windows[2:end]-windows[1:end-1]) 
    end

    W0 = Diagonal(inv(LowerTriangular(ones(Int, bins, bins))) * windows)

    # Window average matrix, A

    A = zeros(bins, size(G)[1])
    for i in 1:bins
        if i == 1
            A[i, 1:windows[1]] .= 1 / diag(W0)[i]
        else
            A[i, (windows[i-1]+1):windows[i]] .= 1 / diag(W0)[i]
        end
    end


    t_direct = A * t_direct # Replace old direct time array with compressed one
    G = A * G # Replace old G with compressed one

    # sanity check plot:
    usv1 = svd(sqrt(W0) * K1) #paper (13)
    # surface(G, camera=(110, 25), xscale=:log10)

end
