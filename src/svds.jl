struct svdstruct
    K::AbstractMatrix
    g::AbstractVector
    s::AbstractVector
    V::AbstractMatrix
    x_dir::AbstractVector
    x_indir::AbstractVector
    sigma_n::AbstractFloat
    SNR::AbstractFloat
end



# Compression of 2D data
function svdcompress(exptype::inversion2D, t_direct::AbstractVector, t_indirect::AbstractVector, Raw::AbstractMatrix;
    rdir=(-5, 1, 100), rindir=(-5, 1, 100), reducesize=false) 
    G = real(Raw)
    ## Determine SNR
    # σ_n is the STD of the latter half of the imaginary signal (former half might contain signal residues as well) 
    # sigma_n = mean(std((imag(Raw)[(round(Int64, (size(imag(Raw),1) / 2))):end, :]),dims=2))
    sigma_n = std((imag(Raw)[(floor(Int64, (size(imag(Raw), 1) / 2))):end, :]))
    SNR = maximum(abs.(G)) / sigma_n

    if SNR < 1000
        @warn("The SNR is $(round(SNR, digits=1)), which is below the recommended value of 1000. Consider running experiment with more scans.")
    end

    # Kernel ranges
    x_direct = exp10.(range(rdir...)) # Range of direct dimension 
    x_indirect = exp10.(range(rindir...)) # Range of indirect dimension


    # Kernel functions
    # K1_func(t, T2) = exp(-t / T2) # Kernel equation
    # K2_func(t, T1) = 1 - 2 * exp(-t / T1) # Kernel equation

    if exptype == IRCPMG 
        K1_func = (t, T2) -> exp(-t / T2)
        K2_func = (t, T1) -> 1 - 2 * exp(-t / T1)
    
    end

    # Generate Kernels
    K1 = K1_func.(t_direct, x_direct')
    K2 = K2_func.(t_indirect, x_indirect')


    ## Data compression

    if reducesize == true
        # Compress direct dimension to the following length
        bins = 64

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
        # surface(G, camera=(110, 25), xscale=:log10)

    end

    ## Perform SVD truncation
    if reducesize == true
        usv1 = svd(sqrt(W0) * K1) #paper (13)
    else
        usv1 = svd(K1) #paper (13)
    end
    usv2 = svd(K2) #paper (14)

    # finding which singular components are contributed from K1 and K2
    S21 = usv2.S * usv1.S' #Using outer product instead of Kronecker, to make indices more obvious
    indices = findall(i -> i .> (sigma_n / findmax(G)[1]), S21) # find elements in S12 above the noise threshold
    s̃ = S21[indices]
    ñ = length(s̃)

    si = (first.(Tuple.(indices)))  # direct dimension singular vector indices
    sj = (last.(Tuple.(indices)))   # indirect dimension singular vector indices

    g̃ = diag(usv2.U[:, si]' * G' * usv1.U[:, sj])

    V1t = repeat(usv1.V[:, sj], size(usv2.V, 1), 1)
    V2t = reshape(repeat(usv2.V[:, si]', size(usv1.V, 1), 1), ñ, size(V1t, 1))'
    Ṽ₀ = V1t .* V2t
    K̃₀ = Diagonal(s̃) * Ṽ₀'


    return svdstruct(K̃₀, g̃, s̃, Ṽ₀, x_direct, x_indirect, sigma_n, SNR)
end



function svdcompress(exptype::String, t::AbstractVector, g::AbstractVector; rdir=(-5, 1, 100))

    if exptype == "T1"
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)

    elseif exptype == "T2"
        kernel_eq = (t, T) -> exp(-t / T)

    elseif exptype == "D"
        kernel_eq = (t, D) -> exp(-t / D)

    end

    T_range = exp10.(range(rdir...))
    kernel = kernel_eq.(t, T_range')

    usv = svd(kernel)

    p1 = scatter([usv.S, abs.(usv.U' * g), abs.(usv.U' * g ./ usv.S)], yscale=:log10, label=["σ" "|U'g|" "|U'g|./σ"])

    display(p1)

    display("How many singular values do you want to keep? ")

    i = parse(Int, readline())
    # i = 15

    # Truncate singular values after i 
    s̃ = usv.S[1:i]
    Ṽ = usv.V[:, 1:i]
    K̃ = Diagonal(s̃) * Ṽ'
    g̃ = usv.U[:, 1:i]' * g

    return svdstruct(K̃, g̃, s̃, Ṽ, T_range, [], 0, 0)

end