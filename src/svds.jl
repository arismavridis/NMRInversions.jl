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

    g = real(raw)

    ## determine snr
    # σ_n is the std of the latter half of the imaginary signal (former half might contain signal residues as well) 
    # sigma_n = mean(std((imag(raw)[(round(int64, (size(imag(raw),1) / 2))):end, :]),dims=2))
    # sigma_n = std((imag(raw)[(floor(int64, (size(imag(raw), 1) / 2))):end, :]))

    n = (imag(raw)[(floor(int64, (size(imag(raw), 1) / 2))):end, :]) # noisy part of the signal
    sigma_n = sqrt(sum((n.- mean(n)).^2) / (length(n) - 1)) # std of the noisy part 
    snr = maximum(abs.(g)) / sigma_n

    if snr < 1000
        @warn("the snr is $(round(snr, digits=1)), which is below the recommended value of 1000. consider running experiment with more scans.")
    end

    # kernel ranges
    x_direct = exp10.(range(rdir...)) # range of direct dimension 
    x_indirect = exp10.(range(rindir...)) # range of indirect dimension


    # kernel functions

    if exptype == ircpmg
        k1_func = (t, t2) -> exp(-t / t2)
        k2_func = (t, t1) -> 1 - 2 * exp(-t / t1)

    end

    # generate kernels
    k1 = k1_func.(t_direct, x_direct')
    k2 = k2_func.(t_indirect, x_indirect')


    ## data compression

    if reducesize == true
        # compress direct dimension to the following length
        bins = 64

        # use logarithmically spaced windows for a window average
        windows = zeros(bins)

        x = 1
        while windows[1] == 0
            windows = (exp10.(range(log10(x), log10(10), bins))) # make log array
            windows = (windows / sum(windows)) * size(g)[1] # normalize it so that elements add up to correct size
            windows = floor.(int, lowertriangular(ones(int, bins, bins)) * windows) # make valid indices
            x += 1
            #sanity check, this sum below should be almost equal to the uncompressed array size
            #sum(windows[2:end]-windows[1:end-1]) 
        end

        w0 = diagonal(inv(lowertriangular(ones(int, bins, bins))) * windows)

        # window average matrix, a

        a = zeros(bins, size(g)[1])
        for i in 1:bins
            if i == 1
                a[i, 1:windows[1]] .= 1 / diag(w0)[i]
            else
                a[i, (windows[i-1]+1):windows[i]] .= 1 / diag(w0)[i]
            end
        end


        t_direct = a * t_direct # replace old direct time array with compressed one
        g = a * g # replace old g with compressed one

        # sanity check plot:
        # surface(g, camera=(110, 25), xscale=:log10)

    end

    ## perform svd truncation
    if reducesize == true
        usv1 = svd(sqrt(w0) * k1) #paper (13)
    else
        usv1 = svd(k1) #paper (13)
    end
    usv2 = svd(k2) #paper (14)

    # finding which singular components are contributed from k1 and k2
    s21 = usv2.s * usv1.s' #using outer product instead of kronecker, to make indices more obvious
    indices = findall(i -> i .> (sigma_n / findmax(g)[1]), s21) # find elements in s12 above the noise threshold
    s̃ = s21[indices]
    ñ = length(s̃)

    si = (first.(tuple.(indices)))  # direct dimension singular vector indices
    sj = (last.(tuple.(indices)))   # indirect dimension singular vector indices

    g̃ = diag(usv2.u[:, si]' * g' * usv1.u[:, sj])

    v1t = repeat(usv1.v[:, sj], size(usv2.v, 1), 1)
    v2t = reshape(repeat(usv2.v[:, si]', size(usv1.v, 1), 1), ñ, size(v1t, 1))'
    ṽ₀ = v1t .* v2t
    k̃₀ = diagonal(s̃) * ṽ₀'


    return svdstruct(k̃₀, g̃, s̃, ṽ₀, x_direct, x_indirect, sigma_n, snr)
end



function svdcompress(exptype::inversion1D, t::AbstractVector, g::AbstractVector; rdir=(-5, 1, 100))

    if exptype == IR
        kernel_eq = (t, T) -> 1 - 2 * exp(-t / T)

    elseif exptype == CPMG
        kernel_eq = (t, T) -> exp(-t / T)

    elseif exptype == PFG
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