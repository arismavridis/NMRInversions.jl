using LinearAlgebra

"""
A. Reci et al. / Journal of Magnetic Resonance 281 (2017) 188–198

The particular choice of σ and r is heuristic. 
A smaller σ will increase the stability while reducing the convergence speed
of the algorithm. A good compromise between the two was found when σ = 0.1 and r = 10. 
The best values of σ and r will depend slightly on the scaling of the signal. 
To avoid this, it is best to normalize the NMR signal t a maximum of 1, a technique which was followed in this study.
"""
function PDHGM(K::AbstractMatrix, s::AbstractVector, α::Real; tol=10, τ=10 , σ=0.1)

    K = K ./ maximum(K)
    s = s ./ maximum(s)
    B = inv(LinearAlgebra.I + τ * α * K' * K)
    Y = zeros(size(K, 2))
    Ỹ = copy(Y)
    f = ones(size(K, 2))
    f̃ = copy(f)
    f_prev = deepcopy(f)
    ε = tol + 1

    while ε > tol
        Ỹ .= Y + σ * f̃
        Y .= Ỹ ./ max.(1, abs.(Ỹ))
        f .= B * (f̃ - τ .* Y + τ * α * K' * s)
        f .= max.(0, f)
        f̃ .= 2 .* f .- f_prev
        ε = norm(f - f_prev) / norm(f_prev)
        f_prev .= f
    end
    return f
end


f = PDHGM(K, y, 1, tol=10, τ=10, σ=0.1); plot(f)