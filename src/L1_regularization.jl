
function PDHGM(K, s, α; tol=1e-6, τ=0.1, σ=0.1)

    B = inv(I + τ * α * K' * K)
    Y = zeros(size(K, 2))
    Ỹ = deepcopy(Y)
    f = ones(size(K, 2))
    f̃ = deepcopy(f)
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

f = PDHGM(rand(10, 1000), rand(10), 100)