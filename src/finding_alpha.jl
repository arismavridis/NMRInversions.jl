"""
Compute the Generalized Cross Validation (GCV) score for a regularization solution, and return the score, as well as the next α to test. 
(look Mitchell 2012 paper https://doi.org/10.1016/j.pnmrs.2011.07.002)

- `α` : smoothing term
- `r` : residuals from the regularization (`r = Kf .- g`)
- `s` : singular values of the kernel
- `x` : `Ṽ₀'f` or `s .* c`

The symbols above follow the notation used in the paper.

"""
function gcv_score(α, r, s, x; next_alpha=true)

    ñ = length(s)
    c = r ./ α
    σ² = s .^ 2
    m̂ = sum(σ² ./ (σ² .+ α))
    φₙ = (α .^ 2 * c' * c .* ñ) ./ ((ñ - m̂) .^ 2)  # GCV score to be minimized

    if next_alpha

        dm̂ = sum(σ² ./ ((σ² .+ α) .^ 2))
        t̂ = sum(x .^ 2 ./ (σ² .+ α))
        αₙ = (α .^ 2 * c' * c * dm̂) / (t̂ * (ñ - m̂))  # α update value, test this one next
        return φₙ, αₙ

    else
        return φₙ
    end
end


"""
Solve repeatedly until the GCV score stops decreasing.
Select the solution with minimum gcv score and return it, along with the residuals.
"""
function solve_gcv(svds::svd_kernel_struct, solver::Type{<:regularization_solver}, order::Int)

    s̃ = svds.S
    ñ = length(s̃)

    #Initial guess (overestimate)
    αₙ = sum(s̃ .^ 2) / ñ
    α = []
    φ = []
    f_star = []

    done = false
    while ~done

        display("Testing α = $(round(αₙ,sigdigits=3))")
        f, r = solve_regularization(svds.K, svds.g, αₙ, solver, order)

        push!(α, αₙ) # Add the just tested α to the array
        φₙ, αₙ = gcv_score(αₙ, r, svds.S, (svds.V' * f)) # Compute φ for current α, and also compute new α 
        push!(φ, φₙ)

        if length(φ) > 1 && φ[end] < φ[end-1]
            f_star = deepcopy(f)

        elseif length(φ) > 1 && φ[end] >= φ[end-1]
            done = true

            display("The optimal α is $(round(α[end-1],digits=3))")
        end

    end

    α = α[end-1]
    f = f_star
    r = svds.g - svds.K * f

    return f, r, α
end

"""
Compute the curvature of the L-curve at a given point.

- `f` : solution vector
- `r` : residuals
- `α` : smoothing term
- `A` : Augmented kernel matrix (`K` and `αI` stacked vertically)

"""
function l_curvature(f, r, α, A)

    ξ = f'f
    ρ = r'r
    λ = √α

    z = solve_nnls(A, r)

    ∂ξ∂λ = 4 / λ * f'z

    ĉ = (2ξ * ρ / ∂ξ∂λ) * (α * ∂ξ∂λ * ρ + 2 * ξ * λ * ρ + λ^4 * ξ * ∂ξ∂λ) / (α * ξ^2 + ρ^2)^(3 / 2)

    return ĉ

end


"""

"""
function solve_l_curve(K, g, lower, upper, n, order)

    alphas = exp10.(range(log10(lower), log10(upper), n))
    curvatures = zeros(length(alphas))

    ξarray = zeros(length(alphas))
    ρarray = zeros(length(alphas))

    for (i, α) in enumerate(alphas)
        A = sparse([K; √(α) .* NMRInversions.Γ(size(K, 2), order)])

        f = vec(nonneg_lsq(A, [y; zeros(size(A, 1) - size(y, 1))], alg=:nnls))
        r = K * f - y

        ξ = f'f
        ρ = r'r

        λ = √α

        ξarray[i] = ξ
        ρarray[i] = ρ

        z = vec(nonneg_lsq(A, [r; zeros(size(A, 1) - size(r, 1))], alg=:nnls))

        ∂ξ∂λ = (4 / λ) * f'z

        ĉ = (2ξ * ρ / ∂ξ∂λ) * (α * ∂ξ∂λ * ρ + 2 * ξ * λ * ρ + λ^4 * ξ * ∂ξ∂λ) / (α * ξ^2 + ρ^2)^(3 / 2)

        curvatures[i] = ĉ

    end
    # plot(ρarray, ξarray, xscale=:log10, yscale=:log10)

    non_inf_indx = findall(!isinf, curvatures)
    argmax(curvatures[non_inf_indx])
    α = alphas[non_inf_indx][argmax(curvatures[non_inf_indx])]
    A = sparse([K; √(α) .* NMRInversions.Γ(size(K, 2), order)])
    f = vec(nonneg_lsq(A, [y; zeros(size(A, 1) - size(y, 1))], alg=:nnls))
    r = K * f - y

    return f, r, α

end
