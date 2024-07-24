"""
Create a finite difference matrix of order `order` for a matrix of size `m x m`.
"""
function Γ(m::Int, order::Int)
    # Eilers, P. H. C. (2003). Analytical Chemistry, 75(14), 3631–3636. (Supporting Information)
    # Taken from RegularizationTools.jl
    if order == 0
        return Array{Float64}(LinearAlgebra.I, (m, m))
    end
    return diff(Γ(m, order - 1), dims=1)
end


"""
Compute the Generalized Cross Validation (GCV) score, 
and return the score, as well as the next α to test. 
In the inputs, r is the residuals of the solution, 
s are the singular values of the kernel and x=Ṽ₀'f or x= s .* c 
(look Mitchell 2012 paper https://doi.org/10.1016/j.pnmrs.2011.07.002)
"""
function gcv_score(α, r, s, x)

    ñ = length(s)
    c = r ./ α
    σ² = s .^ 2
    m̂ = sum(σ² ./ (σ² .+ α))
    dm̂ = sum(σ² ./ ((σ² .+ α) .^ 2))
    t̂ = sum(x .^ 2 ./ (σ² .+ α))
    φₙ = (α .^ 2 * c' * c .* ñ) ./ ((ñ - m̂) .^ 2)  # GCV score to be minimized
    αₙ = (α .^ 2 * c' * c * dm̂) / (t̂ * (ñ - m̂))  # α update value, test this one next

    return φₙ, αₙ
end


"""
Solve repeatedly until the GCV score stops decreasing.
Select the solution with minimum score and return it, along with the residuals.
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

        display("Testing α = $(round(αₙ,digits=3))")
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

    return f, r
end