
function solve_tikhonov(K::AbstractMatrix, g::AbstractVector, α::Real, solver::Symbol, order::Int=0)

    A = sparse([K; √(α) .* Γ(size(K, 2), order)])
    b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    @eval model = Model($(solver).Optimizer)
    # set_silent(model)
    @variable(model, f[1:size(K,2)] >= 0)
    @variable(model, z[1:size(b,1)])
    @constraint(model, z .== A*f - b)
    @objective(model, Min, sum(z .^ 2))

    optimize!(model)

    return value.(f), vec(K * value.(f) .- g)


end


