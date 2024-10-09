module tikhonov_jump_ext

using NMRInversions, JuMP, LinearAlgebra, SparseArrays

# import JuMP: @variable
# import JuMP: @variables
# import JuMP: @constraint
# import JuMP: @constraints
# import JuMP: @objective

function NMRInversions.solve_regularization(K::AbstractMatrix, g::AbstractVector,
        α::Real, solver::Type{jump_nnls})

    A = sparse([K; √(α) .* NMRInversions.Γ(size(K, 2), solver.order)])
    b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    @eval model = JuMP.Model($(solver.solver).Optimizer)
    set_silent(model)
    @variable(model, f[1:size(K, 2)] >= 0)
    @variable(model, z[1:size(b, 1)])
    @constraint(model, z .== A * f - b)
    @objective(model, Min, sum(z .^ 2))

    JuMP.optimize!(model)

    return JuMP.value.(f), vec(K * JuMP.value.(f) .- g)

end


# other qp formulations

#     # set_silent(model)
#     # @variable(model, f[1:size(K, 2)] >= 0)
#     # @variable(model, r[1:size(K, 1)])
#     # @constraint(model, r .== g .- K * f)
#     # @objective(model, Min, LinearAlgebra.dot(r, r) + α * LinearAlgebra.dot(f, f))
#     # optimize!(model)
#     # @assert is_solved_and_feasible(model)
#     # return value.(f), value.(r)


end
