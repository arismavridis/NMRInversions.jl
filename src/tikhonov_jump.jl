
import JuMP: @variable
import JuMP: @variables
import JuMP: @constraint
import JuMP: @constraints
import JuMP: @objective

function solve_tikhonov(K::AbstractMatrix, g::AbstractVector, α::Real, solver::Symbol, order::Int=0)

    A = sparse([K; √(α) .* Γ(size(K, 2), order)])
    b = sparse([g; zeros(size(A, 1) - size(g, 1))])

    @eval model = JuMP.Model($(solver).Optimizer)
    # set_silent(model)
    @variable(model, f[1:size(K, 2)] >= 0)
    @variable(model, z[1:size(b, 1)])
    @constraint(model, z .== A * f - b)
    @objective(model, Min, sum(z .^ 2))

    JuMP.optimize!(model)

    return JuMP.value.(f), vec(K * JuMP.value.(f) .- g)


end


# other qp formulations

#     # if solver == :Gurobi
#     #     model = Model(() -> Gurobi.Optimizer(GRB_ENV))
#     # elseif solver == :HiGHS
#     #     model = Model(HiGHS.Optimizer)
#     # else
#     #     error("Unknown solver")
#     # end

#     # set_silent(model)
#     # @variable(model, f[1:size(K, 2)] >= 0)
#     # @variable(model, r[1:size(K, 1)])
#     # @constraint(model, r .== g .- K * f)
#     # @objective(model, Min, LinearAlgebra.dot(r, r) + α * LinearAlgebra.dot(f, f))
#     # optimize!(model)
#     # @assert is_solved_and_feasible(model)
#     # return value.(f), value.(r)


function solve_tikhonov(K::AbstractMatrix, g::AbstractVector, α::Real, solver::jump_L1_solver, order::Int=0)

    m, n = size(K)
    model = JuMP.Model(Ipopt.Optimizer)

    ## Method 1
    @variable(model, f[1:n] >= 0)
<<<<<<< HEAD
    @objective(model, Min, sum((K * f .- g) .^ 2) + a * sum(abs.(f)))
=======
    @objective(model, Min, sum((K * f .- s) .^ 2) + a * sum(abs.(f)))
>>>>>>> d91c22f436d08919efc57d216a3c464da74f9d3b

    ## Method 2
    # @variables(model, begin
    #     f[1:n] >= 0
    #     residuals[1:m]
    #     l1_terms[1:n] >= 0
    # end)
    # @constraints(model, begin
    #     residuals .== K * f .- g
    #     l1_terms .>= f
    #     l1_terms .>= -f
    # end)
    # @objective(model, Min, sum(residuals .^ 2) + α * sum(l1_terms))


    JuMP.optimize!(model)
    # @assert is_solved_and_feasible(model)
    JuMP.value.(f)

    return JuMP.value.(f), vec(K * JuMP.value.(f) .- g)

end
