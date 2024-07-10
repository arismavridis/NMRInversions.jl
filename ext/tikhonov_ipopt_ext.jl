module tikhonov_ipopt_ext

using NMRInversions, Tulip


function NMRInversions.solve_tikhonov(K::AbstractMatrix, g::AbstractVector, α::Real, solver::NMRInversions.ip_solver, order::Int=0)
    display("Solving Tikhonov with Ipopt")
end

export solve_tikhonov

end