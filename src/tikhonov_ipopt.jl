using JuMP

struct ip_tkh end
IP = ip_tkh()

function solve_tikhonov(K::AbstractMatrix,g::AbstractVector,α ; solver::ip_tkh = IP,order=0)
end
