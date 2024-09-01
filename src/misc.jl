"""
Create a finite difference differentiation matrix of order `order` for a matrix of size `m x m`.
"""
function Γ(m::Int, order::Int)
    # Eilers, P. H. C. (2003). Analytical Chemistry, 75(14), 3631–3636. (Supporting Information)
    # Taken from RegularizationTools.jl
    if order == 0
        return Array{Float64}(LinearAlgebra.I, (m, m))
    end
    return diff(Γ(m, order - 1), dims=1)
end


function selections(inv_results::invres2D)

    dir = inv_results.X_dir
    indir = inv_results.X_indir
    F = inv_results.F

    x = collect(1:length(indir))
    y = collect(1:length(dir))
    z = copy(F)

    points = [[i, j] for i in x, j in y] #important, used by inpolygon later
    mask = zeros(size(points))

    polygons = inv_results.sp
    selections = [zeros(size(F)) for _ in 1:length(polygons)]

    for (i, polygon) in enumerate(polygons)
        mask .= [PolygonOps.inpolygon(p, polygon; in=1, on=1, out=0) for p in points]
        selections[i] = mask .* z 
    end

    return selections

end

function calc_weighted_averages(F::Matrix, dir::Vector, indir::Vector)

    T1 = vec(sum(F, dims=2)) ⋅ indir / sum(F)
    T2 = vec(sum(F, dims=1)) ⋅ dir / sum(F)

    return [T1, T2]
end
