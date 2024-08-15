using NMRInversions
using LinearAlgebra
using Test


function test1D(exptype::Type{<:inversion1D})

    x = exp10.(range(log10(1e-4), log10(5), 32)) # acquisition range
    X = exp10.(range(-5, 1, 128)) # T range

    K = create_kernel(exptype, x, X)
    f_custom = [0.5exp.(-(x)^2 / 3) + exp.(-(x - 1.3)^2 / 0.5) for x in range(-5, 5, length(X))]

    g = K * f_custom
    y = g + 0.001 * maximum(g) .* randn(length(x))

    f_estimated, r = invert(exptype, x, y, α=gcv)

    return norm(f_estimated - f_custom) < 0.5
end

function test2D(exptype::Type{<:inversion2D})

    x_direct = exp10.(range(log10(1e-4), log10(5), 1024)) # acquisition range
    x_indirect = exp10.(range(log10(1e-4), log10(5), 32)) # acquisition range

    X_direct = exp10.(range(-5, 1, 64)) # T range
    X_indirect = exp10.(range(-5, 1, 64)) # T range

    θ = 135
    σ₁ = 1.3
    σ₂ = 0.4
    x₀ = 0
    y₀ = 1.3
    a = ((cosd(θ)^2) / (2 * σ₁^2)) + ((sind(θ)^2) / (2 * σ₂^2))
    b = -((sind(2 * θ)) / (4 * σ₁^2)) + ((sind(2 * θ)) / (4 * σ₂^2))
    c = ((sind(θ)^2) / (2 * σ₁^2)) + ((cosd(θ)^2) / (2 * σ₂^2))
    f_custom = ([exp.(-(a * (x - x₀)^2 + 2 * b * (x - x₀) * (y - y₀) + c * (y - y₀)^2)) for x in range(-5, 5, 100), y in range(-5, 5, 100)])
    K = create_kernel(exptype, x_direct, x_indirect, X_direct, X_indirect, complex.(f_custom, 0.001 * randn(size(f_custom)) .* f_custom)) 

    g = K * f_custom
    y = g + 0.001 * maximum(g) .* randn(length(x_direct) * length(x_indirect))

    f_estimated, r = invert(exptype, x_direct, x_indirect, y, α=gcv)

    return norm(f_estimated - f_custom) < 0.5

end


begin
end

@testset "NMRInversions.jl" begin
    # Write your tests here.
    @test test1D(IR)
    @test test1D(CPMG)

end

