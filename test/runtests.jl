using NMRInversions
using LinearAlgebra
using Test
using Optimization, OptimizationOptimJL


function test1D(exptype::Type{<:inversion1D})

    x = exp10.(range(log10(1e-4), log10(5), 32)) # acquisition range
    X = exp10.(range(-5, 1, 128)) # T range

    K = create_kernel(exptype, x, X)
    f_custom = [0.5exp.(-(x)^2 / 3) + exp.(-(x - 1.3)^2 / 0.5) for x in range(-5, 5, length(X))]

    g = K * f_custom
    y = g + 0.001 * maximum(g) .* randn(length(x))

    # results = invert(exptype, x, y, α=gcv)
    results = invert(exptype, x, y, α=gcv, solver=linear_tikhonov,order=1)

    return norm(results.f - f_custom) < 0.5
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

    results = invert(exptype, x_direct, x_indirect, y, α=gcv)

    return norm(resutls.f - f_custom) < 0.5

end

function test_phase_correction()

    # Create real and imaginary parts
    Re_original = exp.(-range(1, 20, 1000)) + randn(1000) .* 0.01
    Im_original = randn(1000) .* 0.01
    
    # Get them out of phase
    ϕd = rand() * 2π
    Re_shifted, Im_shifted = NMRInversions.phase_shift(Re_original, Im_original, ϕd)
    
    # Correct the phase
    Rₙ, Iₙ, ϕc = NMRInversions.autophase(Re_shifted, Im_shifted, 1)
    
    # Plots for sanity check
    
    # p1 = plot([Re_original, Im_original], label=["Original real" "Original Imaginary"])
    # p2 = plot([Re_shifted, Im_shifted], label=["Dephased real" "Dephased Imaginary"])
    # p3 = plot([Rₙ, Iₙ], label=["Corrected real" "Corrected Imaginary"])
    # ϕ_range = range(0, 2π, 20000)
    # Re1_vs_φ = Re_shifted[1] .* cos.(ϕ_range) - Im_shifted[1] .* sin.(ϕ_range)
    # Im_sum_vs_φ = [im_cost([ϕ], (Re_shifted, Im_shifted)) for ϕ in ϕ_range]
    # p4 = plot(ϕ_range, Re1_vs_φ, xlabel="ϕ", label="Re[1]")
    # p4 = plot!(ϕ_range, (Im_sum_vs_φ ./ maximum(Im_sum_vs_φ)) .* maximum(Re1_vs_φ), xlabel="ϕ", label="sum(im.^2)", legend=:topleft)
    # p4 = vline!(p4, [ϕc], label="corrected phase")
    # display(plot(p1, p2, p3, p4))

    display("The correction error is $(2π - (ϕd + ϕc)) radians")

    return abs(2π - (ϕd + ϕc)) < 0.01
end



@testset "NMRInversions.jl" begin
    # Write your tests here.
    @test test1D(IR)
    @test test1D(CPMG)
    @test test_phase_correction()

end

