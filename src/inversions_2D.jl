struct invres2D
    ps::String
    dir::Vector
    indir::Vector
    f::Vector
    r::Vector
    SNR::Real
    alpha::Real
end


function invert(exptype::Type{<:inversion2D}, directory::String=pick_folder(pwd()); kwargs...)

    cd(directory)
    invert(exptype, import_spinsolve(pwd())...; kwargs...)
end


function invert(
    exptype::Type{<:inversion2D}, x_direct::AbstractVector, x_indirect::AbstractVector, Raw::AbstractMatrix;
    α=:gcv, rdir=(-5, 1, 100), rindir=(-5, 1, 100),
    solver=song,  order=0, savedata::Bool=true, plot::Bool=true)

    X_direct = exp10.(range(rdir...))
    X_indirect = exp10.(range(rindir...))

    svds = create_kernel(exptype, x_direct, x_indirect, X_direct, X_indirect, Raw)

    if isa(α, Real)
        f, r = solve_regularization(svds.K, svds.g, α, solver, order)

    elseif α == :gcv
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

    end


    #Save .txt file
    if savedata == true
        open("inversion_results.txt", "w") do io
            write(io, "Pulse Sequence : " * "IR-CPMG" * "\n")
            write(io, "SNR : " * string(calc_snr(Raw)) * "\n")
            write(io, "alpha : " * string(α) * "\n")
            write(io, "Direct Dimension : " * join(X_direct, ", ") * "\n")
            write(io, "Indirect Dimension : " * join(X_indirect, ", ") * "\n")
            write(io, "Inversion Results : " * join(f, ", ") * "\n")
            write(io, "Residuals : " * join(r, ", ") * "\n")

        end
        display("Data saved as inversion_results.txt")
    end

    if plot == true
        return select_peaks(joinpath(pwd(), "inversion_results.txt"))
    end

    return f, r

end




