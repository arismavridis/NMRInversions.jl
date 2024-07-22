struct invres2D
    ps::String
    dir::Vector
    indir::Vector
    f::Vector
    r::Vector
    SNR::Real
    alpha::Real
end


function invert(exptype::inversion2D,directory::String=pick_folder(pwd());kwargs...)

    invert(exptype, import_spinsolve(directory)... ;kwargs...)

end


function invert(
    exptype::inversion2D, t_direct::AbstractVector, t_indirect::AbstractVector, Raw::AbstractMatrix;
    α=:gcv, rdir=(-5, 1, 100), rindir=(-5, 1, 100),
    solver=brd, aopt=:none, order=0, savedata::Bool=true, plot::Bool=true)

    svds = svdcompress(exptype, t_direct, t_indirect, Raw, rdir=rdir, rindir=rindir)

    if isa(α, Real)
        f, r = solve_regularization(svds.K, svds.g, α, solver, order)

    elseif α == :gcv
        s̃ = svds.s
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
            φₙ, αₙ = gcv_score(αₙ, r, svds.s, (svds.V' * f)) # Compute φ for current α, and also compute new α 
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
            write(io, "SNR : " * string(svds.SNR) * "\n")
            write(io, "alpha : " * string(α) * "\n")
            write(io, "Direct Dimension : " * join(svds.x_dir, ", ") * "\n")
            write(io, "Indirect Dimension : " * join(svds.x_indir, ", ") * "\n")
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




