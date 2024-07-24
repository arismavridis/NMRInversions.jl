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
    α=gcv, rdir=(-5, 1, 100), rindir=(-5, 1, 100),
    solver=song,  order=0, savedata::Bool=true, plot::Bool=true)

    X_direct = exp10.(range(rdir...))
    X_indirect = exp10.(range(rindir...))

    ker_struct = create_kernel(exptype, x_direct, x_indirect, X_direct, X_indirect, Raw)

    if isa(α, Real)
        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver, order)

    elseif α == gcv
        f, r , α= solve_gcv(ker_struct, solver, order)

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




