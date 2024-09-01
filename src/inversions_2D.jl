
function invert(input::input2D; kwargs...)

    invert(input.seq, input.x_direct, input.x_indirect, input.data; kwargs...)

end


function invert(
    seq::Type{<:pulse_sequence2D}, x_direct::AbstractVector, x_indirect::AbstractVector, Data::AbstractMatrix;
    α=gcv, rdir=(-5, 1, 100), rindir=(-5, 1, 100),
    solver=brd, order=0, savedata::Bool=true)

    X_direct = exp10.(range(rdir...))
    X_indirect = exp10.(range(rindir...))

    ker_struct = create_kernel(seq, x_direct, x_indirect, X_direct, X_indirect, Data)

    if isa(α, Real)
        f, r = solve_regularization(ker_struct.K, ker_struct.g, α, solver, order)

    elseif α == gcv
        f, r, α = solve_gcv(ker_struct, solver, order)

    else
        error("α must be a real number or gcv")

    end

    F = reshape(f, length(X_direct), length(X_indirect))

    results = invres2D(seq, X_direct, X_indirect, F, r, calc_snr(Data), α,[],[])

    if savedata == true
        writeresults("inversion_results.txt", results)
        display("Data saved as inversion_results.txt in $(pwd())")
    end

    return results

end
