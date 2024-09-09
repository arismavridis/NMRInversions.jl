module gui_ext

using NMRInversions, GLMakie, PolygonOps, LinearAlgebra, NativeFileDialog

# TESTING
# res = invert(import_spinsolve(), alpha=0.7)
# plot(res)

function Makie.plot(results::NMRInversions.inv_out_2D; kwargs...)

    f = Figure(size=(500, 500))
    plot!(f, results; kwargs...)

    return f

end

function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::NMRInversions.inv_out_2D;
    title="", clmap=:tempo)

    x = res.X_indir
    y = res.X_dir

    xplot = collect(range(0, 1, length(x)))
    yplot = collect(range(0, 1, length(y)))

    # Make axes

    axmain = Axis(fig[3:10, 1:8])
    axtop = Axis(fig[1:2, 1:8], limits=((xplot[1], xplot[end]), (0, nothing)), title=title, titlesize=17)
    axright = Axis(fig[3:10, 9:10], limits=((0, nothing), (yplot[1], yplot[end])))

    hidedecorations!(axtop)
    hidedecorations!(axright)
    hidedecorations!(axmain)

    axmain_values = Axis(fig[3:10, 1:8],
        xscale=log10, yscale=log10,
        xlabel=L"T_1 \, \textrm{(s)}", ylabel=L"T_2 \,\textrm{(s)}",
        xlabelsize=23, ylabelsize=23,
        limits=(x[1], x[end], y[1], y[end])
    )

    Makie.deactivate_interaction!(axmain_values, :rectanglezoom) # Disable zoom

    linkxaxes!(axmain, axtop)
    linkyaxes!(axmain, axright)

    draw_on_axes(axmain, axtop, axright, res, clmap)

end

function draw_on_axes(axmain, axtop, axright, res, clmap)

    empty!(axmain)
    empty!(axtop)
    empty!(axright)

    z = res.F' .* res.filter'
    x = range(0, 1, size(z, 1))
    y = range(0, 1, size(z, 2))

    # Plots
    contourf!(axmain, x, y, z, colormap=clmap, levels=50)
    lines!(axtop, x, vec(sum(z, dims=2)), colormap=clmap, colorrange=(1, 10), color=5, alpha=0.5)
    lines!(axright, vec(sum(z, dims=1)), y, colormap=clmap, colorrange=(1, 10), color=5, alpha=0.5)

    # Plot diagonal line
    lines!(axmain, [(0, 0), (1, 1)], color=:black, linewidth=1)

    #Create a matrix for all the discrete points in the space
    points = [[i, j] for i in x, j in y]
    mask = zeros(size(points))

    for (i, polygon) in enumerate(res.selections)

        # Selected points only
        mask .= [PolygonOps.inpolygon(p, polygon; in=1, on=1, out=0) for p in points]
        spo = mask .* z

        indir_dist = vec(sum(spo, dims=2))
        dir_dist = vec(sum(spo, dims=1))

        #draw current polygon
        lines!(axmain, Point2f.(polygon), linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.9)

        T1 = indir_dist ⋅ res.X_indir / sum(spo)
        T2 = dir_dist ⋅ res.X_dir / sum(spo)

        # Plot peak label
        xc = vec(sum(spo, dims=2)) ⋅ x / sum(spo)
        yc = vec(sum(spo, dims=1)) ⋅ y / sum(spo)

        sc = scatter!(axmain, xc, yc, markersize=15,
            marker=collect('a':'z')[i],
            colormap=:tab10, colorrange=(1, 10),
            color=i,
            glowcolor=:white, glowwidth=4,
            label=" : T₁/T₂ = $(round(T1/T2 , digits=1)) \n    Volume = $(round(sum(spo)/sum(z) *100 ,digits = 1))%")


        text!(axmain, 2.5 * (1 / 100), (99 - 13 * (i - 1)) * (1 / 100),
            text=
            collect('a':'z')[i] *
            " : T₁/T₂ = $(round(T1/T2 , digits=1)) \n     Volume = $(round(sum(spo)/sum(z) *100 ,digits = 1))%",
            align=(:left, :top), colormap=:tab10, colorrange=(1, 10), color=i)

        # lines!(axtop, indir_dist[], linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.5)
        # lines!(axright, dir_dist[], 1:length(res.X_dir), linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.5)
        # axislegend( axmain,  framevisible=false)

    end

end

function dynamic_plots(axmain, axtop, axright, selection, polygon)

    scatter!(axmain, selection, colormap=:tab10, colorrange=(1, 10), color=8)
    lines!(axmain, polygon, colormap=:tab10, colorrange=(1, 10), color=8)

end

function draw_3d(ax3d, res, sp)

    surface!(ax3d, res.F' .* res.filter', colormap=:tempo)

    surface!(ax3d, spo, colormap=:blues,
        colorrange=(minimum(filter(x -> x != 0, vec(z))), maximum(z)),
        lowclip=:transparent)

end

function NMRInversions.select_peaks(res::NMRInversions.inv_out_2D)
    # begin
    gui = Figure(size=(900, 500))
    Makie.plot!(gui[2:10, 1:9], res)

    axmain = gui.content[1]
    axtop = gui.content[2]
    axright = gui.content[3]
    # Legend(gui[5, 10:19], axmain)

    # Buttons
    labelb = Button(gui[1, 10:15]; label="Label current selection")
    clearb = Button(gui[2, 10:15]; label="Clear current selection")
    resetb = Button(gui[1, 15:19]; label="Reset everything")
    filterb = Button(gui[2, 15:19]; label="Filter-out unselected")
    saveb = Button(gui[3, 15:19]; label="Save plot (WIP)")

    # Title textbox
    tb = Textbox(gui[1, 2:7], placeholder="Insert a title for the plot, then press enter.", width=300, reset_on_defocus=true)

    z = res.F' .* res.filter'
    x = collect(range(0, 1, size(z, 1)))
    y = collect(range(0, 1, size(z, 2)))

    #Create a matrix for all the discrete points in the space
    points = [[i, j] for i in x, j in y] #important, used by inpolygon later
    mask = Observable(zeros(size(points)))

    # Selected points only
    spo = @lift(($mask .* $z))

    # Selected points Direct dimension distribution
    dir_dist = @lift(vec(sum($spo, dims=1)))
    indir_dist = @lift(vec(sum($spo, dims=2)))

    # Initialise Selection vectors
    selection = Observable(Point{2,Float32}[])
    polygon = Observable(Point{2,Float32}[])

    # Dynamic plots
    dynamic_plots(axmain, axtop, axright, selection, polygon)


    begin ## SELECTING POINTS

        # Pick points by clicking the plot
        spoint = select_point(gui.content[1].scene, marker=:circle)
        # Update selection Observable by pushing the selected point to it
        on(spoint) do _
            push!(selection[], spoint[])
            selection[] = selection[] # Update explicitly so that the plots are updated automatically

            if size(selection[], 1) > 2
                polygon[] = vcat(selection[], [selection[][1]])
                mask[] = [PolygonOps.inpolygon(p, polygon[]; in=1, on=1, out=0) for p in points]
            end
        end

    end ## SELECTING POINTS


    begin ## BUTTON CLICKS

        on(labelb.clicks) do _

            if size(selection[], 1) < 3
                @warn("You need to make a selection.")
            else

                push!(res.selections, polygon[])

                draw_on_axes(axmain, axtop, axright, res, :tempo)
                dynamic_plots(axmain, axtop, axright, selection, polygon)

                # Clear for next selection
                selection[] = Point{2,Float32}[]
                polygon[] = Point{2,Float32}[]
                mask[] = zeros(size(points))

                reset_limits!(axmain)
            end

        end

        on(clearb.clicks) do _

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))

        end

        on(filterb.clicks) do _

            if size(selection[], 1) < 3
                @warn("You need to make a selection.")
            else

                res.filter .= res.filter .* [PolygonOps.inpolygon(p, polygon[]; in=1, on=1, out=0) for p in points]'
            end

            draw_on_axes(axmain, axtop, axright, res, :tempo)
            dynamic_plots(axmain, axtop, axright, selection, polygon)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))

        end


        on(resetb.clicks) do _
            res.filter .= ones(size(res.F))
            empty!(res.selections)

            draw_on_axes(axmain, axtop, axright, res, :tempo)
            dynamic_plots(axmain, axtop, axright, selection, polygon)

            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))
        end

        on(saveb.clicks) do _

            ttl = ""
            if !isnothing(tb.stored_string[])
                ttl = tb.stored_string[]
            end

            f = plot(res, title=ttl)

        end

    end ## BUTTON CLICKS

    gui
end

end


function Makie.plot(res::NMRInversions.inv_out_1D; kwargs...)

    f = Figure(size=(500, 500))
    plot!(f, res; kwargs...)

    return f

end

function Makie.plot!(fig::Union{Makie.Figure,Makie.GridPosition}, res::NMRInversions.inv_out_1D;
    title="")

    # Make axes
    if res.seq in [NMRInversions.IR, NMRInversions.CPMG]
        ax1 = Axis(fig[:, 1], xlabel="time", ylabel="Signal")
        ax2 = Axis(fig[:, 2], xlabel="T", xscale=log10)

    elseif res.seq in [NMRInversions.PFG]
        ax1 = Axis(fig[:, 1], xlabel="b factor", ylabel="Signal")
        ax2 = Axis(fig[:, 2], xlabel="D (m²/s)", xscale=log10)
    end

    ax1 = Axis(fig[:, 1], xlabel="time", ylabel="Signal")
    ax2 = Axis(fig[:, 2], xscale=log10)

    draw_on_axes(ax1, ax2, res)

end

function draw_on_axes(ax1, ax2, res::NMRInversions.inv_out_1D)

    scatter!(ax1, res.x, res.y)
    lines!(ax1, res.xfit, res.yfit)
    lines!(ax2, res.X, res.f)

end