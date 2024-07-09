function select_peaks(file::string = pick_file(pwd()))

    invres = readresults(file)
    dir = invres.dir
    indir = invres.indir
    f = invres.f
    F = collect(reshape(f, length(dir), length(indir))')

    x = collect(1:length(indir))
    y = collect(1:length(dir))
    z = copy(F)

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

    # Peak selection counter
    peak_counter = 0

    # Create figure and axes
    GLMakie.activate!(float=true) # Keep plot window on top
    f = Figure(size=(700, 700))
    m = GridLayout(f[1, 1:2])
    n = GridLayout(f[2:3, 2])
    b = GridLayout(f[2, 1])
    # ax3d = Axis3(n[1, 1])

    # Title textbox
    tb = Textbox(f[3, 1], placeholder="Insert a title for the plot, then press enter.",
        width=300, reset_on_defocus=true)


    axmain = Axis(m[2, 1])
    axtop = Axis(m[1, 1])
    axright = Axis(m[2, 2])


    # Create another axis just to show the actual limits of the plot, 
    # since array indices are used on the previous axis
    abs(log10(dir[1])) + abs(log10(dir[end]))

    axmainr = Axis(m[2, 1],
        xscale=log10, yscale=log10,
        xlabel=L"T_1 \, \textrm{(s)}", ylabel=L"T_2 \,\textrm{(s)}",
        xlabelsize=23, ylabelsize=23,
        limits=(indir[1], indir[end], dir[1], dir[end]))

    colsize!(m, 1, Fixed(300))
    rowsize!(m, 2, Fixed(300))
    rowsize!(m, 1, Aspect(1, 0.25))
    colsize!(m, 2, Aspect(1, 1))

    Makie.deactivate_interaction!(axmain, :rectanglezoom) # Disable zoom
    Makie.deactivate_interaction!(axmainr, :rectanglezoom) # Disable zoom
    linkyaxes!(axmain, axright)
    linkxaxes!(axmain, axtop)
    # linkxaxes!(axmain, ax3d)
    # linkyaxes!(axmain, ax3d)


    hidedecorations!(axtop)
    hidedecorations!(axright)
    hidedecorations!(axmain)
    # hidedecorations!(ax3d)

    colgap!(m, 0)
    rowgap!(m, 0)

    function drawplots()

        empty!(axmain)
        empty!(axtop)
        empty!(axright)
        # empty!(ax3d)

        # Static plots
        contourf!(axmain, z, colormap=:tempo, levels=50)
        lines!(axtop, vec(sum(z, dims=2)), colormap=:tempo, colorrange=(1, 10), color=5, alpha=0.5)
        lines!(axright, vec(sum(z, dims=1)), 1:length(dir), colormap=:tempo, colorrange=(1, 10), color=5, alpha=0.5)
        # Plot diagonal line
        lines!(axmain, [(0, 0), (length(indir), length(dir))], color=:black, linewidth=1)

        # surface!(ax3d, x, y, z, colormap=:tempo)

        # Dynamic plots
        contourf!(axmain, spo, mode=:relative, levels=range(0.01, 1, 20), extendlow=:transparent, extendhigh=:transparent, colormap=:blues)
        lines!(axtop, indir_dist, colormap=:tab10, colorrange=(1, 10), color=8)
        lines!(axright, dir_dist, 1:length(dir), colormap=:tab10, colorrange=(1, 10), color=8)

        # surface!(ax3d, x, y, spo, colormap=:blues,
            # colorrange=(minimum(filter(x -> x != 0, vec(z))), maximum(z)),
            # lowclip=:transparent)

        # Plot the points in the selection vector
        scatter!(axmain, selection, colormap=:tab10, colorrange=(1, 10), color=8)
        lines!(axmain, polygon, colormap=:tab10, colorrange=(1, 10), color=8)
    end

    drawplots()

    # Pick points by clicking the plot
    spoint = select_point(axmain.scene, marker=:circle)

    # Update selection Observable by pushing the selected point to it
    on(spoint) do _
        push!(selection[], spoint[])
        selection[] = selection[] # Update explicitly so that the plots are updated automatically

        if size(selection[], 1) > 2
            polygon[] = vcat(selection[], [selection[][1]])
            mask[] = [inpolygon(p, polygon[]; in=1, on=1, out=0) for p in points]
        end

    end


    # Label current selection button
    labelb = Button(b[1, 1]; label="Label current \nselection")

    on(labelb.clicks) do _

        if size(selection[], 1) < 3
            @warn("You need to make a selection.")
        else

            peak_counter += 1

            #draw current polygon
            lines!(axmain, polygon[], linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=peak_counter, alpha=0.9)

            # Plot peak label
            xc = vec(sum(spo[], dims=2)) ⋅ collect(1:size(spo[], 1)) / sum(spo[])
            yc = vec(sum(spo[], dims=1)) ⋅ collect(1:size(spo[], 2)) / sum(spo[])
            scatter!(axmain, xc, yc, markersize=15,
                marker=collect('a':'z')[peak_counter],
                colormap=:tab10, colorrange=(1, 10),
                color=peak_counter,
                glowcolor=:white, glowwidth=4)

            T1 = vec(sum(spo[], dims=2)) ⋅ indir / sum(spo[])
            T2 = vec(sum(spo[], dims=1)) ⋅ dir / sum(spo[])

            text!(axmain, 2.5 * (length(indir) / 100), (99 - 13 * (peak_counter - 1)) * (length(dir) / 100),
                text=
                collect('a':'z')[peak_counter] *
                " : T₁/T₂ = $(round(T1/T2 , digits=1)) \n     Volume = $(round(sum(spo[])/sum(z) *100 ,digits = 1))%",
                align=(:left, :top), colormap=:tab10, colorrange=(1, 10), color=peak_counter)

            lines!(axtop, indir_dist[], linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=peak_counter, alpha=0.5)
            lines!(axright, dir_dist[], 1:length(dir), linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=peak_counter, alpha=0.5)

            # Add info to the file
            open(file, "a") do file
                println(file, "Selection " * collect('a':'z')[peak_counter] *
                              " : T1=$T1 , T2=$T2 , T1/T2=$(T1/T2) , VolFraction=$(sum(spo[])/sum(z)) , Polygon Selection : " *
                              join(reduce(vcat, [[polygon[][i][1], polygon[][i][2]] for i in axes(polygon[], 1)]), ", ")
                )
            end


            # Clear for next selection
            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))

            reset_limits!(axmain)
        end

    end

    # Collect points button
    singlepeakb = Button(b[2, 1]; label="Select \nsingle peak")

    on(singlepeakb.clicks) do _

        if size(selection[], 1) < 3
            spo[] = F
        end

        T1 = vec(sum(spo[], dims=2)) ⋅ indir / sum(spo[])
        T2 = vec(sum(spo[], dims=1)) ⋅ dir / sum(spo[])

        drawplots()
        peak_counter = 0

        text!(axmain, 20, 90,
            text="T₁/T₂ = $(round(T1/T2 , digits=1))",
            align=(:left, :top), offset=(4, 0))


        # Clear for next selection
        selection[] = Point{2,Float32}[]
        polygon[] = Point{2,Float32}[]
        mask[] = zeros(size(points))

    end

    # Clear current selection button
    clearb = Button(b[1, 2]; label="Clear current \nselection")
    on(clearb.clicks) do _
        selection[] = Point{2,Float32}[]
        polygon[] = Point{2,Float32}[]
        mask[] = zeros(size(points))

    end

    # Reset Selection button
    resetb = Button(b[2, 3]; label="Reset all")

    on(resetb.clicks) do _
        z = copy(F)
        spo = @lift(($mask .* $z))
        selection[] = Point{2,Float32}[]
        polygon[] = Point{2,Float32}[]
        mask[] = zeros(size(points))
        drawplots()
        peak_counter = 0

        # Rewrite the file, keeping only the first 7 lines
        lines = readlines(file)
        lines_to_keep = lines[1:min(7, end)]
        open(file, "w") do io
            for line in lines_to_keep
                write(io, line * "\n")
            end
        end

    end

    # Clear non-selected
    artclearb = Button(b[1, 3]; label="Delete \nnon-selected")

    on(artclearb.clicks) do _

        if size(selection[], 1) < 3
            @warn("You need to make a selection.")
        else
            # Rewrite the file, keeping only the first 7 lines plus polygon selection
            lines = readlines(file)
            lines_to_keep = lines[1:min(7, end)]
            open(file, "w") do io
                for line in lines_to_keep
                    write(io, line * "\n")
                end

                write(io, "Keep data within: ")
                write(io,
                    join(reduce(vcat, [[polygon[][i][1], polygon[][i][2]] for i in axes(polygon[], 1)]), ", ") * "\n"
                )
            end

            z = spo[]
            spo = @lift(($mask .* $z))
            selection[] = Point{2,Float32}[]
            polygon[] = Point{2,Float32}[]
            mask[] = zeros(size(points))
            drawplots()
            peak_counter = 0
        end
    end

    # Export figure button
    exportb = Button(b[2, 2]; label="Export")

    on(exportb.clicks) do _
        if isnothing(tb.stored_string[])
            pubfig(title="")
        else
            pubfig(title=tb.stored_string[])
        end
    end



    reset_limits!(axmain)
    f

end


function pubfig(file="inversion_results.txt"; title="", ppu=2)

    invres = readresults(file)
    dir = invres.dir
    indir = invres.indir
    f = invres.f

    x = collect(1:length(indir))
    y = collect(1:length(dir))
    z = collect(reshape(f, length(dir), length(indir))')

    #Look if there is any deleted selection
    open(file) do io

        readuntil(io, "Keep data within: ")
        if !eof(io) # if you don't hit the end of the file...
            # ...the rest of the line is the polygon selection, read it and parse it into Point2f format
            a = parse.(Float32, split(readline(io), ", "))
            apolygon = [Point2f(a[i], a[i+1]) for i in 1:2:length(a)]

            #Create a matrix for all the discrete points in the space
            points = [[i, j] for i in x, j in y]

            mask = [inpolygon(p, apolygon; in=1, on=1, out=0) for p in points]

            # Keep selected points only in the z matrix
            z .= mask .* z

        end
    end

    m = Figure(size=(450, 460))
    axmain = Axis(m[2, 1], limits=(x[1], x[end], y[1], y[end]))
    axtop = Axis(m[1, 1], limits=((x[1], x[end]), (0, nothing)), title=title, titlesize=17)
    axright = Axis(m[2, 2], limits=((0, nothing), (y[1], y[end])))

    axmainr = Axis(m[2, 1],
        xscale=log10, yscale=log10,
        xlabel=L"T_1 \, \textrm{(s)}", ylabel=L"T_2 \,\textrm{(s)}",
        xlabelsize=23, ylabelsize=23,
        limits=(indir[1], indir[end], dir[1], dir[end])
    )

    rowsize!(m.layout, 2, 300)
    colsize!(m.layout, 1, 300)
    rowsize!(m.layout, 1, Aspect(1, 0.25))
    colsize!(m.layout, 2, Aspect(1, 1))

    hidedecorations!(axtop)
    hidedecorations!(axright)
    hidedecorations!(axmain)
    colgap!(m.layout, 0)
    rowgap!(m.layout, 0)

    # Plots
    contourf!(axmain, x, y, z, colormap=:tempo, levels=50)
    lines!(axtop, vec(sum(z, dims=2)), colormap=:tempo, colorrange=(1, 10), color=5, alpha=0.5)
    lines!(axright, vec(sum(z, dims=1)), 1:length(dir), colormap=:tempo, colorrange=(1, 10), color=5, alpha=0.5)
    # Plot diagonal line
    lines!(axmain, [(0, 0), (length(indir), length(dir))], color=:black, linewidth=1)


    # Look for selections in the text file and print them on the plot
    open(file) do io
        i = 1
        while true
            readuntil(io, "Selection " * collect('a':'z')[i] * " : ")

            if eof(io)
                break
            end

            # Read useful data from file
            line = readuntil(io, "Polygon Selection : ")
            matches = match(r"T1=([\d\.e-]+)\s*,\s*T2=([\d\.e-]+)\s*,.*VolFraction=([\d\.e-]+)", line)
            T1 = parse(Float64, matches.captures[1])
            T2 = parse(Float64, matches.captures[2])
            VolFraction = parse(Float64, matches.captures[3])

            # The rest of the line is the polygon selection, read it and parse it into Point2f format
            x = parse.(Float32, split(readline(io), ", "))
            apolygon = [Point2f(x[i], x[i+1]) for i in 1:2:length(x)]

            #draw current polygon
            lines!(axmain, apolygon, linestyle=:dash, colormap=:tab10, colorrange=(1, 10), color=i, alpha=0.9)


            #add peak "center of mass" label
            scatter!(axmainr, T1, T2, markersize=15,
                marker=collect('a':'z')[i],
                colormap=:tab10, colorrange=(1, 10),
                color=i,
                glowcolor=:white, glowwidth=4)

            #add text to plot
            text!(axmain, 2.5 * (length(indir) / 100), (99 - 13 * (i - 1)) * (length(dir) / 100),
                text=
                collect('a':'z')[i] *
                " : T₁/T₂ = $(round(T1/T2 , digits=1)) \n     Volume = $(round(VolFraction*100 ,digits = 1))%",
                align=(:left, :top), colormap=:tab10, colorrange=(1, 10), color=i
            )

            # Look for the next selection (if it exists)
            i += 1
        end
    end

    Makie.save("Results.png", m, px_per_unit=ppu)

end