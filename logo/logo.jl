using CairoMakie

begin
    x = range(1, 10, 50)
    y = copy(x)
    C1 = [5.15, 4.2]
    C2 = [7, 6.4]
    C3 = C1 + [cosd(-60) -sind(-60); sind(-60) cosd(-60)] * (C2 - C1)
    C = [C1, C2, C3]
end

# Define the colormaps
begin
    jpurple = [RGBAf(0.584, 0.345, 0.698, 0.0), RGBAf(0.584, 0.345, 0.698, 1.0)]
    jgreen = [RGBAf(0.220, 0.569, 0.149 , 0.0), RGBAf(0.220, 0.569, 0.149, 1.0)]
    jred = [RGBAf(0.796, 0.235, 0.200, 0.0), RGBAf(0.796, 0.235, 0.200, 1.0)]
    jblue = RGBAf(0.251 , 0.388, 0.847, 1.0)
    jcolors = [jred, jgreen, jpurple]
end

begin
    f = Figure(size=(500, 500), backgroundcolor=:transparent)
    ax = Axis(f[1, 1], aspect=1, spinewidth=7.5,bottomspinecolor=jblue, leftspinecolor=jblue, topspinecolor=jblue, rightspinecolor=jblue)
    lines!(ax, -1:length(x)+1, -1:length(y)+1, color=jblue, linewidth=7.5)

    ax.backgroundcolor = :transparent
    ax.yticks = [-10]
    ax.xticks = [-10]

    ax.xlabel = "𝑇₁"
    ax.ylabel = "𝑇₂"
    ax.xlabelcolor = jblue
    ax.ylabelcolor = jblue
    ax.xlabelsize = 60
    ax.ylabelsize = 60
    ax.ylabelpadding = -17
    ax.xlabelpadding = -17


    for i in 1:3
        A = [exp.((-((x - C[i][1])^2 + (y - C[i][2])^2)) / 1.75) / 100 for x in x, y in y]
        heatmap!(ax, A, colormap=jcolors[i])
    end

    f
end
