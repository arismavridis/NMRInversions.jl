using CairoMakie

x = range(1, 10, 50)
y = copy(x)

C1 = [5.15, 4.2]
C2 = [7, 6.4]
C3 = C1 + [cosd(-60) -sind(-60); sind(-60) cosd(-60)] * (C2 - C1)
C = [C1, C2, C3]

# Define the colormaps
purple = [RGBAf(0.514, 0.227, 0.584, 0.0), RGBAf(0.514, 0.227, 0.584, 1.0)];
green = [RGBAf(0.298, 0.569, 0.255, 0.0), RGBAf(0.298, 0.569, 0.255, 1.0)];
red = [RGBAf(0.882, 0.141, 0.188, 0.0), RGBAf(0.882, 0.141, 0.188, 1.0)];
jblue =  RGBAf(0.235, 0.443, 0.698, 1.0);
jcolors = [red, green, purple]


begin
    f = Figure(size = (500,500),backgroundcolor=:transparent)
    ax = Axis(f[1, 1], aspect=1,spinewidth=4)
    lines!(ax, -1:length(x)+1,-1:length(y)+1,color=jblue,linewidth=4)

    ax.backgroundcolor =:transparent
    ax.yticks= [-10]
    ax.xticks= [-10]

    for i in 1:3
        A = [exp.((-((x - C[i][1])^2 + (y - C[i][2])^2))/1.75) / 100 for x in x, y in y]
        heatmap!(ax, A, colormap=jcolors[i])
    end

    f
end
