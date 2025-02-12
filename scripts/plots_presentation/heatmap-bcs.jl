include("../intro.jl")

file2D = datadir("output/ais/hyster/aqef/20K/1/yelmo2D.nc")
nx, ny = 191, 191
ghf = reshape(ncread(file2D, "Q_geo", start=[1, 1, 1], count=[-1, -1, 1]), nx, ny)
T_srf = reshape(ncread(file2D, "T_srf", start=[1, 1, 1], count=[-1, -1, 1]), nx, ny)
T_shlf = reshape(ncread(file2D, "T_shlf", start=[1, 1, 1], count=[-1, -1, 1]), nx, ny)

fs = 20
set_theme!(theme_latexfonts())
fig = Figure(size = (500, 800), fotsize = fs)
axs = [Axis(fig[i, 1], aspect = AxisAspect(1)) for i in 1:2]
[hidedecorations!(ax) for ax in axs]
heatmap!(axs[1], T_srf; cmaps["T_srf"]...)
heatmap!(axs[2], T_shlf; cmaps["T_shlf"]...)
Colorbar(fig[1, 2], label = "Surface temperature (K)", height = Relative(0.7); cmaps["T_srf"]...)
Colorbar(fig[2, 2], label = "Shelf temperature (K)", height = Relative(0.7); cmaps["T_shlf"]...)
rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
fig
save(plotsdir("bc.png"), fig)