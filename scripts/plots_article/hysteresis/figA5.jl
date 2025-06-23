include("../../intro.jl")
file_t0 = datadir("output/ais/hyster/16km/retreat/aqef/"*
    "pmpt-normvisc-normforcing-withrestarts/0/yelmo2D.nc")
file_t0_hr = datadir("BedMachineAntarctica-v3.nc")

x, y = ncread(file_t0, "xc"), ncread(file_t0, "yc")
nx, ny = length(x), length(y)
nt1 = length(ncread(file_t0, "time"))
z_bed_t0 = ncread(file_t0, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])
z_bed_t1 = ncread(file_t0, "z_bed", start = [1, 1, nt1], count = [-1, -1, 1])
u = reshape(z_bed_t1 .- z_bed_t0, nx, ny)
itp = linear_interpolation((x, y), u, extrapolation_bc = (Flat(), Flat()))

x_hr = ncread(file_t0_hr, "x") .* 1f-3
y_hr = ncread(file_t0_hr, "y") .* 1f-3
y_hr = reverse(y_hr)
X, Y = ndgrid(x_hr, y_hr)
z_bed_t0_hr = ncread(file_t0_hr, "bed", start = [1, 1], count = [-1, -1])
z_bed_t0_hr = reverse(z_bed_t0_hr, dims = 2)
z_bed_t1_hr = z_bed_t0_hr .+ itp.(X, Y)

xy_max = 2800
i1_ais = findfirst(x_hr .> -xy_max)
i2_ais = findlast(x_hr .< xy_max)
j1_ais = findfirst(y_hr .> -xy_max)
j2_ais = findlast(y_hr .< xy_max)

i1_wais = findfirst(x_hr .> -2000)
i2_wais = findlast(x_hr .< -400)
j1_wais = findfirst(y_hr .> -1200)
j2_wais = findlast(y_hr .< 400)

i1_asb = findfirst(x_hr .> 200)
i2_asb = findlast(x_hr .< 2700)
j1_asb = findfirst(y_hr .> -2300)
j2_asb = findlast(y_hr .< 200)

mask = zeros(eltype(z_bed_t0_hr), size(z_bed_t0_hr))
view(mask, i1_wais:i2_wais, j1_wais:j2_wais) .= 1
view(mask, i1_asb:i2_asb, j1_asb:j2_asb) .= 2

bsl_0 = 0
bsl_1 = 58
lw1 = 2
lw2 = 5

framecolors = [:red, :orange, :gray20, :gray70]
nrows = 2
ncols = 3
set_theme!(theme_latexfonts())
fig = Figure(size = (1500, 1120), fontsize = 32)
axs = [Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
[hidedecorations!(ax) for ax in axs]
axs[1, 1].title = "Antarctica"
axs[1, 2].title = "West-Antarctica"
axs[1, 3].title = "Aurora and Wilkes"
axs[1, 1].ylabel = "Present day (DPR)"
axs[2, 1].ylabel = "Without ice (UPL)"
axs[1, 1].ylabelvisible = true
axs[2, 1].ylabelvisible = true
heatmap!(axs[1, 1], view(z_bed_t0_hr, i1_ais:i2_ais, j1_ais:j2_ais); cmaps["z_bed5"]...)
heatmap!(axs[2, 1], view(z_bed_t1_hr, i1_ais:i2_ais, j1_ais:j2_ais); cmaps["z_bed5"]...)
contour!(axs[1, 1], view(mask, i1_ais:i2_ais, j1_ais:j2_ais), levels = [0.01, 1.99], linewidth = lw2,
    colormap = cgrad([framecolors[1], framecolors[2]]))
contour!(axs[2, 1], view(mask, i1_ais:i2_ais, j1_ais:j2_ais), levels = [0.01, 1.99], linewidth = lw2,
    colormap = cgrad([framecolors[3], framecolors[4]]))
heatmap!(axs[1, 2], view(z_bed_t0_hr, i1_wais:i2_wais, j1_wais:j2_wais); cmaps["z_bed5"]...)
heatmap!(axs[2, 2], view(z_bed_t1_hr, i1_wais:i2_wais, j1_wais:j2_wais); cmaps["z_bed5"]...)
heatmap!(axs[1, 3], view(z_bed_t0_hr, i1_asb:i2_asb, j1_asb:j2_asb); cmaps["z_bed5"]...)
heatmap!(axs[2, 3], view(z_bed_t1_hr, i1_asb:i2_asb, j1_asb:j2_asb); cmaps["z_bed5"]...)
contour!(axs[1, 2], view(z_bed_t0_hr, i1_wais:i2_wais, j1_wais:j2_wais), levels = [58],
    colormap = cgrad([:gray20, :gray20]), linewidth = lw1)
contour!(axs[2, 2], view(z_bed_t1_hr, i1_wais:i2_wais, j1_wais:j2_wais), levels = [58],
    colormap = cgrad([:gray20, :gray20]), linewidth = lw1)
contour!(axs[1, 3], view(z_bed_t0_hr, i1_asb:i2_asb, j1_asb:j2_asb), levels = [58],
    colormap = cgrad([:gray20, :gray20]), linewidth = lw1)
contour!(axs[2, 3], view(z_bed_t1_hr, i1_asb:i2_asb, j1_asb:j2_asb), levels = [58],
    colormap = cgrad([:gray20, :gray20]), linewidth = lw1)

for (ax, c) in zip(vec([axs[i, j] for j in 2:ncols, i in 1:nrows]), framecolors)
    ax.bottomspinecolor = c
    ax.topspinecolor = c
    ax.leftspinecolor = c
    ax.rightspinecolor = c
    ax.spinewidth = lw2
end

Colorbar(fig[nrows+1, 1:2], label = "Bed elevation (km)", width = Relative(0.5), flipaxis = false,
    vertical = false, ticks = latexifyticks(-4:2, 1e3); cmaps["z_bed5"]...)
elem_1 = LineElement(color = :gray20, linewidth = 2)
Legend(fig[nrows+1, 2:3], [elem_1], [L"$z_b = 58 \, \mathrm{m}$ isoline"], valign = :top,
    framevisible = false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)
save(plotsdir("16km/bedrock/bedmachine3.png"), fig)
save(plotsdir("16km/bedrock/bedmachine3.pdf"), fig)