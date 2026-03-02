file_yelmo = datadir("output/ais/v2/ramps/steps-sigmarange/2/yelmo2Dsm.nc")
h_ice = ncslice(file_yelmo, "H_ice", 1)
ice_mask = h_ice .> 0.5
xc, yc = ncread(file_yelmo, "xc"), ncread(file_yelmo, "yc")
nx, ny = length(xc), length(yc)

file = "/p/projects/megarun/jan/rtip-ais/data/output/ais/v2/hyster/retreat/aqef/elra/0/isostasy_data/earth_structure/yelmo/ANT-16KM_GIA_HR24.nc"

sigma = ncread(file, "stddev_log10_visc")
Te = ncread(file, "litho_thickness")
zc = ncread(file, "zc")
re = maximum(zc)
sigma2D = fill(NaN, nx, ny)

for i in 1:nx, j in 1:ny
    @show sigma[i, j, :]
    sigma[i, j] = mean(sigma[i, j, re-Te[i, j] .> zc .>= re-400f3])
    sigma2D[i, j] = sigma[i, j]
end

set_theme!(theme_latexfonts())
stdfig = Figure(size = (800, 700), fontsize = 24)
ax = Axis(stdfig[1, 1], aspect = DataAspect())
hidedecorations!(ax)
hm = heatmap!(ax, xc, yc, sigma2D, colormap=:viridis)
contour!(ax, xc, yc, ice_mask, levels = [0.5], color = :red, linewidth = 2)
xlims!(ax, extrema(xc))
ylims!(ax, extrema(yc))
vlines!(ax, 1f10, color = :red, label = "Modelled PD ice extent", linewidth = 2)
Colorbar(stdfig[1, 2], hm, height = Relative(0.5), 
    label = L"\mathrm{log}_{10}(\sigma) \: \mathrm{(Pa \, s)}")
axislegend(ax, position = :lb)
save(plotsdir("v2/rtip/stddev_visc.png"), stdfig)