include("../intro.jl")

fn_m2 = datadir("output/ais/ramps/16km/steps-sigmarange/0/fastisostasy.nc")
fn_m1 = datadir("output/ais/ramps/16km/steps-sigmarange/1/fastisostasy.nc")
fn_nom = datadir("output/ais/ramps/16km/steps-sigmarange/2/fastisostasy.nc")
fn_p1 = datadir("output/ais/ramps/16km/steps-sigmarange/3/fastisostasy.nc")
fn_p2 = datadir("output/ais/ramps/16km/steps-sigmarange/4/fastisostasy.nc")
filenames = [fn_m2, fn_m1, fn_nom, fn_p1, fn_p2]

T = Float32
xc = ncread(fn_nom, "xc")
yc = ncread(fn_nom, "yc")
eta = [ncread(fn, "log10_eta_eff") for fn in filenames]
Tlitho = ncread(fn_nom, "He_lith")

file_yelmo = datadir("output/ais/ramps/16km/steps-sigmarange/2/yelmo2D.nc")
h_ice = ncslice(file_yelmo, "H_ice", 1)
ice_mask = h_ice .> 0.5
xc_yelmo = ncread(file_yelmo, "xc")
yc_yelmo = ncread(file_yelmo, "yc")

ncols = 4
set_theme!(theme_latexfonts())
fig = Figure(size = (2000, 650), fontsize = 26)
axs = [Axis(fig[1, j], aspect = DataAspect()) for j in 1:ncols]
[hidedecorations!(ax) for ax in axs]
copts_litho = (colorrange = (0, 300), colormap = cgrad(:inferno,
    range(0, stop = 1, length = 13), rev = true, categorical = true))
copts_visco = (colorrange = (19, 22.5), colormap = cgrad(:jet,
    range(0, stop = 1, length = 15), categorical = true, rev = true))
Colorbar(fig[2, 1:2], vertical = false, width = Relative(0.4),
    ticks = ([0, 100, 200, 300], latexify.([0, 100, 200, 300])),
    label = "Lithospheric thickness (km)", flipaxis=false; copts_litho...)
Colorbar(fig[2, 3:ncols], vertical = false, width = Relative(0.4),
    ticks = ([19, 20, 21, 22], latexify.([19, 20, 21, 22])),
    label = "Log10 effective viscosity (Pa s)", flipaxis=false; copts_visco...)

lw = 4
lin1 = lines!(axs[1], -1:1, -1:1, color = :black, linewidth = lw,
    label = "Modelled PD ice extent")
heatmap!(axs[1], xc, yc, Tlitho; copts_litho...)
heatmap!(axs[2], xc, yc, eta[1]; copts_visco...)
heatmap!(axs[3], xc, yc, eta[3]; copts_visco...)
heatmap!(axs[4], xc, yc, eta[5]; copts_visco...)

for i in 1:ncols
    contour!(axs[i], xc_yelmo, yc_yelmo, ice_mask, levels = [0.5], linewidth = lw,
        color = :black)
    xlims!(axs[i], extrema(xc_yelmo) .+ (100, -100))
    ylims!(axs[i], extrema(yc_yelmo) .+ (100, -100))
end
Legend(fig[2, 2:3], axs[1], valign = :top, nbanks = 2, fontsize = 26)

axs[1].title = L"\textbf{a} \quad Lithospheric thickness $\,$"
axs[2].title = L"\textbf{b} \quad Lowest effective viscosity ($-2 \, σ$)"
axs[3].title = L"\textbf{c} \quad Nominal effective viscosity $\,$"
axs[4].title = L"\textbf{d} \quad Highest effective viscosity ($+2 \, σ$)"
rowgap!(fig.layout, 1)
rowsize!(fig.layout, 2, 20)
save(plotsdir("16km/solidearth.png"), fig)
