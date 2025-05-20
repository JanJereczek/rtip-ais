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

nrows, ncols = 2, 2
set_theme!(theme_latexfonts())
fig = Figure(size = (1050, 1250), fontsize = 40)
axs = [Axis(fig[i+1, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
[hidedecorations!(ax) for ax in axs]
copts_litho = (colorrange = (0, 300), colormap = cgrad(:inferno,
    range(0, stop = 1, length = 13), rev = true, categorical = true))
copts_visco = (colorrange = (19, 22.5), colormap = cgrad(:jet,
    range(0, stop = 1, length = 15), categorical = true, rev = true))
Colorbar(fig[1, 1], vertical = false, width = Relative(0.8), height = 20,
    ticks = ([0, 100, 200, 300], latexify.([0, 100, 200, 300])),
    label = "Lithospheric thickness (km)", flipaxis=true; copts_litho...)
Colorbar(fig[4, 1:ncols], vertical = false, width = Relative(0.4), height = 20,
    ticks = ([19, 20, 21, 22], latexify.([19, 20, 21, 22])),
    label = "Log10 effective viscosity (Pa s)", flipaxis=false; copts_visco...)

lw = 4
heatmap!(axs[1, 1], xc, yc, Tlitho; copts_litho...)
heatmap!(axs[1, 2], xc, yc, eta[1]; copts_visco...)
heatmap!(axs[2, 1], xc, yc, eta[3]; copts_visco...)
heatmap!(axs[2, 2], xc, yc, eta[5]; copts_visco...)

for i in 1:nrows, j in 1:ncols
    contour!(axs[i, j], xc_yelmo, yc_yelmo, ice_mask, levels = [0.5], linewidth = lw,
        color = :black)
    xlims!(axs[i, j], extrema(xc_yelmo) .+ (100, -100))
    ylims!(axs[i, j], extrema(yc_yelmo) .+ (100, -100))
end
elem_1 = LineElement(color = :black, linewidth = lw)
Legend(fig[0, 2], [elem_1], ["Modelled PD ice extent"], valign = :top,
    framevisible = true)

# axs[1, 1].title = L"\textbf{d} \quad Lithospheric thickness $\,$"
# axs[1, 2].title = L"\textbf{e} \quad Lowest effective viscosity ($-2 \, σ$)"
# axs[2, 1].title = L"\textbf{f} \quad Nominal effective viscosity $\,$"
# axs[2, 2].title = L"\textbf{g} \quad Highest effective viscosity ($+2 \, σ$)"

text!(axs[1, 1], -2600, 2400, text="d", font = :bold)
text!(axs[1, 2], -2600, 2400, text="e", font = :bold)
text!(axs[2, 1], -2600, 2400, text="f", font = :bold)
text!(axs[2, 2], -2600, 2400, text="g", font = :bold)

rowgap!(fig.layout, 15)
colgap!(fig.layout, 0)
rowgap!(fig.layout, 1, -40)
rowgap!(fig.layout, 2, -5)
rowsize!(fig.layout, 1, 60)
rowsize!(fig.layout, 0, 1)
colsize!(fig.layout, 1, 500)
colsize!(fig.layout, 2, 500)
save(plotsdir("16km/solidearth.png"), fig)
