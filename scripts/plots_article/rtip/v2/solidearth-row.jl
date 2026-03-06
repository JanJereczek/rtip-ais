include("../../../intro.jl")

fn_m2 = datadir("output/ais/v2/hyster/retreat/aqef/refm2slow/0/fastisostasy.nc")
fn_m1 = datadir("output/ais/v2/hyster/retreat/aqef/refm1slow/0/fastisostasy.nc")
fn_nom = datadir("output/ais/v2/hyster/retreat/aqef/refnomslow/0/fastisostasy.nc")
fn_p1 = datadir("output/ais/v2/hyster/retreat/aqef/refp1slow/0//fastisostasy.nc")
fn_p2 = datadir("output/ais/v2/hyster/retreat/aqef/refp2slow/0//fastisostasy.nc")
filenames = [fn_m2, fn_m1, fn_nom, fn_p1, fn_p2]

T = Float32
xc = ncread(fn_nom, "xc")
yc = ncread(fn_nom, "yc")
eta = [ncread(fn, "log10_eta_eff") for fn in filenames]
Tlitho = ncread(fn_nom, "He_lith")

file_yelmo = datadir("output/ais/v2/hyster/retreat/aqef/refnomslow/0/yelmo2D.nc")
h_ice = ncslice(file_yelmo, "H_ice", 1)
ice_mask = h_ice .> 0.5
xc_yelmo = ncread(file_yelmo, "xc")
yc_yelmo = ncread(file_yelmo, "yc")

nrows, ncols = 1, 3
set_theme!(theme_latexfonts())
fig = Figure(size = (1500, 700), fontsize = 30)
axs = [Axis(fig[1, j], aspect = DataAspect()) for j in 1:ncols]
[hidedecorations!(ax) for ax in axs]
copts_visco = (colorrange = (19, 22.5), colormap = cgrad(:jet,
    range(0, stop = 1, length = 15), categorical = true, rev = true))
Colorbar(fig[2, 1:2], vertical = false, width = Relative(0.4), height = 20,
    ticks = ([19, 20, 21, 22], latexify.([19, 20, 21, 22])),
    label = "Log10 effective viscosity (Pa s)", flipaxis=false; copts_visco...)

lw = 4
heatmap!(axs[1], xc, yc, eta[1]; copts_visco...)
heatmap!(axs[2], xc, yc, eta[3]; copts_visco...)
heatmap!(axs[3], xc, yc, eta[5]; copts_visco...)

for j in 1:ncols
    contour!(axs[j], xc_yelmo, yc_yelmo, ice_mask, levels = [0.5], linewidth = lw,
        color = :black)
    xlims!(axs[j], extrema(xc_yelmo) .+ (100, -100))
    ylims!(axs[j], extrema(yc_yelmo) .+ (100, -100))
end
elem_1 = LineElement(color = :black, linewidth = lw)
Legend(fig[2, 2:3], [elem_1], ["Modelled PD ice extent"], valign = :top,
    framevisible = true)

axs[1].title = L"$-2 \, \sigma$"
axs[2].title = L"$0 \, \sigma$"
axs[3].title = L"$+2 \, \sigma$"

rowgap!(fig.layout, 15)
colgap!(fig.layout, 5)
fig
save(plotsdir("v2/rtip/solidearth-row.png"), fig)
