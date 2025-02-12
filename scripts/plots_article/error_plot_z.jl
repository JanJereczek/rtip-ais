include("../intro.jl")

fn = datadir("output/ais/spinup/lloyd2024-largepad/0/yelmo2D.nc")
fn_ref = datadir("output/ais/spinup/lloyd2024-largepad/0/ice_data/Antarctica/ANT-32KM/ANT-32KM_TOPO-RTOPO-2.0.1.nc")
# z_srf_ref = ncread(fn, "z_srf", start = [1, 1, 1], count = [-1, -1, 1])
z_srf_now = ncread(fn, "z_srf", start = [1, 1, 20], count = [-1, -1, 1])
f_grnd = ncread(fn, "f_grnd", start = [1, 1, 20], count = [-1, -1, 1])
z_srf_ref = ncread(fn_ref, "z_srf")
mask_ref = ncread(fn_ref, "mask")
z_bed_ref = ncread(fn, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])
z_bed_now = ncread(fn, "z_bed", start = [1, 1, 20], count = [-1, -1, 1])

nrows, ncols = 2, 2
rw = 0.7
ms = 2
fs = 20
fs_off = 2
fig = Figure(size=(1100, 1200), fontsize = fs)
axs = [Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
[hidedecorations!(ax) for ax in [axs[1, 1], axs[1, 2], axs[2, 2]]]
diff_cmap = (colormap = cgrad([:darkred, :white, :darkblue], 11, categorical = true),
    colorrange = (-5, 5))

heatmap!(axs[1, 1], view(z_bed_ref, :, :, 1); cmaps["z_bed"]...)
heatmap!(axs[1, 1], view(z_srf_ref, :, :, 1); cmaps["z_srf"]...)
contour!(axs[1, 1], mask_ref .== 2, levels = [0.5], color = :red, linewidth = 3)

heatmap!(axs[1, 2], view(z_bed_now, :, :, 1); cmaps["z_bed"]...)
heatmap!(axs[1, 2], view(z_srf_now, :, :, 1); cmaps["z_srf"]...)
contour!(axs[1, 2], view(f_grnd, :, :, 1), levels = [0.5], color = :orange, linewidth = 3)

scatter!(axs[2, 1], vec(z_srf_ref), vec(z_srf_now), markersize = ms,
    color = :gray5)
axs[2, 1].xlabel = L"PD observed ice surface elevation (km) $\,$"
axs[2, 1].ylabel = L"PD modelled ice surface elevation (km) $\,$"
axs[2, 1].xticks = latexifyticks(0:4, 1f3)
axs[2, 1].yticks = latexifyticks(0:4, 1f3)

heatmap!(axs[2, 2], view(z_srf_now .- z_srf_ref, :, :, 1) * 1f-2; diff_cmap...)
contour!(axs[2, 2], mask_ref .== 2, levels = [0.5], color = :red, linewidth = 3)
contour!(axs[2, 2], view(f_grnd, :, :, 1), levels = [0.5], color = :orange, linewidth = 3)

Colorbar(fig[0, 1], vertical = false, width = Relative(rw), flipaxis = true,
    label = L"Ice surface elevation (km) $\,$", ticks = latexifyticks(0:4, 1f3),
    halign = :left; cmaps["z_srf"]...)
Colorbar(fig[0, 2], vertical = false, width = Relative(rw), flipaxis = true,
    label = L"Bedrock elevation (km) $\,$", ticks = latexifyticks(-4:2, 1f3),
    halign = :right; cmaps["z_bed"]...)
Colorbar(fig[3, 2], vertical = false, width = Relative(rw), flipaxis = false,
    label = L"Surface elevation error w.r.t. present day ($\times$ 100 m)",
    ticks = latexifyticks(-5:5); diff_cmap...)

l1 = LineElement(color = :red, linewidth = 3)
l2 = LineElement(color = :orange, linewidth = 3)
Legend(fig[0, :], [l1, l2], [L"PD grounding line $\,$", L"Modeled grounding line $\,$"],
    nbanks = 1, valign = :bottom)

text!(axs[1, 1], 10, 170, text = L"\textbf{(a)}", fontsize = fs + fs_off, color = :white)
text!(axs[1, 2], 10, 170, text = L"\textbf{(b)}", fontsize = fs + fs_off, color = :white)
text!(axs[2, 1], 0, 3770, text = L"\textbf{(c)}", fontsize = fs + fs_off, color = :gray5)
text!(axs[2, 2], 10, 170, text = L"\textbf{(d)}", fontsize = fs + fs_off, color = :gray5)

rowgap!(fig.layout, 1, 10)
rowgap!(fig.layout, 2, 10)
rowgap!(fig.layout, 3, -40)
colgap!(fig.layout, 10)
fig

save(plotsdir("error_plot_z.pdf"), fig)
save(plotsdir("error_plot_z.png"), fig)