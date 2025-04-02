include("../intro.jl")

dx = 16
ncdir = datadir("output/ais/spinup/$(dx)km/corrected-hydro")
filepaths = recursive_global(ncdir, "yelmo2D.nc", 5)
target_dir = plotsdir("$(dx)km/errormap_ensemble/")

function plot_error_z(fn, hash, dx, target_dir)

    isdir(target_dir) || mkdir(target_dir)
    fn_ref = "/p/projects/megarun/ice_data/Antarctica/ANT-$(dx)KM/ANT-$(dx)KM_TOPO-RTOPO-2.0.1.nc"

    nt = ncread(fn, "time") |> length
    z_srf_now = ncread(fn, "z_srf", start = [1, 1, nt], count = [-1, -1, 1])[:, :, 1]
    f_grnd = ncread(fn, "f_grnd", start = [1, 1, nt], count = [-1, -1, 1])[:, :, 1]
    z_srf_ref = ncread(fn_ref, "z_srf")
    mask_ref = ncread(fn_ref, "mask")
    z_bed_ref = ncread(fn, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])[:, :, 1]
    z_bed_now = ncread(fn, "z_bed", start = [1, 1, nt], count = [-1, -1, 1])[:, :, 1]
    xc, yc = ncread(fn, "xc"), ncread(fn, "yc")
    diff = z_srf_now .- z_srf_ref
    rmse = sqrt(mean(diff[mask_ref .> 0] .^ 2))

    mask_ref[250:300, 150:200] .= 2

    nrows, ncols = 2, 2
    rw = 0.65
    ms = 2
    fs = 24
    fs_off = 2
    fig = Figure(size=(1100, 1200), fontsize = fs)
    axs = [Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
    [hidedecorations!(ax) for ax in [axs[1, 1], axs[1, 2], axs[2, 2]]]
    diff_cmap = (colormap = cgrad([:darkred, :white, :darkblue],
        range(0, stop = 1, length = 12),
        categorical = true), colorrange = (-5, 5))

    heatmap!(axs[1, 1], xc, yc, z_bed_ref; cmaps["z_bed"]...)
    heatmap!(axs[1, 1], xc, yc, z_srf_ref; cmaps["z_srf"]...)
    contour!(axs[1, 1], xc, yc, mask_ref .== 2, levels = [0.5],
        color = :orange, linewidth = 3)

    heatmap!(axs[1, 2], xc, yc, z_bed_now; cmaps["z_bed"]...)
    heatmap!(axs[1, 2], xc, yc, z_srf_now; cmaps["z_srf"]...)
    contour!(axs[1, 2], xc, yc, f_grnd, levels = [0.5], color = :red, linewidth = 3)

    scatter!(axs[2, 1], vec(z_srf_ref), vec(z_srf_now), markersize = ms,
        color = :gray5)
    axs[2, 1].xlabel = L"Observed ice surface elevation (km) $\,$"
    axs[2, 1].ylabel = L"Modelled ice surface elevation (km) $\,$"
    axs[2, 1].xticks = latexifyticks(0:4, 1f3)
    axs[2, 1].yticks = latexifyticks(0:4, 1f3)

    heatmap!(axs[2, 2], xc, yc, diff * 1f-2; diff_cmap...)
    contour!(axs[2, 2], xc, yc, f_grnd, levels = [0.5], color = :red, linewidth = 3)
    contour!(axs[2, 2], xc, yc, mask_ref .== 2, levels = [0.5],
        color = :orange, linewidth = 3)

    Colorbar(fig[0, 1], vertical = false, width = Relative(rw), flipaxis = true,
        label = L"Ice surface elevation (km) $\,$", ticks = latexifyticks(0:4, 1f3),
        halign = :left; cmaps["z_srf"]...)
    Colorbar(fig[0, 2], vertical = false, width = Relative(rw), flipaxis = true,
        label = L"Bedrock elevation (km) $\,$", ticks = latexifyticks(-4:2, 1f3),
        halign = :right; cmaps["z_bed"]...)
    Colorbar(fig[3, 2], vertical = false, width = Relative(rw), flipaxis = false,
        label = L"Surface elevation error at PD ($\times$ 100 m)",
        ticks = latexifyticks(-4:2:4); diff_cmap...)

    l1 = LineElement(color = :orange, linewidth = 3)
    l2 = LineElement(color = :red, linewidth = 3)
    Legend(fig[0, :], [l1, l2], [L"Observed grounding line $\,$", L"Modelled grounding line $\,$"],
        nbanks = 1, valign = :bottom)

    xt, yt = -2900, 2600
    text!(axs[1, 1], xt, yt, text = L"\textbf{(1)}", fontsize = fs + fs_off, color = :white)
    text!(axs[1, 2], xt, yt, text = L"\textbf{(2)}", fontsize = fs + fs_off, color = :white)
    text!(axs[2, 1], 100, 3900, text = L"\textbf{(3)}", fontsize = fs + fs_off, color = :gray5)
    text!(axs[2, 1], 2000, 100.0, color = :gray5, fontsize = fs + fs_off,
        text = L"RMSE = %$(round(rmse; digits = 2)) $\mathrm{m}$")
    text!(axs[2, 2], xt, yt, text = L"\textbf{(4)}", fontsize = fs + fs_off, color = :gray5)
    xlims!(axs[2, 1], -0.1, 4200)
    ylims!(axs[2, 1], -0.1, 4200)
    lines!(axs[2, 1], -1000:1000:5000, -1000:1000:5000, color = :green,
        linewidth = 3)

    rowgap!(fig.layout, 1, 10)
    rowgap!(fig.layout, 2, 10)
    rowgap!(fig.layout, 3, -45)
    colgap!(fig.layout, 10)
    save("$target_dir/z-$hash.png", fig)
end
for (i, fn) in enumerate(filepaths)
    plot_error_z(fn, i, dx, target_dir)
end