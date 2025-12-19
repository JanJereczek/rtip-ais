include("../../../intro.jl")

dx = 16
ncdir = datadir("output/ais/v2/spinup/16km/minvisc/3")
filepaths = recursive_global(ncdir, "yelmo2D.nc", 5)
target_dir = plotsdir("v2/hysteresis/")

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

    set_theme!(theme_latexfonts())
    nrows, ncols = 2, 4
    rw = 0.8
    ms = 2
    fs = 23
    fs_off = 2
    fig = Figure(size=(1600, 900), fontsize = fs)
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
    axs[2, 1].xlabel = L"Observed $z_\mathrm{srf}$ (km)"
    axs[2, 1].ylabel = L"Modelled $z_\mathrm{srf}$ (km)"
    axs[2, 1].xticks = latexifyticks(0:4, 1f3)
    axs[2, 1].yticks = latexifyticks(0:4, 1f3)

    heatmap!(axs[2, 2], xc, yc, diff * 1f-2; diff_cmap...)
    contour!(axs[2, 2], xc, yc, f_grnd, levels = [0.5], color = :red, linewidth = 3)
    contour!(axs[2, 2], xc, yc, mask_ref .== 2, levels = [0.5],
        color = :orange, linewidth = 3)

    Colorbar(fig[0, 1], vertical = false, width = Relative(rw), flipaxis = true,
        label = L"$z_\mathrm{srf}$ (km)", ticks = latexifyticks(0:4, 1f3); cmaps["z_srf"]...)
    # Colorbar(fig[0, 2], vertical = false, width = Relative(rw), flipaxis = true,
    #     label = L"$z_\mathrm{b}$ (km) $\,$", ticks = latexifyticks(-4:2, 1f3),
    #     halign = :right; cmaps["z_bed"]...)
    Colorbar(fig[3, 2], vertical = false, width = Relative(rw), flipaxis = false,
        label = L"$z_\mathrm{srf}$ error at PD ($\times$ 100 m)",
        ticks = latexifyticks(-4:2:4); diff_cmap...)

    l1 = LineElement(color = :orange, linewidth = 3)
    l2 = LineElement(color = :red, linewidth = 3)
    Legend(fig[0, 1:2], [l1, l2], [L"Observed $\,$", L"Modelled $\,$"], halign = 0.8,
        nbanks = 1, valign = :bottom)

    xt, yt = -2900, 2500
    text!(axs[1, 1], xt, yt, text = L"\textbf{(a)}", fontsize = fs + fs_off, color = :white)
    text!(axs[1, 2], xt, yt, text = L"\textbf{(b)}", fontsize = fs + fs_off, color = :white)
    text!(axs[2, 1], 100, 3800, text = L"\textbf{(e)}", fontsize = fs + fs_off, color = :gray5)
    text!(axs[2, 1], 1500, 100.0, color = :gray5, fontsize = fs + fs_off,
        text = "RMSE = $(round(rmse; digits = 2)) m")
    text!(axs[2, 2], xt, yt, text = L"\textbf{(f)}", fontsize = fs + fs_off, color = :gray5)
    xlims!(axs[2, 1], -0.1, 4200)
    ylims!(axs[2, 1], -0.1, 4200)
    lines!(axs[2, 1], -1000:1000:5000, -1000:1000:5000, color = :green,
        linewidth = 3)

    # rowgap!(fig.layout, 1, 10)
    # rowgap!(fig.layout, 2, 10)
    # rowgap!(fig.layout, 3, -60)
    # colgap!(fig.layout, 10)
    # save("$target_dir/figA2.png", fig)
    # save("$target_dir/figA2.pdf", fig)

    fn_ref = "/p/projects/megarun/ice_data/Antarctica/ANT-$(dx)KM/ANT-$(dx)KM_VEL-R11-2.nc"
    nt = ncread(fn, "time") |> length
    uxy_s_now = ncread(fn, "uxy_s", start = [1, 1, nt], count = [-1, -1, 1])
    uxy_s_ref = ncread(fn_ref, "uxy_srf")
    z_bed_ref = ncread(fn, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])
    z_bed_now = ncread(fn, "z_bed", start = [1, 1, nt], count = [-1, -1, 1])
    z_srf_ref = ncread(fn, "z_srf", start = [1, 1, 1], count = [-1, -1, 1])
    xc, yc = ncread(fn, "xc"), ncread(fn, "yc")
    X, Y = ncread(fn, "x2D"), ncread(fn, "y2D")
    diff = view(uxy_s_now, :, :, 1) .- uxy_s_ref
    idx = CartesianIndices(diff)
    mask_now = idx[uxy_s_now .!= 0]
    mask_ref = idx[z_srf_ref .!= 0]
    rmse = sqrt(mean(diff[mask_now] .^ 2))

    tol = 1f-8
    [hidedecorations!(ax) for ax in [axs[1, 3], axs[1, 4], axs[2, 3]]]
    diff_cmap = (colormap = cgrad(
        [:darkred, :lightsalmon, :white, :white, :white, :skyblue, :darkblue],
        11, categorical = true, rev = true),
        lowclip = :darkblue, highclip = :darkred, colorrange = (-3, 3))

    uxy_s_ref[uxy_s_ref .== 0] .= NaN
    heatmap!(axs[1, 3], xc, yc, view(z_bed_ref, :, :, 1); cmaps["z_bed"]...)
    heatmap!(axs[1, 3], X[mask_ref], Y[mask_ref],
        log10.(view(uxy_s_ref, mask_ref, 1) .+ tol); cmaps["uxy_s"]...)
    heatmap!(axs[1, 4], xc, yc, view(z_bed_now, :, :, 1); cmaps["z_bed"]...)
    heatmap!(axs[1, 4], X[mask_now], Y[mask_now],
        log10.(view(uxy_s_now, mask_now, 1) .+ tol); cmaps["uxy_s"]...)
    heatmap!(axs[2, 3], xc, yc, sign.(diff) .* log10.(abs.(diff .+ tol)); diff_cmap...)

    scatter!(axs[2, 4], log10.(vec(uxy_s_ref) .+ tol), log10.(vec(uxy_s_now) .+ tol),
        markersize = ms, color = :gray5)
    axs[2, 4].xlabel = L"Observed $u_\mathrm{s}$ ($\mathrm{m \, yr^{-1}}$)"
    axs[2, 4].ylabel = L"Modelled $u_\mathrm{s}$ ($\mathrm{m \, yr^{-1}}$)"
    axs[2, 4].xticks = (0:3, latexify.([1, 10, 100, 1000]))
    axs[2, 4].yticks = (0:3, latexify.([1, 10, 100, 1000]))
    xlims!(axs[2, 4], (-1, 3.5))
    ylims!(axs[2, 4], (-1, 3.5))
    lines!(axs[2, 4], -1:0.5:3.5, -1:0.5:3.5, color = :green, linewidth = 3)
    axs[2, 4].yaxisposition = :right

    Colorbar(fig[0, 3], vertical = false, width = Relative(rw), flipaxis = true,
        label = L"$u_\mathrm{s}$ ($\mathrm{m \, yr^{-1}}$)",
        ticks = (0:3, latexify.([1, 10, 100, 1000])); cmaps["uxy_s"]...)
    Colorbar(fig[0, 4], vertical = false, width = Relative(rw), flipaxis = true,
        label = L"$z_\mathrm{b}$ (km) $\,$", ticks = latexifyticks(-4:2, 1f3);
        cmaps["z_bed"]...)
    Colorbar(fig[3, 3], vertical = false, width = Relative(rw + 0.2), flipaxis = false,
        label = L"$u_\mathrm{s}$ error at PD ($\mathrm{m \, yr^{-1}}$)",
        ticks = (-3:3, latexify.([-1000, -100, -10, 0, 10, 100, 1000])); diff_cmap...)

    px, py = -2900, 2450
    text!(axs[1, 3], px, py, text = L"\textbf{(c)}", color = :white, fontsize = fs + fs_off)
    text!(axs[1, 4], px, py, text = L"\textbf{(d)}", color = :white, fontsize = fs + fs_off)
    text!(axs[2, 4], -0.9, 3.0, text = L"\textbf{(h)}", color = :gray5, fontsize = fs + fs_off)
    text!(axs[2, 4], 0.5, -0.8, color = :gray5, fontsize = fs,
        text = "RMSE=$(round(rmse; digits = 2)) m/yr")
    text!(axs[2, 3], px, py, text = L"\textbf{(g)}", color = :gray5, fontsize = fs + fs_off)
    rowgap!(fig.layout, 1, 10)
    rowgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 3, -50)
    colgap!(fig.layout, 5)

    save("$target_dir/figA2.png", fig)
    save("$target_dir/figA2.pdf", fig)
end
plot_error_z(filepaths[1], 1, dx, target_dir)