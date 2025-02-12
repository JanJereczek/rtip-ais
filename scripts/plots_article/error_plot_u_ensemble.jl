include("../intro.jl")

dx = 16
ncdir = datadir("output/ais/spinup/$(dx)km/")
filepaths = recursive_global(ncdir, "yelmo2D.nc", 5)
target_dir = plotsdir("errormap_ensemble/$(dx)km/")

function plot_error_u(fn, hash, dx, target_dir)

    isdir(target_dir) || mkdir(target_dir)
    fn_ref = "/p/projects/megarun/ice_data/Antarctica/ANT-$(dx)KM/ANT-$(dx)KM_VEL-R11-2.nc"

    nt = ncread(fn, "time") |> length
    uxy_s_now = ncread(fn, "uxy_s", start = [1, 1, nt], count = [-1, -1, 1])
    uxy_s_ref = ncread(fn_ref, "uxy_srf")
    z_bed_ref = ncread(fn, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])
    z_bed_now = ncread(fn, "z_bed", start = [1, 1, nt], count = [-1, -1, 1])
    z_srf_ref = ncread(fn, "z_srf", start = [1, 1, 1], count = [-1, -1, 1])
    diff = view(uxy_s_now, :, :, 1) .- uxy_s_ref
    idx = CartesianIndices(diff)
    mask_now = idx[uxy_s_now .!= 0]
    mask_ref = idx[z_srf_ref .!= 0]
    x = axes(diff, 1)
    y = axes(diff, 2)
    X = repeat(x, 1, length(y))
    Y = repeat(y', length(x), 1)

    tol = 1f-8
    nrows, ncols = 2, 2
    rw = 0.7
    ms = 2
    fs = 18
    fs_off = 2
    fig = Figure(size=(1100, 1200), fontsize = fs)
    axs = [Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
    [hidedecorations!(ax) for ax in [axs[1, 1], axs[1, 2], axs[2, 2]]]
    diff_cmap = (colormap = cgrad([:darkred, :white, :darkblue], 11, categorical = true,
        rev = true), lowclip = :darkblue, highclip = :darkred, colorrange = (-3, 3))

    heatmap!(axs[1, 1], view(z_bed_ref, :, :, 1); cmaps["z_bed"]...)
    heatmap!(axs[1, 1], X[mask_ref], Y[mask_ref], log10.(view(uxy_s_ref, mask_ref, 1) .+ tol);
        cmaps["uxy_s"]...)
    heatmap!(axs[1, 2], view(z_bed_now, :, :, 1); cmaps["z_bed"]...)
    heatmap!(axs[1, 2], X[mask_now], Y[mask_now], log10.(view(uxy_s_now, mask_now, 1) .+ tol);
        cmaps["uxy_s"]...)
    heatmap!(axs[2, 2], sign.(diff) .* log10.(abs.(diff .+ tol)); diff_cmap...)

    scatter!(axs[2, 1], log10.(vec(uxy_s_ref) .+ tol), log10.(vec(uxy_s_now) .+ tol),
        markersize = ms, color = :gray5)
    axs[2, 1].xlabel = L"PD observed ice surface velocity ($\mathrm{m \, yr^{-1}}$)"
    axs[2, 1].ylabel = L"PD modelled ice surface velocity ($\mathrm{m \, yr^{-1}}$)"
    axs[2, 1].xticks = (0:3, latexify.([1, 10, 100, 1000]))
    axs[2, 1].yticks = (0:3, latexify.([1, 10, 100, 1000]))
    xlims!(axs[2, 1], (-1, 3.5))
    ylims!(axs[2, 1], (-1, 3.5))

    Colorbar(fig[0, 1], vertical = false, width = Relative(rw), flipaxis = true,
        label = L"Ice surface velocity ($\mathrm{m \, yr^{-1}}$)",
        ticks = (0:3, latexify.([1, 10, 100, 1000])); cmaps["uxy_s"]...)
    Colorbar(fig[0, 2], vertical = false, width = Relative(rw), flipaxis = true,
        label = L"Bedrock elevation (km) $\,$", ticks = latexifyticks(-4:2, 1f3);
        cmaps["z_bed"]...)
    Colorbar(fig[3, 2], vertical = false, width = Relative(rw), flipaxis = false,
        label = L"Surface velocity error w.r.t. PD ($\mathrm{m \, yr^{-1}}$)",
        ticks = (-3:3, latexify.([-1000, -100, -10, 0, 10, 100, 1000])); diff_cmap...)


    text!(axs[1, 1], 10, 170, text = L"\textbf{(a)}", color = :white, fontsize = fs + fs_off)
    text!(axs[1, 2], 10, 170, text = L"\textbf{(b)}", color = :white, fontsize = fs + fs_off)
    text!(axs[2, 1], -0.8, 3, text = L"\textbf{(c)}", color = :gray5, fontsize = fs + fs_off)
    text!(axs[2, 2], 10, 170, text = L"\textbf{(d)}", color = :gray5, fontsize = fs + fs_off)
    rowgap!(fig.layout, 1, 10)
    rowgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 3, -40)
    colgap!(fig.layout, 5)

    save("$target_dir/u-$hash.png", fig)
end
for (i, fn) in enumerate(filepaths)
    plot_error_u(fn, i, dx, target_dir)
end

writedlm("$target_dir/summary.txt", [eachindex(filepaths) filepaths])