include("../intro.jl")

pfx = datadir("output/ais/ramps/")
cases = ["wais-sparse-nomvisc", "wais-sparse-hivisc", "wais-sparse-maxvisc"]
dirs = [joinpath(pfx, case) for case in cases]
ndirs = length(dirs)

set_theme!(theme_latexfonts())
fs = 24
fig = Figure(size = (1500, 700), fontsize = fs)
axs = [Axis(fig[2, j], aspect = AxisAspect(1)) for j in 1:ndirs]
cmap = (colormap = cgrad([:darkred, :white, :darkblue]), colorrange = (0.5, 4.5))
Colorbar(fig[1, 2], label = L"$V_\mathrm{WAIS}$ (m SLE)", vertical = false,
    width = Relative(0.8), flipaxis = true; cmap...)

logxticks = (-4:0, [L"10^{%$n} $\,$" for n in -4:0])

axs[1].title = "Nominal viscosity"
axs[2].title = "High viscosity (10*Nominal)"
axs[3].title = "Max viscosity"

axs[2].yticksvisible = false
axs[2].yticklabelsvisible = false
axs[1].ylabel = L"$f_\mathrm{max}$ (K)"
axs[3].ylabel = L"$f_\mathrm{max}$ (K)"
axs[3].yaxisposition = :right
fig

for j in eachindex(dirs)
    dir = dirs[j]
    params = readdlm(joinpath(dir, "params.txt"))
    fmax_idx = findfirst(params[1, :] .== "hyster.f_max")
    dfdt_idx = findfirst(params[1, :] .== "hyster.df_dt_max")

    T = Float32
    fmax = T.(params[2:end, fmax_idx])
    dfdt = T.(params[2:end, dfdt_idx])
    if occursin("wais-sparse-nomvisc", dir)
        params = readdlm(joinpath(dir, "wais-sparse-high", "params.txt"))
        fmax_idx = findfirst(params[1, :] .== "hyster.f_max")
        dfdt_idx = findfirst(params[1, :] .== "hyster.df_dt_max")

        fmax = vcat(fmax, T.(params[2:end, fmax_idx]))
        dfdt = vcat(dfdt, T.(params[2:end, dfdt_idx]))
    end

    V_wais = T[]

    for (root, dirs, files) in walkdir(dir)
        for file in files
            if occursin("yelmo1D_WAIS.nc", file)
                push!(V_wais, ncread(joinpath(root, file), "V_sle")[end])
            end
        end
    end
    heatmap!(axs[j], log10.(dfdt), fmax, V_wais; cmap...)
    scatter!(axs[j], log10.(dfdt), fmax, markersize = 10, color = :white)
end

for ax in axs
    ax.xlabel = L"$\frac{\mathrm{d} f}{\mathrm{d}t}$ (K/yr)"
    ax.xticks = logxticks
    xlims!(ax, (-4.5, 0.5))
    ylims!(ax, (2.3, 4.1))
    hlines!(ax, [3.5], color = :black, linestyle = :dash, linewidth = 6)
end
colgap!(fig.layout, 5)
save(plotsdir("rtip/ramp_heatmap_viscs.png"), fig)