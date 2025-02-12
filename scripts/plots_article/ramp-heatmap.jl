include("../intro.jl")

dir = datadir("output/ais/ramps/16km/coarse-nomax")
params = readdlm(joinpath(dir, "info.txt"))
i_dirs = findfirst(params[1, :] .== "rundir")
i_df_dt_max = findfirst(params[1, :] .== "hyster.df_dt_max")
i_f_max = findfirst(params[1, :] .== "hyster.f_max")
i_visc_method = findfirst(params[1, :] .== "isos.viscosity_scaling_method")
n_files = size(params, 1) - 1

visc_methods = unique(params[2:end, i_visc_method])
n_visc_methods = length(visc_methods)

T = Float32
f = [ T[] for _ in 1:n_visc_methods ]
dfdt = [ T[] for _ in 1:n_visc_methods ]
V = [ T[] for _ in 1:n_visc_methods ]
t_end = [ T[] for _ in 1:n_visc_methods ]
for j in 1:n_files
    file = joinpath(dir, "$(params[j + 1, i_dirs])", "yelmo1D.nc")
    visc_method = params[j + 1, i_visc_method]
    i = findfirst(visc_methods .== visc_method)
    t_e = ncread(file, "time")[end]
    if t_e > 10
        push!(t_end[i], t_e)
        push!(f[i], params[j + 1, i_f_max])
        push!(dfdt[i], params[j + 1, i_df_dt_max])
        V[i] = vcat(V[i], ncread(file, "V_sle")[end])
    end
end

@show length(f[1]), length(f[2])


set_theme!(theme_latexfonts())
fs = 24
fig = Figure(size = (1500, 700), fontsize = fs)
axs = [Axis(fig[2, j], aspect = AxisAspect(1)) for j in 1:n_visc_methods]
cmap = (colormap = cgrad([:darkred, :white, :darkblue]), colorrange = (54, 58))
Colorbar(fig[1, 2], label = L"$V_\mathrm{WAIS}$ (m SLE)", vertical = false,
    width = Relative(0.8), flipaxis = true; cmap...)

logxticks = (-4:0, [L"10^{%$n} $\,$" for n in -4:0])

[axs[j].title = visc_methods[j] for j in 1:n_visc_methods]

# axs[2].yticksvisible = false
# axs[2].yticklabelsvisible = false
axs[1].ylabel = L"$f_\mathrm{max}$ (K)"
axs[n_visc_methods].ylabel = L"$f_\mathrm{max}$ (K)"
axs[n_visc_methods].yaxisposition = :right

for j in eachindex(visc_methods)
    if length(f[j]) > 1
        heatmap!(axs[j], log10.(dfdt[j]), f[j], V[j]; cmap...)
        scatter!(axs[j], log10.(dfdt[j]), f[j], markersize = 10, color = :white)
    end
end
for ax in axs
    ax.xlabel = L"$\dot{f}$ (K/yr)"
    ax.xticks = logxticks
    xlims!(ax, (-4.5, 0.5))
    # ylims!(ax, (2.3, 4.1))
    # hlines!(ax, [3.5], color = :black, linestyle = :dash, linewidth = 6)
end
colgap!(fig.layout, 5)
fig

save(plotsdir("16km/rtip/ramp_heatmap_viscs.png"), fig)