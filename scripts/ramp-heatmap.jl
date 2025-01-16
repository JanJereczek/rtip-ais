include("intro.jl")

dir = datadir("output/ais/ramps/wais-sparse-intermediate")
params = readdlm(joinpath(dir, "params.txt"))
fmax_idx = findfirst(params[1, :] .== "hyster.f_max")
dfdt_idx = findfirst(params[1, :] .== "hyster.df_dt_max")

T = Float32
fmax = T.(params[2:end, fmax_idx])
dfdt = T.(params[2:end, dfdt_idx])
A_wais = T[]
V_wais = T[]

for (root, dirs, files) in walkdir(dir)
    for file in files
        if occursin("yelmo1D_WAIS.nc", file)
            push!(A_wais, ncread(joinpath(root, file), "A_ice")[end])
            push!(V_wais, ncread(joinpath(root, file), "V_sle")[end])
        end
    end
end

cmap = cgrad([:cornflowerblue, :midnightblue])
fig = Figure()
axs = [Axis(fig[1, j]) for j in 1:2]
hm1 = heatmap!(axs[1], log10.(dfdt), fmax, A_wais, colormap = cmap)
hm2 = heatmap!(axs[2], log10.(dfdt), fmax, V_wais, colormap = cmap)
Colorbar(fig[2, 1], hm1, label = "A_ice", vertical = false, flipaxis = false)
Colorbar(fig[2, 2], hm2, label = "V_sle (m)", vertical = false, flipaxis = false)

logxticks = (-4:0, ["10^$n" for n in -4:0])
axs[1].xlabel = "df/dt (K/yr)"
axs[2].xlabel = "df/dt (K/yr)"
axs[1].xticks = logxticks
axs[2].xticks = logxticks
axs[1].ylabel = "f_max (K)"
axs[2].yaxisposition = :right
fig
save(plotsdir("rtip/ramp_heatmap.png"), fig)
# nrows, ncols = 2, 2

# fig = Figure(size=(1600, 980), fontsize = 24)
# ax_forcing = Axis(fig[1, 1])
# ax_state = Axis(fig[2, 1])
# ax_rtip = Axis(fig[:, 2])