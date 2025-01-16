include("intro.jl")

dir = datadir("output/ais/ramps/isos1245-steps")
params = readdlm(joinpath(dir, "params.txt"))
fmax_idx = findfirst(params[1, :] .== "hyster.f_max")
visc = findfirst(params[1, :] .== "isos.viscosity_scaling_method")

T = Float32
fmax = T.(params[2:end, fmax_idx])
visc = T.(params[2:end, visc])
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
hm1 = heatmap!(axs[1], visc, fmax, A_wais, colormap = cmap)
hm2 = heatmap!(axs[2], visc, fmax, V_wais, colormap = cmap)
Colorbar(fig[2, 1], hm1, label = "A_ice", vertical = false, flipaxis = false)
Colorbar(fig[2, 2], hm2, label = "V_sle", vertical = false, flipaxis = false)
fig

# nrows, ncols = 2, 2

# fig = Figure(size=(1600, 980), fontsize = 24)
# ax_forcing = Axis(fig[1, 1])
# ax_state = Axis(fig[2, 1])
# ax_rtip = Axis(fig[:, 2])