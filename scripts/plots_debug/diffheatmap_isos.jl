include("../intro.jl")

# Load the data
prefix = datadir("output/test")
dirs = ["isos_20yr", "isos_restart10yr_20yr", "isos_yelmorestart10yr_20yr", "isos_10yr"]
files_isos = ["$prefix/$fn/0/isos_restart.nc" for fn in dirs]
files_1D = ["$prefix/$fn/0/yelmo1D.nc" for fn in dirs]

w_range = (-0.2, 0.2)
we_range = (-2, 2)
dz_ss_range = (-0.2, 0.2)
z_bed_range = (-2, 2)
bsl_range = (-0.1, 0.1)

plot_diffheatmap(files_isos[1], files_isos[2], "w", dirs[2], w_range)
plot_diffheatmap(files_isos[1], files_isos[2], "we", dirs[2], we_range)
plot_diffheatmap(files_isos[1], files_isos[2], "dz_ss", dirs[2], dz_ss_range)
plot_diffheatmap(files_isos[1], files_isos[2], "z_bed", dirs[2], z_bed_range)
plot_diffheatmap(files_isos[1], files_isos[2], "bsl", dirs[2], bsl_range)

plot_diffheatmap(files_isos[1], files_isos[3], "w", dirs[3], w_range)
plot_diffheatmap(files_isos[1], files_isos[3], "we", dirs[3], we_range)
plot_diffheatmap(files_isos[1], files_isos[3], "dz_ss", dirs[3], dz_ss_range)
plot_diffheatmap(files_isos[1], files_isos[3], "z_bed", dirs[3], z_bed_range)
plot_diffheatmap(files_isos[1], files_isos[3], "bsl", dirs[3], bsl_range)

fig = Figure(size = (900,800))
ax = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
u = [load_1Dvar(fn, "uxy_s") for fn in files_1D]
V = [load_1Dvar(fn, "V_sle") for fn in files_1D]

function load_1Dvar(fn, varname)
    ds = NCDataset(fn)
    var = ds[varname][:]
    close(ds)
    return var
end

lw = 3
lines!(ax, u[1], label = dirs[1], linewidth = lw)
lines!(ax, 3:5, u[2], label = dirs[2], linewidth = lw)
lines!(ax, 3:5, u[3], label = dirs[3], linewidth = lw)
lines!(ax, u[4], label = dirs[4], linestyle = :dash, linewidth = lw + 2)
ax.ylabel = "uxy_s"
ax.xticklabelsvisible = false
axislegend(ax, location = :lb)

lines!(ax2, V[1], label = dirs[1], linewidth = lw)
lines!(ax2, 3:5, V[2], label = dirs[2], linewidth = lw)
lines!(ax2, 3:5, V[3], label = dirs[3], linewidth = lw)
lines!(ax2, V[4], label = dirs[4], linestyle = :dash, linewidth = lw + 2)
ax2.xlabel = "time"
ax2.ylabel = "V_sle"
axislegend(ax2, location = :lb)

fig