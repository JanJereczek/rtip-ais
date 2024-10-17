include("intro.jl")

filepaths = String[]
for (root, dirs, files) in walkdir(datadir("output/ais/spinup"))
    for file in files
        if occursin("yelmo2D.nc", file)
            push!(filepaths, joinpath(root, file))
        end
    end
end

fig = Figure(size = (1600, 900), fontsize = 30)
axs = [Axis(fig[1, j], aspect = DataAspect()) for j in 1:2]
[hidedecorations!(ax) for ax in axs]

H_opts = (colormap = :balance, colorrange = (-1500, 1500))
u_opts = (colormap = cgrad([:orange, :white, :purple]), colorrange = (-1000, 1000))
Colorbar(fig[2, 1], vertical = false, width = Relative(0.5); H_opts...)
Colorbar(fig[2, 2], vertical = false, width = Relative(0.5); u_opts...)

t_e = []
e_u = []
e_H = []
for i in eachindex(filepaths)
    filepath = filepaths[i]
    if occursin("16km", filepath)
        println("Not processing 16km file")
        push!(t_e, Inf32)
        push!(e_H, Inf32)
        push!(e_u, Inf32)
    else
        ds = NCDataset(filepath)
        heatmap!(axs[1], ds["H_ice_pd_err"][:, :, end]; H_opts...)
        heatmap!(axs[2], ds["uxy_s_pd_err"][:, :, end]; u_opts...)
        
        for ax in axs
            contour!(ax, ds["f_grnd"][:, :, 1]; levels = [0.5],
                linewidth = 2, color = :gray)
            contour!(ax, ds["f_grnd"][:, :, end]; levels = [0.5],
                linewidth = 2, color = :black)
        end

        push!(t_e, ds["time"][:][end])
        push!(e_H, mean(abs.(ds["H_ice_pd_err"][:, :, end])))
        push!(e_u, mean(abs.(ds["uxy_s_pd_err"][:, :, end])))
        close(ds)
        save(plotsdir("spinup/pd_err_$i.png"), fig)
    end
end

writedlm(plotsdir("spinup/t_e.csv"), t_e, ',')
writedlm(plotsdir("spinup/filepaths.csv"), filepaths, ',')
writedlm(plotsdir("spinup/e_H.csv"), e_H, ',')
writedlm(plotsdir("spinup/e_u.csv"), e_u, ',')