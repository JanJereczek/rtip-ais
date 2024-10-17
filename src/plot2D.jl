function plot_diffheatmaps(xp1, xp2, prefix)
    ds1 = NCDataset(xp1)
    ds2 = NCDataset(xp2)
    # iterate over all variables
    for (varname,var) in ds1
        if ndims(var) == 3
            fig = Figure(size = (900,800))
            ax = Axis(fig[1,1], aspect = DataAspect())
            hm = heatmap!(ax, ds2[varname][:, :, 1] .- ds1[varname][:, :, end], colormap = :viridis)
            Colorbar(fig[1,2], hm, label = "Δ $varname", height = Relative(0.5))
            save(plotsdir("diffheatmap/$prefix/$(varname).png"), fig)
        end
    end
    close(ds1)
    close(ds2)
end

function plot_heatmaps(xp, path_to_folder, idx)
    ds = NCDataset(xp)
    # iterate over all variables
    for (varname,var) in ds
        if ndims(var) == 3
            fig = Figure(size = (900,800))
            ax = Axis(fig[1,1], aspect = DataAspect())
            hm = heatmap!(ax, ds[varname][idx], colormap = :viridis)
            Colorbar(fig[1,2], hm, label = "Δ $varname", height = Relative(0.5))
            save(plotsdir("$path_to_folder/$varname.png"), fig)
        end
    end
    close(ds)
end