function plot_diffheatmaps(xp1, xp2, prefix; symmetric = true, tol = 1e-8)
    ds1 = NCDataset(xp1)
    ds2 = NCDataset(xp2)
    # iterate over all variables
    for (varname,var) in ds1
        if ndims(var) == 3
            fig = Figure(size = (900,800))
            ax = Axis(fig[1,1], aspect = DataAspect())
            diff = ds2[varname][:, :, 1] .- ds1[varname][:, :, end]
            amp = maximum(abs.(diff))
            if symmetric
                cmap = (colormap = cgrad([:cornflowerblue, :white, :darkorange]),
                    colorrange = (-amp-tol, amp+tol))
            else
                cmap = (colormap = :viridis, colorrange = extrema(diff))
            end
            hm = heatmap!(ax, diff; cmap...)
            Colorbar(fig[1,2], hm, label = "Δ $varname", height = Relative(0.5))
            save(plotsdir("diffheatmap/$prefix/$(varname).png"), fig)
        end
    end
    close(ds1)
    close(ds2)
end

function plot_diffheatmap(xp1, xp2, varname::String, prefix, crange)

    fig = Figure(size = (450, 400))
    ax = Axis(fig[1,1], aspect = DataAspect())
    hidedecorations!(ax)
    hm = plot_diffheatmap!(ax, xp1, xp2, varname, crange)
    Colorbar(fig[1,2], hm, label = "Δ $varname", height = Relative(0.5))

    isdir(plotsdir("diffheatmap/$prefix")) || mkdir(plotsdir("diffheatmap/$prefix"))
    save(plotsdir("diffheatmap/$prefix/$(varname).png"), fig)
end

function plot_diffheatmap!(ax, xp1, xp2, varname::String, crange)
    ds1 = NCDataset(xp1)
    ds2 = NCDataset(xp2)

    diff = ds2[varname][:, :, 1] .- ds1[varname][:, :, 1]
    cmap = (colormap = cgrad([:cornflowerblue, :white, :darkorange]), colorrange = crange)

    close(ds1)
    close(ds2) 

    return heatmap!(ax, diff; cmap...)
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