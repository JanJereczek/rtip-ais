function plot_lines(xp, path_to_folder; stride = 1)
    ds = NCDataset(xp)
    # iterate over all variables
    for (varname,var) in ds
        if ndims(var) == 1
            fig = Figure(size = (900,800))
            ax = Axis(fig[1,1])
            if length(ds[varname][:]) > 1_000_000
                lines!(ax, ds[varname][:][1:stride:end])
            else
                lines!(ax, ds[varname][:])
            end
            save(plotsdir("$path_to_folder/$varname.png"), fig)
        end
    end
    close(ds)
end

function rectangle!(ax, p0, w, h)
    poly!(
        ax,
        Point2f[
            (p0[1], p0[2]),
            (p0[1] + w, p0[2]),
            (p0[1] + w, p0[2] + h),
            (p0[1], p0[2] + h),
        ],
        color = :white,
        strokecolor = :black,
        strokewidth = 1,
    )
end