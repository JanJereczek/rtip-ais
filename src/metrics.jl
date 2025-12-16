function get_area_perimeter_ratio(aqef)
    apr = [Vector{Float32}(undef, length(t)) for t in aqef.t_2D]
    for (i, xp) in enumerate(aqef.xps)
        mask = Bool.(round.(ncread(joinpath([xp, "0/yelmo2D.nc"]), "f_ice")))
        marg = similar(mask)
        for k in axes(mask, 3)
            for i in axes(mask, 1)[2:end-1]
                for j in axes(mask, 2)[2:end-1]
                    marg[i, j, k] = mask[i, j, k] && (
                        mask[i-1, j, k] < 1 || mask[i+1, j, k] < 1 ||
                        mask[i, j-1, k] < 1 || mask[i, j+1, k] < 1)
                end
            end
        end
        apr[i] .= [count(view(mask, :, :, k)) / count(view(marg, :, :, k)) for
            k in axes(mask, 3)]
    end
    return apr
end


function get_slope_flow_product(
    xp;
    tol = 10, # min velocity
    dslp = 3,
    dx = 16f3,  # grid spacing in meters
    )

    # grline or lakes
    grline = (0.5 .< ncread(joinpath([xp, "0/yelmo2D.nc"]), "mask_ocn") .< 1.5)
        # .|| (4.5 .< ncread(joinpath([xp, "0/yelmo2D.nc"]), "mask_ocn") .< 5.5)
    ux = ncread(joinpath([xp, "0/yelmo2D.nc"]), "ux_s")
    uy = ncread(joinpath([xp, "0/yelmo2D.nc"]), "uy_s")
    zb = ncread(joinpath([xp, "0/yelmo2D.nc"]), "z_bed")
    flslp = similar(zb)

    for k in axes(grline, 3)
        for i in axes(grline, 1)[dslp+1:end-dslp]
            for j in axes(grline, 2)[dslp+1:end-dslp]
                if (grline[i,j,k] && (abs(ux[i,j,k]) + abs(uy[i,j,k]) > tol))
                    # flslp[i, j, k] = (zb[i+dslp, j, k] - zb[i-dslp, j, k]) * ux[i,j,k] +
                    #     (zb[i, j+dslp, k] - zb[i, j-dslp, k]) * uy[i,j,k]
                    zb_x = 0
                    zb_y = 0
                    for dij in 1:dslp
                        zb_x += zb[i+dij, j, k] - zb[i-dij, j, k]
                        zb_y += zb[i, j+dij, k] - zb[i, j-dij, k]
                    end
                    flslp[i, j, k] = (zb_x * ux[i,j,k] + zb_y * uy[i,j,k]) /
                        (2 * dslp * dx)   # account for grid spacing
                end
            end
        end
    end
    return grline, flslp
end


function Loess.loess(xout, x, y; span = 0.1)
    model = loess(x, y, span = span)
    return predict(model, xout)
end

function rollmean(x, y; nx = 2)
    yout = fill(NaN, length(x))
    for i in nx+1:length(x)-nx
        yout[i] = mean(y[i-nx:i+nx])
    end
    return yout
end

function sym_mean_diff(x, n, dt)
    y = fill(NaN, length(x))
    for i in n+1:length(x)-n
        y[i] = 0
        for j in 1:n
            y[i] = (x[i+j] - x[i-j]) / (2*j*dt)
        end
    end
    return y
end

function sym_diff(x, n, dt)
    y = fill(NaN, length(x))
    for i in n+1:length(x)-n
        y[i] = (x[i+n] - x[i-n]) / (2*n*dt)
    end
    return y
end