struct MISIMetrics{V}
    mini::V
    maxi::V
end

function MISIMetrics(t)
    T = eltype(t)
    maxi = Vector{T}(undef, length(t))
    mini = Vector{T}(undef, length(t))
    return MISIMetrics(mini, maxi)
end

function update!(ms::MISIMetrics, field)
    (; mini, maxi) = ms
    for k in axes(field, 3)
        maxi[k] = maximum(view(field, :, :, k))
        mini[k] = minimum(view(field, :, :, k))
    end
end

function get_area_perimeter_ratio(aqef)
    apr = [Vector{Float32}(undef, length(t)) for t in aqef.t_2D]
    for (l, xp) in enumerate(aqef.xps)
        mask = ncread(joinpath([xp, "0/yelmo2D.nc"]), "f_ice") .> 0.5
        grline = (0.5 .< ncread(joinpath([xp, "0/yelmo2D.nc"]), "mask_ocn") .< 1.5)
        apr[l] .= [count(view(mask, :, :, k)) / count(view(grline, :, :, k)) for
            k in axes(mask, 3)]
    end
    return apr
end

function rate_velocity_slope_product(
    xp;
    tol = 1,    # min velocity
    dslp = 1,
    dx = 16f3,  # grid spacing in meters
    T = Float32,
    k_end = nothing,
    time_difference = :central,
    field = "H_grnd",
    )

    grline = (0.5 .< ncread(joinpath([xp, "0/yelmo2D.nc"]), "mask_ocn") .< 1.5)
    t = ncread(joinpath([xp, "0/yelmo2D.nc"]), "time")
    ux = ncread(joinpath([xp, "0/yelmo2D.nc"]), "ux_s")
    uy = ncread(joinpath([xp, "0/yelmo2D.nc"]), "uy_s")
    uxy = sqrt.(ux .^ 2 .+ uy .^ 2)
    X = ncread(joinpath([xp, "0/yelmo2D.nc"]), field)
    zb = ncread(joinpath([xp, "0/yelmo2D.nc"]), "z_bed")
    if isnothing(k_end)
        k_end = length(t) - 1
    end

    X_dt = Matrix{T}(undef, size(ux, 1), size(ux, 2))
    rvsp = zeros(T, size(zb)...)

    for k in 2:k_end
        mod(k, 100) == 0 && @show k
        if time_difference == :backward
            X_dt .= (view(X, :, :, k) .- view(X, :, :, k-1)) ./ (t[k] - t[k-1])
        elseif time_difference == :forward
            X_dt .= (view(X, :, :, k+1) .- view(X, :, :, k)) ./ (t[k+1] - t[k])
        elseif time_difference == :central
            X_dt .= (view(X, :, :, k+1) .- view(X, :, :, k-1)) ./ (t[k+1] - t[k-1])
        else
            error("Invalid time_difference option")
        end

        for i in axes(grline, 1)[dslp+1:end-dslp]
            for j in axes(grline, 2)[dslp+1:end-dslp]
                if (grline[i,j,k] &&
                    zb[i, j, k] < 0 &&
                    (abs(ux[i,j,k]) + abs(uy[i,j,k]) > tol)
                )
                    zb_dx = 0
                    zb_dy = 0
                    for dij in 1:dslp
                        zb_dx += (zb[i+dij, j, k] - zb[i-dij, j, k]) / (2 * dij * dx)
                        zb_dy += (zb[i, j+dij, k] - zb[i, j-dij, k]) / (2 * dij * dx)
                    end
                    zb_dx /= dslp
                    zb_dy /= dslp
                    rvsp[i, j, k] = max(zb_dx * ux[i,j,k] + zb_dy * uy[i,j,k], 0) * X_dt[i,j]
                end
            end
        end
    end
    return grline, rvsp
end


# function get_area_perimeter_ratio(aqef)
#     apr = [Vector{Float32}(undef, length(t)) for t in aqef.t_2D]
#     for (i, xp) in enumerate(aqef.xps)
#         mask = Bool.(round.(ncread(joinpath([xp, "0/yelmo2D.nc"]), "f_ice")))
#         marg = similar(mask)
#         for k in axes(mask, 3)
#             mod(k, 100) == 1 && @show k
#             for i in axes(mask, 1)[2:end-1]
#                 for j in axes(mask, 2)[2:end-1]
#                     marg[i, j, k] = mask[i, j, k] && (
#                         mask[i-1, j, k] < 1 || mask[i+1, j, k] < 1 ||
#                         mask[i, j-1, k] < 1 || mask[i, j+1, k] < 1)
#                 end
#             end
#             apr[i][k] = count(view(mask, :, :, k)) / count(view(marg, :, :, k))
#         end
#     end
#     return apr
# end


function velocity_slope_product(
    xp;
    tol = 1,    # min velocity
    dslp = 3,
    dx = 16f3,  # grid spacing in meters
    T = Float32,
    k_end = nothing,
    )

    # grline or lakes
    grline = (0.5 .< ncread(joinpath([xp, "0/yelmo2D.nc"]), "mask_ocn") .< 1.5)
    ux = ncread(joinpath([xp, "0/yelmo2D.nc"]), "ux_s")
    uy = ncread(joinpath([xp, "0/yelmo2D.nc"]), "uy_s")
    uxy = sqrt.(ux .^ 2 .+ uy .^ 2)
    zb = ncread(joinpath([xp, "0/yelmo2D.nc"]), "z_bed")
    vsp = zeros(T, size(zb)...)

    for k in axes(grline, 3)
        mod(k, 100) == 0 && @show k
        for i in axes(grline, 1)[dslp+1:end-dslp]
            for j in axes(grline, 2)[dslp+1:end-dslp]
                if (grline[i, j, k] && zb[i, j, k] < 0 && uxy[i,j,k] > tol)
                    zb_dx = 0
                    zb_dy = 0
                    for dij in 1:dslp
                        zb_dx += (zb[i+dij, j, k] - zb[i-dij, j, k]) / (2 * dij * dx)
                        zb_dy += (zb[i, j+dij, k] - zb[i, j-dij, k]) / (2 * dij * dx)
                    end
                    zb_dx /= dslp
                    zb_dy /= dslp
                    vsp[i, j, k] = max(zb_dx * ux[i,j,k] + zb_dy * uy[i,j,k], 0)
                end
            end
        end
    end
    return grline, vsp
end


# function Loess.loess(xout, x, y; span = 0.1)
#     model = loess(x, y, span = span)
#     return predict(model, xout)
# end

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


struct FullMISIMetrics{T, V}
    p_stat::T
    p_norm::T
    mini::V
    maxi::V
    normed::V
    lo_percentile::V
    hi_percentile::V
    lo_mean::V
    hi_mean::V
end

function FullMISIMetrics(t, p_stat, p_norm)
    T = eltype(t)
    maxi = Vector{T}(undef, length(t))
    mini = Vector{T}(undef, length(t))
    normed = Vector{T}(undef, length(t))
    lo_percentile = Vector{T}(undef, length(t))
    hi_percentile = Vector{T}(undef, length(t))
    lo_mean = Vector{T}(undef, length(t))
    hi_mean = Vector{T}(undef, length(t))
    return FullMISIMetrics(p_stat, p_norm, mini, maxi, normed,
        lo_percentile, hi_percentile, lo_mean, hi_mean)
end

function update!(ms::FullMISIMetrics, field, mask)
    (; p_stat, p_norm, mini, maxi) = ms
    for k in axes(field, 3)
        maxi[k] = maximum(view(field, :, :, k))
        mini[k] = minimum(view(field, :, :, k))
        # normed[k] = mean(view(field, :, :, k)[view(mask, :, :, k)] .^ p_norm) ^ (1 / p_norm)
        # lo_percentile[k] = percentile(field[:, :, k][view(mask, :, :, k)], p_stat)
        # hi_percentile[k] = percentile(field[:, :, k][view(mask, :, :, k)], 1 - p_stat)
        # lo_mean[k] = mean(field[:, :, k][view(field, :, :, k) .< lo_percentile[k]])
        # hi_mean[k] = mean(field[:, :, k][view(field, :, :, k) .> hi_percentile[k]])
    end
end
