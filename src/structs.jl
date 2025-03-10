
struct AQEFResults{T<:AbstractFloat}
    dir::String
    xps::Vector{String}
    n_xps::Int
    t_1D::Vector{Vector{T}}
    f::Vector{Vector{T}}
    V_sle::Vector{Vector{T}}
    t_2D::Vector{Vector{T}}
    bmb::Vector{Vector{T}}
    smb::Vector{Vector{T}}
    cmb::Vector{Vector{T}}
end

function AQEFResults(T, dir, xps::Vector{String})
    n_xps = length(xps)
    
    t_1D = Vector{Vector{T}}(undef, n_xps)
    f = Vector{Vector{T}}(undef, n_xps)
    V_sle = Vector{Vector{T}}(undef, n_xps)

    t_2D = Vector{Vector{T}}(undef, n_xps)
    bmb = Vector{Vector{T}}(undef, n_xps)
    smb = Vector{Vector{T}}(undef, n_xps)
    cmb = Vector{Vector{T}}(undef, n_xps)

    for (i, xp) in enumerate(xps)
        xpdir = joinpath([dir, xp, "0"])

        t_1D[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "time")
        f[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "hyst_f_now")
        V_sle[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "V_sle")

        t_2D[i] = ncread(joinpath(xpdir, "yelmo2D.nc"), "time")
        bmb[i] = vec(mean(ncread(joinpath(xpdir, "yelmo2D.nc"), "bmb"), dims = (1, 2)))
        smb[i] = vec(mean(ncread(joinpath(xpdir, "yelmo2D.nc"), "smb"), dims = (1, 2)))
        cmb[i] = vec(mean(ncread(joinpath(xpdir, "yelmo2D.nc"), "cmb"), dims = (1, 2)))
    end
    AQEFResults(dir, xps, n_xps, t_1D, f, V_sle, t_2D, bmb, smb, cmb)
end

struct EquilResults{T<:AbstractFloat}
    dir::String
    files1D::Vector{String}
    files2D::Vector{String}
    f::Vector{T}
    V_sle::Vector{T}
    nt::Vector{Int}
end

function EquilResults(T, dir; level = 3)
    files1D = recursive_global(dir, "yelmo1D.nc", level)
    files2D = recursive_global(dir, "yelmo2D.nc", level)
    n_xps = length(files1D)
    f = zeros(T, n_xps)
    V_sle = zeros(T, n_xps)
    nt = zeros(Int, n_xps)
    for i in eachindex(files1D)
        f[i] = ncread(files1D[i], "hyst_f_now")[end]
        V_sle[i] = ncread(files1D[i], "V_sle")[end]
        nt[i] = length(ncread(files2D[i], "time"))
    end
    idx = sortperm(f)
    return EquilResults(dir, files1D[idx], files2D[idx], f[idx], V_sle[idx], nt[idx])
end

function hyst!(ax, eql::EquilResults; kwargs...)
    lines!(ax, eql.f, eql.V_sle; kwargs...)
    return nothing
end

function hm!(ax, eql::EquilResults, f; kwargs...)
    i = findfirst(eql.f .>= f)
    z_bed, z_srf, f_grnd, f_ice = load_netcdf_2D(eql.files2D[i],
        ["z_bed", "z_srf", "f_grnd", "f_ice"], eql.nt[i])
    heatmap!(ax, z_bed; kwargs...)
    heatmap!(ax, z_srf .* f_ice; kwargs...)
    contour!(ax, f_grnd .+ f_ice, levels = [1.9], color = :red, linewidth = 2)
    return nothing
end

struct Transect
    i1::Int
    i2::Int
    j1::Int
    j2::Int
end
makeline(ax, tr::Transect; kwargs...) = lines!(ax, [tr.i1, tr.i2], [tr.j1, tr.j2]; kwargs...)


struct PrecomputedTransect{T<:AbstractFloat}
    tr::Transect
    x::Vector{T}
    y::Vector{T}
    dist::Vector{T}
end

function line_on(tr::Nothing, X, Y)
    return nothing
end

function line_on(tr::Transect, X, Y)

    ii = [tr.i1, tr.i2]
    jj = [tr.j1, tr.j2]
    i1 = minimum(ii)
    i2 = maximum(ii)
    j1 = minimum(jj)
    j2 = maximum(jj)

    x1, x2 = X[tr.i1, tr.j1], X[tr.i2, tr.j2]
    y1, y2 = Y[tr.i1, tr.j1], Y[tr.i2, tr.j2]

    if abs(tr.i1 - tr.i2) > abs(tr.j1 - tr.j2)
        x = X[i1:i2, 1]
        m = (y2 - y1) / (x2 - x1)
        p = y1 - m * x1
        y = m .* x .+ p
    else
        y = Y[1, j1:j2]
        m = (x2 - x1) / (y2 - y1)
        p = x1 - m * y1
        x = m .* y .+ p
    end

    dist = sqrt.((x .- X[tr.i1, tr.j1]).^2 .+ (y .- Y[tr.i1, tr.j1]).^2)
    idx = sortperm(dist)
    
    return PrecomputedTransect(tr, x[idx], y[idx], dist[idx])
end