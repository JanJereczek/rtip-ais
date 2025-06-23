
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

function crossection(x, ii, jj)
    crossx = zeros(eltype(x), length(ii), size(x, 3))
    for i in eachindex(ii)
        crossx[i, :] = x[ii[i], jj[i], :]
    end
    return vec(crossx)
end

function plot_transect(ax, z_bed, z_bed_ref, f_grnd, H_ice, z_sl, ii, jj, rr)
    # Preprocess fields
    mask_continent = z_bed .> 0
    mask_grounded = f_grnd .>= 1
    mask_ice = H_ice .> 0
    mask_ocean = .!(mask_continent .| mask_grounded)
    mask_shelves = mask_ocean .& mask_ice
    mask_open_ocean = mask_ocean .& .!(mask_ice)

    zsbottom = fill(0f0, size(H_ice))
    zstop = fill(0f0, size(H_ice))
    zsbottom[mask_ocean] .= (z_sl[mask_ocean] .- H_ice[mask_ocean]) .* (931 / 1028)
    zstop[mask_ocean] .= zsbottom[mask_ocean] .+ H_ice[mask_ocean]

    z_bed_cross = crossection(z_bed, ii, jj)
    H_ice_cross = crossection(H_ice, ii, jj)
    mask_grounded_cross = crossection(mask_grounded, ii, jj)
    z_sl_cross = crossection(z_sl, ii, jj)
    zsbottom_cross = crossection(zsbottom, ii, jj)
    zstop_cross = crossection(zstop, ii, jj)
    mask_open_ocean_cross = crossection(mask_open_ocean, ii, jj)

    zz = z_bed_cross + H_ice_cross .* mask_grounded_cross
    bb = z_bed_cross
    gg = z_sl_cross
    ssmask = mask_open_ocean_cross
    zzloclip = min.(zz, bb)
    zzhiclip = max.(zz, bb)
    bbloclip = min.(1.1*zzloclip, bb)
    bbhiclip = max.(1.1*zzloclip, bb)
    sshiclip = max.(ssmask, bb)
    zshelfbottom = zsbottom_cross
    zshelftop = zstop_cross

    lines!(ax, rr, ssmask, color = :dodgerblue4)
    band!(ax, rr, bb, sshiclip, color = :royalblue, label = L"Open ocean $\,$")
    band!(ax, rr, zzloclip, zzhiclip, color = :skyblue1, label = L"Grounded ice $\,$")
    band!(ax, rr, zshelfbottom, zshelftop, color = :lightskyblue1, label = L"Floating ice $\,$")
    lines!(ax, rr, zz, color = :deepskyblue1, linewidth = 3)

    lines!(ax, rr, z_bed_cross, color = :saddlebrown, linewidth = 1, linestyle = :dash)

    lines!(ax, rr, bb, color = :saddlebrown, linewidth = 3)
    band!(ax, rr, bbloclip, bbhiclip, color = :tan1, label = L"Bedrock $\,$")
    if isnothing(z_bed_ref)
        println("not plotting ref bedrock")
    else
        z_bed_2_cross = crossection(z_bed_ref, ii, jj)
        lines!(ax, rr, z_bed_2_cross, color = :saddlebrown, linewidth = 3, linestyle = :dash)
    end
    axislegend(ax, position = :lt, nbanks = 1, framevisible = false)
    # vlines!(ax, [rr[r_grline]], color = :black, linestyle = :dash, linewidth = 4)
    # rsl = round(z_sl_cross[r_grline] - z_bed_cross[r_grline], digits = 1)
    # text!(ax, rr[r_grline], z_bed_cross[r_grline], text = "Relative sea level = $rsl m",
    #     offset = (-300, -100), font = :bold)
end


function plot_simple_transect(ax, z_bed, f_grnd, H_ice, z_sl, ii, jj, rr)
    # Preprocess fields
    mask_continent = z_bed .> z_sl
    mask_grounded = f_grnd .>= 1
    mask_ice = H_ice .> 0
    mask_ocean = .!(mask_continent .| mask_grounded)
    mask_shelves = mask_ocean .& mask_ice
    mask_open_ocean = mask_ocean .& .!(mask_ice)

    zsbottom = fill(0f0, size(H_ice))
    zstop = fill(0f0, size(H_ice))
    zsbottom[mask_ocean] .= (z_sl[mask_ocean] .- H_ice[mask_ocean]) .* (931 / 1028)
    zstop[mask_ocean] .= zsbottom[mask_ocean] .+ H_ice[mask_ocean]

    z_bed_cross = crossection(z_bed, ii, jj)
    H_ice_cross = crossection(H_ice, ii, jj)
    mask_grounded_cross = crossection(mask_grounded, ii, jj)
    z_sl_cross = crossection(z_sl, ii, jj)
    zsbottom_cross = crossection(zsbottom, ii, jj)
    zstop_cross = crossection(zstop, ii, jj)
    mask_open_ocean_cross = crossection(mask_open_ocean, ii, jj)
    lolim = fill(minimum(z_bed_cross)-100, size(z_bed_cross)...)

    zz = z_bed_cross + H_ice_cross .* mask_grounded_cross
    ssmask = mask_open_ocean_cross
    zzloclip = min.(zz, z_bed_cross)
    zzhiclip = max.(zz, z_bed_cross)

    bbloclip = min.(lolim, z_bed_cross)
    bbhiclip = max.(lolim, z_bed_cross)
    sshiclip = max.(ssmask, z_bed_cross)
    zshelfbottom = zsbottom_cross
    zshelftop = zstop_cross

    band!(ax, rr, z_bed_cross, z_sl_cross, color = :royalblue, label = L"Ocean $\,$")
    # band!(ax, rr, zzloclip, zzhiclip, color = :skyblue1, label = L"Grounded ice $\,$")
    band!(ax, rr, zzloclip, zzhiclip, color = :gray80, label = L"Grounded ice $\,$")
    band!(ax, rr, zshelfbottom, zshelftop, color = :lightskyblue1, label = L"Floating ice $\,$")
    # lines!(ax, rr, zz, color = :deepskyblue1, linewidth = 3)
    lines!(ax, rr, zz, color = :gray50, linewidth = 3)
    lines!(ax, rr, z_bed_cross, color = :saddlebrown, linewidth = 1, linestyle = :dash)
    lines!(ax, rr, z_bed_cross, color = :saddlebrown, linewidth = 3)
    band!(ax, rr, bbloclip, bbhiclip, color = :tan1, label = L"Bedrock $\,$")
    lines!(ax, rr, z_sl_cross, color = :dodgerblue4, linestyle = :dash)

    xlims!(ax, (minimum(rr), maximum(rr)))
    ylims!(ax, (minimum(z_bed_cross) - 100, maximum(zz) + 100))
end

function indices_from_mask(mask, X, Y, begin_transect)
    idx = CartesianIndices(mask)
    mask_idx = idx[mask]
    remaining_idx = copy(mask_idx)
    ordered_idx = CartesianIndex[]
    xt, yt, rt = T[], T[], T[0]
    filler = CartesianIndex(100_000, 100_000)

    if begin_transect == :left
        i1 = argmin([idx.I[1] for idx in mask_idx])
    elseif begin_transect == :right
        i1 = argmax([idx.I[1] for idx in mask_idx])
    elseif begin_transect == :top
        i1 = argmax([idx.I[2] for idx in mask_idx])
    elseif begin_transect == :bottom
        i1 = argmin([idx.I[2] for idx in mask_idx])
    end

    push!(ordered_idx, mask_idx[i1])
    push!(xt, X[mask_idx[i1]])
    push!(yt, Y[mask_idx[i1]])
    remaining_idx[i1] = filler

    for i in eachindex(mask_idx)[1:end]
        I = ordered_idx[end]
        diff = (remaining_idx .- I)
        inext = argmin( [x[1]^2 + x[2]^2 for x in diff] )
        if mask_idx[inext] != filler
            push!(ordered_idx, mask_idx[inext])
            remaining_idx[inext] = filler
            push!(xt, X[mask_idx[inext]])
            push!(yt, Y[mask_idx[inext]])
        end
    end

    for i in eachindex(xt)[2:end]
        push!(rt, sqrt( (xt[i] - xt[i-1])^2 + (yt[i] - yt[i-1])^2 ) + rt[end])
    end

    return ordered_idx, rt
end


