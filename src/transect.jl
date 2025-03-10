
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
    bbloclip = min.(1.1*zlolim, bb)
    bbhiclip = max.(1.1*zlolim, bb)
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
    vlines!(ax, [rr[r_grline]], color = :black, linestyle = :dash, linewidth = 4)
    rsl = round(z_sl_cross[r_grline] - z_bed_cross[r_grline], digits = 1)
    text!(ax, rr[r_grline], z_bed_cross[r_grline], text = "Relative sea level = $rsl m",
        offset = (-300, -100), font = :bold)
end
