z_srf_map = cgrad([:gray10, :gray95], 0:0.1:1, categorical = true)
z_bed_map = cgrad(:oleron, [2/3])
z_bed_map2 = cgrad(:oleron, [3/4])
z_bed_map3 = cgrad(:oleron, [1/3])
z_bed_map5 = cgrad(:oleron)

p = 6
nl_bathym = (range(0.0, stop = (2/3)^p, length = 12)) .^ (1/p)
lin_oro = range(2/3 + nl_bathym[end] - nl_bathym[end-1], stop = 1, length = 5)
z_bed_map4 = cgrad(:oleron, vcat(nl_bathym, lin_oro))

# z_bed_map4 = cgrad(:oleron, [2/3])

uxy_s_map = cgrad(:inferno, 0:0.1:1, categorical = true)
# uxy_s_map2 = cgrad([:white, :slateblue4], sqrt.(0:0.1:1), categorical = true)
uxy_s_map2 = cgrad([:white, :darkred], sqrt.(0:0.1:1), categorical = true)
uxy_s_map3 = cgrad([:white, :darkred], range(0, stop = 1, length = 9), categorical = true)
u_bed_map = cgrad([:orange, :white, :midnightblue], [1/9])
z_sl_map = cgrad(:BrBg)
smb_map = cgrad(:PuOr, [6/8])
T_shlf_map = cgrad(:jet)

cmaps = Dict(
    "z_srf" => (colormap = z_srf_map, colorrange = (1e-1, 4e3),
        lowclip = :transparent, highclip = :white),
    "z_bed" => (colormap = z_bed_map, colorrange = (-4e3, 2e3),
        lowclip = z_bed_map[1], highclip = z_bed_map[end]),
    "z_bed2" => (colormap = z_bed_map2, colorrange = (-6e3, 2e3),
        lowclip = z_bed_map2[1], highclip = z_bed_map2[end]),
    "z_bed3" => (colormap = z_bed_map3, colorrange = (-1e3, 2e3),
        lowclip = z_bed_map3[1], highclip = z_bed_map3[end]),
    "z_bed4" => (colormap = z_bed_map4, colorrange = (-4e3, 2e3),
        lowclip = z_bed_map4[1], highclip = z_bed_map4[end]),
    "z_bed5" => (colormap = z_bed_map5, colorrange = (-2.5e3, 2.5e3),
        lowclip = z_bed_map5[1], highclip = z_bed_map[end]),
    "z_bed6" => (colormap = cgrad(z_bed_map5,
        range(0, stop = 1, length = 20), categorical = true),
        colorrange = (-1.6e3, 1.6e3),
        lowclip = z_bed_map5[1], highclip = z_bed_map5[end]),
    "uxy_s" => (colormap = cgrad(:inferno, range(0, stop = 1, length = 11), categorical = true),
        colorrange = (-1, 4), lowclip = uxy_s_map[1], highclip = uxy_s_map[end]),
    "uxy_s2" => (colormap = uxy_s_map2, colorrange = (-1, 4),
        lowclip = uxy_s_map2[1], highclip = :transparent),
    "uxy_s3" => (colormap = uxy_s_map3, colorrange = (0.5, 4.5),
        lowclip = uxy_s_map3[1], highclip = :transparent),
    "u_bed" => (colormap = u_bed_map, colorrange = (-100, 800),
        lowclip = u_bed_map[1], highclip = u_bed_map[end]),
    "z_sl" => (colormap = z_sl_map, colorrange = (-10, 10),
        lowclip = z_sl_map[1], highclip = z_sl_map[end]),
    "smb" => (colormap = smb_map, colorrange = (-60, 20)),
    "T_shlf" => (colormap = T_shlf_map, colorrange = (271, 276)),
    "T_srf" => (colormap = cgrad(:jet), colorrange = (240, 290)),
    "f_grnd" => (levels = [0.5], color = :red, linewidth = 3),
)