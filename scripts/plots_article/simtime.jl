rsb_rtip_xps = 77 * 5
rsb_rtip_xps_time = 45f3
rsb_time = rsb_rtip_xps_time * rsb_rtip_xps

wais_rtip_xps = 7 * 8 * 5
wais_rtip_xps_time = 30f3
wais_time = wais_rtip_xps_time * wais_rtip_xps

hyster_time = 500f3

step_time = 360 * 30f3

time = sum([rsb_time, wais_time, hyster_time, step_time])
@show time