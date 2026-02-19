
visc_cases = [
    "m2stddev",
    "m1stddev",
    "nominal",
    "p1stddev",
    "p2stddev",
]
visc_labels = [
    L"$-2 \, \sigma$",
    L"$-1 \, \sigma$",
    L"$0 \, \sigma$",
    L"$+1 \, \sigma$",
    L"$+2 \, \sigma$",
]
visc_num_labels = [
    L"$\textbf{a} \quad -2 \, \sigma$",
    L"$\textbf{b} \quad -1 \, \sigma$",
    L"$\textbf{c} \quad 0 \, \sigma$",
    L"$\textbf{d} \quad +1 \, \sigma$",
    L"$\textbf{e} \quad +2 \, \sigma$",
]
f_pd = 1.2
f_to = 0.25
f2020 = 1.2
polar_amplification = 1.8
pa = polar_amplification

# # Compute SSP slopes
# ssps = load_ssps(polar_amplification; wrt = :pd)
# mean_slope(t, x) = [ (x[i] - x[1])/(t[i] - t[1]) for i in 2:length(x) ]
# mean_slope(x) = mean_slope(view(x, :, 1), view(x, :, 2))
# dssps = [mean_slope(ssp) for ssp in ssps]
# ssp_labels = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
# ssp_colors = [:darkblue, :lightblue, :darkorange, :darkred]
# dfdt_min = minimum([minimum(dssp) for dssp in dssps])
# dfdt_max = maximum([maximum(dssp) for dssp in dssps])
# dfdt_min = 10 ^ floor(log10(dfdt_min))
# dfdt_max = 10 ^ ceil(log10(dfdt_max)) # Let's take some margin

ssp126_2100 = 1.8
ssp245_2100 = 2.7
ssp370_2100 = 3.6
ssp585_2100 = 4.4

ssp126_2100_hi = 2.4
ssp245_2100_hi = 3.5
ssp370_2100_hi = 4.6
ssp585_2100_hi = 5.7

# ssp126_2300

ssp_dfdt_min = (ssp126_2100 - f2020) / 0.8
ssp_dfdt_max = (ssp585_2100_hi - f2020) / 0.8