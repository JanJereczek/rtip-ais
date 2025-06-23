using DrWatson
@quickactivate "rtip-ais"

using CairoMakie
using Colors
using DelimitedFiles
using Glob
using Interpolations
using LazyGrids
using Loess
using NCDatasets
using NetCDF
using Random
using Statistics
using StatsBase
using FileIO

# Here you may include files from the source directory
[include(srcdir("$src_file")) for src_file in readdir(srcdir())]

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)
