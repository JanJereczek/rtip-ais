using DrWatson
@quickactivate "rtip-ais"

using CairoMakie
using DelimitedFiles
using Interpolations
using NCDatasets
using NetCDF
using Statistics

# Here you may include files from the source directory
[include(srcdir("$src_file")) for src_file in readdir(srcdir())]

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)
