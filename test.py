# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext_format_version: '1.1'
#   kernelspec:
#     display_name: Julia 0.6.4
#     language: julia
#     name: julia-0.6
#   language_info:
#     file_extension: .jl
#     mimetype: application/julia
#     name: julia
#     version: 0.6.4
#   toc:
#     nav_menu: {}
#     number_sections: true
#     sideBar: false
#     skip_h1_title: true
#     toc_cell: false
#     toc_position:
#       height: 40px
#       left: 1459.97px
#       right: 20px
#       top: 52.5962px
#       width: 181.571px
#     toc_section_display: none
#     toc_window_display: true
# ---

# + {"language": "julia", "active": "ipynb"}
# include("kriging.jl")
# gr(size=(1.2, 0.7).*1000)

# + {"language": "julia", "active": "ipynb"}
# using Interact

# + {"scrolled": false, "language": "julia", "active": "ipynb"}
# @manipulate for year_value=1941:2017
#     interpolate_stations(all_station_data, year_value)
# end

# + {"language": "julia", "active": "ipynb"}
# interpolate_stations(all_station_data, 1990)

# + {"language": "julia", "active": "ipynb"}
# using GeoJSON
# using JSON

# + {"language": "julia", "active": "ipynb"}
# uk_coast = JSON.parsefile("Data/uk.json")

# + {"language": "julia", "active": "ipynb"}
# uk_coast = GeoJSON.parse(uk_coast)

# + {"language": "julia", "active": "ipynb"}
# print(dict2geo(uk_coast))

# + {"language": "julia", "active": "ipynb"}
# using Shapefile

# + {"scrolled": true, "language": "julia", "active": "ipynb"}
# coastline = read(open("Data/shapes/subunits.shp"), Shapefile.Handle)

# + {"language": "julia", "active": "ipynb"}
# coastlines

# + {"language": "julia", "active": "ipynb"}
# GeoInterface.coordinates(coastlines.shapes)

# + {"language": "julia", "active": "ipynb"}
# station_data = get_station_data_for_year(all_station_data, 1990)

# + {"language": "julia", "active": "ipynb"}
# @manipulate for nbins=1:50
#     γ = EmpiricalVariogram(GeoDataFrame(station_data, [:x, :y]), :tmean_C, nbins=nbins)
#     plot(γ)
# end

# + {"language": "julia", "active": "ipynb"}
# γemp = EmpiricalVariogram(GeoDataFrame(station_data, [:x, :y]), :tmean_C, nbins=50)

# + {"language": "julia", "active": "ipynb"}
# γtheo

# + {"language": "julia", "active": "ipynb"}
# γtheo = fit(SphericalVariogram, γemp)
#
# plot(γemp, label="empirical")
# plot!(γtheo, maxlag=1e6, label="theoretical")

# + {"language": "julia", "active": "ipynb"}
# @manipulate for s=0.0:0.1/50:.1,
#                 r=0.0:2:100,
#                 n=0.0:0.1/50:.1
#     # theoretical variogram
#     γtheo = SphericalVariogram(sill=s, range=r, nugget=n)
#
#     plot(γ)
#     plot!(γtheo, maxlag=1e6, label="theoretical")
# end

# + {"language": "julia", "active": "ipynb"}
#
