using GeoStats
# using PyPlot, PyCall
# @pyimport mpl_toolkits.basemap as basemap
using Proj4
using DataFrames
using DataFramesMeta
using Feather
using Arrow
using Plots; gr(size=(2000,1000))
using Reactive
using Dates
using ArchGDAL
using Images

# Convert values not in right datetime format
Dates.year(dt::Arrow.Datestamp) = Dates.year(convert(Date, dt))

# Load weather data
all_station_data = Feather.read("all_stations.feather")

# Get year from date
all_station_data[:year] = map(d -> year(d), all_station_data[:date])
all_station_data = @where(all_station_data, !ismissing.(:year))
all_station_data[:elevation_] = all_station_data[:altitude_m]

# Convert subdataframe to dataframe
DataFrame(x::SubDataFrame) = DataFrames.DataFrame(Any[view(col, x.rows) for col in x.parent.columns], names(x))

# Define long lat airy coordinates
longlat_airy = Proj4.Projection("+proj=longlat +ellps=airy +towgs84=0,0,0 +no_defs")
# Define Albers Equal Area coordinates
aea = Projection("+proj=aea +ellps=WGS84 +lat_1=51.29427447728234 +lat_2=58.02853263411481 +lon_0=-4.482421875")

longlat_wgs84 = Projection("+proj=longlat +datum=WGS84 +no_defs")
# osgb36 = Projection("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.999601 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894 ")

scale_between(array, min_value, max_value) = (array - minimum(array))/(maximum(array) - minimum(array)) * max_value + min_value

# Get elevation data from GeoTiff
elevation_data, geotransform = ArchGDAL.registerdrivers() do
   ArchGDAL.read("Data/uk_srtm/uk.tif") do dataset
       ArchGDAL.read(dataset, 1),
       ArchGDAL.getgeotransform(dataset)
   end
end

# These are the bounds and pixel sizes for the geotiff in the albers projection
x_min = geotransform[1]
x_res = geotransform[2]
y_max = geotransform[4]
y_res = geotransform[6]

# These helper functions extract the elevation data from the GeoTiff
row_index(y::Number) = Int(fld(y - y_max, y_res)) + 1
column_index(x::Number) = Int(fld(x - x_min, x_res)) + 1

function get_elevation(elevation_data, x::Number, y::Number)
    c = column_index(x)
    r = row_index(y)
    try
        elevation_data[c, r]
    catch
        println((x, y))
    end
end

function get_elevation(elevation_data, x::Array{Number, 1}, y::Array{Number, 1})
    map(d -> get_elevation(elevation_data, d[1], d[2]), zip(x, y))
end


# This is a reactive signal
group_index_signal = Signal(1)

# This function does all the interpolation between stations for a particular year
function interpolate_stations(all_station_data, group_index)
    station_data = @by(@where(all_station_data, :year .== 1941), :station, rain_mm = mean(:rain_mm), lat = first(:lat), long=first(:long), elevation_ = first(:elevation_), tmean_C = mean(:tmean_C))

    coords = Proj4.transform(longlat_wgs84, aea, convert(Array, station_data[[:long, :lat]]))

    station_data[:x] = coords[:, 1]
    station_data[:y] = coords[:, 2]

    station_data[:elevation] = map(d -> get_elevation(elevation_data, d[1], d[2]), zip(station_data[:x], station_data[:y]))

    # return station_data
    station_data[:x] = scale_between(coords[:, 1], 1, size(elevation_data)[1]/10)
    station_data[:y] = scale_between(coords[:, 2], 1, size(elevation_data)[2]/10)

    # Convert station data to a geo data frame
    geodata = GeoDataFrame(dropmissing(station_data), [:x, :y])

    # Domain for kriging interpolation
    domain = RegularGrid{Float64}(div.(size(elevation_data), 10))
    # domain = RegularGrid(collect(div.(size(elevation_data), 10)), [x_min, y_max + size(elevation_data)[2]*y_res], 10*[x_res,-y_res])

    problem = EstimationProblem(geodata, domain, :tmean_C)
    solver = Kriging(:tmean_C => @NT(variogram=GaussianVariogram(range=350.), drifts=[x -> get_elevation(elevation_data, x[1]*10, x[2]*10)]))
    solution = solve(problem, solver)
    digest(solution)

    display(Plots.plot(solution))
    # display(Plots.scatter!(station_data[:x], station_data[:y]))
    solution
end

# This transforms a reactive signal into a resulting graph
map(value -> interpolate_stations(all_station_data, value), group_index_signal)

# for i in 1:200
#        push!(number, i)
#        sleep(2)
# end

coords = interpolate_stations(all_station_data, 1)
coords[:row] = row_index.(coords[:, :y])
coords[:col] = column_index.(coords[:, :x])
# Plots.heatmap(rotl90(elevation_data), aspect_ratio=1)
Plots.heatmap(max(rotl90(imresize(elevation_data, div.(size(elevation_data), 2))), 0), aspect_ratio=1)
# Plots.scatter!(coords[:col], size(elevation_data)[2] - coords[:row])
Plots.scatter!(coords[:col]./2, (size(elevation_data)[2] - coords[:row])./2)
