
import numpy as np
import rasterio
from shapely.geometry import LineString
from shapely.ops import transform
import pandas as pd
import pyproj
from shapely.ops import unary_union, nearest_points
import geopandas as gpd


def MopData(csv_path,MopsRange):
    """ Returns dataframe of MOP lines from Mop.csv given a range.

    Keyword Arguments:
    MopsRange -- A range of numbers corresponding to desired lines of the csv.

    Future Upgrades:
    Will take in region name/code (e.g "San Diego" or corresponding code "D")
        as well as range of actual MOP numbers.
    """
    dat = pd.read_csv(csv_path)
    MopsRegion = dat.loc[MopsRange]
    return MopsRegion

def WSG2UTM(shapely_input):
    wgs84 = pyproj.CRS('EPSG:4326')
    utm = pyproj.CRS('EPSG:26911')
    project = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform
    utm_output = transform(project, shapely_input)
    return utm_output

def CliffDelineaPts(dem_path):
    """ Creates point set for Zuzanna's CliffDelineaTool based off bounds and values of DEM.

    Keyword Arguments:
    dem_path -- The path to the DEM file.

    Returns:
    Point_df -- A Geodataframe of the points, which can then be exported to a shapefile.

    The Steps are aws follows:
        - Define seaward and landward extent lines from DEM bounding box, plus 50m buffer.
        - Create evenly (5 m) spaced points along seaward transect, and find nearest points
        to them along landward transect. These will form DEM transect start/end points.
        - Create evenly (1 m) spaced points along each transect.
        - Extract DEM elevation values for each point 
    """
    # Create parallel Seaward/Landward extent lines
    cliffbackX = []
    cliffbackY = []
    offshoreX = []
    offshoreY = []

    # Base Seaward/Landward extents off DEM bounding box
    src = rasterio.open(dem_path)
    buffer = 50
    offshoreX = [src.bounds[0]-buffer]*2
    offshoreY = [src.bounds[3]+buffer,src.bounds[1]-buffer]
    cliffbackX = [src.bounds[2]+buffer]*2
    cliffbackY = [src.bounds[3]+buffer,src.bounds[1]-buffer]
    
    # Convert Seaward and Landward extent lines to UTM
    Seaward_utm = LineString(zip(offshoreX,offshoreY))
    Landward_utm = LineString(zip(cliffbackX,cliffbackY))
    # Create points every (distance_delta) along Seaward extent line
    distance_delta = 5
    distances = np.arange(0, Seaward_utm.length, distance_delta)
    points = [Seaward_utm.interpolate(distance) for distance in distances] + [Seaward_utm.boundary.geoms[1]]
    multipoint_sea = list(unary_union(points).geoms)
    # Find nearest point along Landward extent to each new point on Seaward line
    # Then, build transect lines from seaward to landward points
    multipoint_land = []
    transect_ID = []
    transect_line = []
    transect_num = 0
    for seapoint in multipoint_sea:
        transect_ID.append(transect_num)
        p1,_ = nearest_points(Landward_utm,seapoint)
        multipoint_land.append(p1)
        transect_line.append(LineString([seapoint,p1]))
        transect_num = transect_num + 1
    # Create point every 1m (DEM Resolution) along transect lines
    distance_delta = 1
    final_pts = []
    final_transectID = []
    Dist_From_Sea = []
    for tID in range(len(transect_line)):
        distances = np.arange(0, transect_line[tID].length, distance_delta)
        dist_list = list(distances)
        dist_list.append(dist_list[-1]+1)
        points = [transect_line[tID].interpolate(distance) for distance in distances] + [transect_line[tID].boundary.geoms[1]]
        transect_pts = list(unary_union(points).geoms)
        tID_templist = [transect_ID[tID]]*len(transect_pts)
        final_pts.extend(transect_pts)
        final_transectID.extend(tID_templist)
        Dist_From_Sea.extend(dist_list)
    # Store Point data in shapefile
    pt_ID = list(range(len(final_pts)))
    Point_Dict = {
        'Point_ID':pt_ID,
        'Transect_ID':final_transectID,
        'Dist_Sea':Dist_From_Sea,
        'geometry':final_pts
    }
    Point_df = gpd.GeoDataFrame(Point_Dict,crs='epsg:26911')

    # Extract Raster Values
    src = rasterio.open(dem_path)
    coords = [(x,y) for x, y in zip(Point_df.geometry.x, Point_df.geometry.y)]
    Point_df['Elevation'] = [x for x in src.sample(coords)]
    Point_df['Elevation'] = Point_df.apply(lambda x: x['Elevation'][0], axis=1)

    return Point_df


dem_path = r"C:\Users\g4thomas\Documents\CliffDelineation\20211103_00518_00568_NoWaves_Blacks_beach_cliff_ground.tif"
csv_path = r"C:\Users\g4thomas\Documents\CliffDelineation\Mop.csv"
SDMops = range(931)
Blacks = range(516,571)
MopsRange = Blacks
MopsRegion = MopData(csv_path, MopsRange)
Point_df = CliffDelineaPts(dem_path)

Point_df.to_file('CliffDelineaPts.shp')