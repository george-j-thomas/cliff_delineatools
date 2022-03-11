import numpy as np
import rasterio
import shapely
from shapely.geometry import LineString, shape, Point
from shapely.ops import transform
import pyproj
from shapely.ops import unary_union, nearest_points
import geopandas as gpd
from rasterio import features
import os

def polygonize(dem_path):
    """ Takes raster input and creates polygons for connected pixels of valid raster data.

    Keyword Arguments:
    dem_path -- The path to desired DEM raster.

    Returns:
    src -- Raster read from dem_path by rasterio.
    partials -- List of polygons for valid raster pixels.

    Creates binary raster mask of input DEM (data / no data) and then creates a polygon
        for each connected group of valid pixels.
    """
    src = rasterio.open(dem_path)
    partials = []
    image = src.read(1)
    nodata = -9999
    is_valid = (image != nodata).astype(np.uint8)
    for coords, value in features.shapes(is_valid, transform=src.transform):
        # ignore polygons corresponding to nodata
        if value != 0:
            # convert geojson to shapely geometry
            geom = shape(coords)
            partials.append(geom)
    return src, partials

def wgs2utm(shapely_input,direction="forward"):
    """ Translates shapely objects between WGS 1984 and NAD 1983 UTM Zone 11N 
    coordinate reference systems.

    Keyword Arguments:
    shapely_input -- The shapely object to translate.
    direction -- Either "forward" (WGS->UTM) or "backward" (UTM->WGS). (default="forward")

    Returns:
    utm_output OR wgs_output -- Translated shapely output depending on direction.
    """
    wgs84 = pyproj.CRS('EPSG:4326')
    utm = pyproj.CRS('EPSG:26911')
    if direction == "forward":
        project = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform
        utm_output = transform(project, shapely_input)
        return utm_output
    elif direction == "backward":
        project = pyproj.Transformer.from_crs(utm, wgs84, always_xy=True).transform
        wgs_output = transform(project, shapely_input)
        return wgs_output

def dem_bestfit(dem_path):
    """ Takes raster input and finds linear best fit line through a group of representative points.

    Keyword Arguments:
    dem_path -- The path to desired DEM raster.

    Returns:
    x -- X values of linear best fit line endpoints.
    y -- Y values of linear best fit line endpoints.

    Creates representative points from polygons output by polygonize(), then draws best 
        fit line through them. Line is drawn from x bounds of raster and then extended if 
        line does not cross the full raster y bounds.
    """
    src, partials = polygonize(dem_path)

    rep_points = []
    for poly in partials:   
        rep_points.append(poly.representative_point())
    rep_points_gs = gpd.GeoSeries(rep_points, crs='epsg:26911')
    m,b = np.polyfit(rep_points_gs.x,rep_points_gs.y,1)
    # first = Point(rep_points_gs.x.iloc[0],m*rep_points_gs.x.iloc[0] + b)
    # last = Point(rep_points_gs.x.iloc[-1],m*rep_points_gs.x.iloc[-1] + b)
    x = [src.bounds[0]-200,src.bounds[2]+200]
    y = [m*a + b for a in x]
    # bestline = LineString(list(zip(x,y)))
    first = Point(x[0],y[0])
    last = Point(x[-1],y[-1])
    first_new = first
    last_new = last
    x = [first_new.x,last_new.x]
    y = [first_new.y,last_new.y]
    return x,y

def dem_extents(dem_path,beachdir):
    """ Creates parallel shapely Seaward and Landward extent lines based off bounds and 
    values of DEM.

    Keyword Arguments:
    dem_path -- The path to the DEM file.
    beach_dir -- General cardinal direction that the beach faces.

    Returns:
    Seaward_utm -- Shapely LineString in UTM coordinates of seaward extent of DEM, oriented 
        generally in same direction as alongshore beach direction.
    Landward_utm -- Shapely LineString in UTM coordinates of landward extent of DEM, oriented 
        generally in same direction as alongshore beach direction.

    The Steps are aws follows:
        - polygonize raster into many small polygons with polygonize() (see function description).
        - Create representative point in each polygon and define best fit line through them
            with dem_bestfit() (see function description).
        - Shift best fit line depending on beach orientation and output shifted extent lines.
    """
    
    src, partials = polygonize(dem_path)
    mult = unary_union(partials)
    a = mult.minimum_rotated_rectangle
    edge = a.boundary
    coords = [c for c in edge.coords]
    segments = [shapely.geometry.LineString([a, b]) for a, b in zip(coords,coords[1:])]

    if beachdir == "W":
        Landward_utm = segments[0]
        Seaward_utm = segments[2]
    elif beachdir == "S":
        if segments[0].coords[0][1] > segments[2].coords[1][1]:
            Landward_utm = segments[0]
            Seaward_utm = segments[2]
        else:
            Landward_utm = segments[2]
            Seaward_utm = segments[0]

    return Seaward_utm, Landward_utm

def cliff_delinea_pts(dem_path,beachdir):
    """ Creates point set for Zuzanna's CliffDelineaTool based off bounds and values of DEM.

    Keyword Arguments:
    dem_path -- The path to the DEM file.
    beach_dir -- General cardinal direction that the beach faces.

    Returns:
    Point_df -- A Geodataframe of the points, which can then be exported to a shapefile.

    The Steps are aws follows:
        - Define seaward and landward extent lines from DEM extent and orientation using dem_extents().
        - Create evenly (5 m) spaced points along seaward transect, and find nearest points
        to them along landward transect. These will form DEM transect start/end points.
        - Create evenly (1 m) spaced points along each transect.
        - Extract DEM elevation values for each point 
    """

    Seaward_utm, Landward_utm = dem_extents(dem_path,beachdir)
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
        'PointID':pt_ID,
        'TransectID':final_transectID,
        'Distance':Dist_From_Sea,
        'geometry':final_pts
    }
    Point_df = gpd.GeoDataFrame(Point_Dict,crs='epsg:26911')

    # Extract Raster Values
    src = rasterio.open(dem_path)
    coords = [(x,y) for x, y in zip(Point_df.geometry.x, Point_df.geometry.y)]
    Point_df['Elevation'] = [x for x in src.sample(coords)]
    Point_df['Elevation'] = Point_df.apply(lambda x: x['Elevation'][0], axis=1)
    Point_df['UTM_E'] = Point_df['geometry'].x
    Point_df['UTM_N'] = Point_df['geometry'].y

    # Arrange Columns in order desired by CliffDelineaTool()
    Point_df = Point_df[['PointID','TransectID','Elevation','Distance','UTM_E','UTM_N','geometry']]
    Point_df = Point_df[Point_df.Elevation != -9999]
    return Point_df

dem_fol = r"C:\Users\g4thomas\Documents\CliffDelineation\Files\1_DEMs"
dem_name = "20210324_02619_02669_NoWaves_WillRogers.tif"
beachdir = "S"
dem_path = os.path.join(dem_fol,dem_name)
Point_df = cliff_delinea_pts(dem_path,beachdir)
header = ['PointID','TransectID','Elevation','Distance','UTM_E','UTM_N']
filebase = os.path.basename(dem_path)[:-4]
filename = filebase+'_CliffPoints.txt'
outfolder =  r"C:\Users\g4thomas\Documents\CliffDelineation\Files\2_To_Delineate_txt"
outpath = os.path.join(outfolder,filename)
Point_df.to_csv(outpath,columns = header,index=False)