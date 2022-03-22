"""
Call the build_polys() function to use cliff base points from cliff_delineatool.py to construct 
beach- and cliff-size MOP area tiles.

"""

from shapely.geometry import Polygon, Point, LineString, LinearRing
from shapely.validation import make_valid
import pandas as pd
import pyproj
import geopandas as gpd
from shapely.geometry.polygon import LinearRing
from shapely.ops import nearest_points
import os
from shapely.ops import split
from shapely.geometry import Point
from cliff_delineatool import fix_line_shp

def basepoints_to_line(basepoints_path):
    """
    Converts data from csv of cliff base points into a shapefile of a line.

    """
    base_points = pd.read_csv(basepoints_path,header=None)
    coordinates = list(zip(base_points.iloc[:,2],base_points.iloc[:,3]))
    line = LineString(coordinates)
    gdr = gpd.GeoSeries(line,crs='epsg:26911')

    folder, name = os.path.split(basepoints_path)
    outpath = os.path.join(folder,name[:-4]+'.shp')
    gdr.to_file(outpath)
    print('The path to the cliff base line shapefile is ',outpath)
    return outpath

def mop_data(mop_csv_path,cliffline_path):
    """ 
    Returns dataframe of MOP lines from Mop.csv given the extent of a cliff line.

    Keyword Arguments:
    mop_csv_path -- Path to the csv containing MOP line data.
    cliffline_path -- Path to shapefile of cliffbase line created by CliffDelineaTool.

    Future Upgrades:
    Will take in region name/code (e.g "San Diego" or corresponding code "D")
        as well as range of actual MOP numbers.
    """

    # Read in csv MOP data
    dat = pd.read_csv(mop_csv_path)
    # Get just backbeach points and and create gdf
    dat_gpd = gpd.GeoDataFrame(dat,geometry=gpd.points_from_xy(dat['BackXutm'],dat['BackYutm']))
    backpts = dat_gpd.geometry.unary_union
    cliffline = gpd.read_file(cliffline_path)

    # If gaps exist, fill them
    if len(cliffline) > 1:
        full_line = fix_line_shp(cliffline)
    else:
        full_line = cliffline.geometry[0]
    full_line = gpd.GeoSeries([full_line])

    # Find mop tiles from extent of cliffline
    first,last = full_line.geometry.boundary[0].geoms
    _, nearest_geom = nearest_points(first, backpts)
    first_i = int(dat_gpd.loc[dat_gpd.geometry== nearest_geom].index[0])
    _, nearest_geom = nearest_points(last, backpts)
    last_i = int(dat_gpd.loc[dat_gpd.geometry== nearest_geom].index[0])

    # Add a few extra mops at the end
    if first_i > last_i:
        MopsRange = range(last_i, first_i+3)
        getrange = range(last_i-1, first_i+3)
    elif first_i == 0:
        MopsRange = range(first_i,last_i+3)
        getrange = range(first_i,last_i+3)
    else:
        MopsRange = range(first_i,last_i+3)
        getrange = range(first_i-1,last_i+3)
            
    MopsRegion = dat.loc[getrange]
    MopsRegion['MOP_index'] = list(getrange)
    return MopsRegion, MopsRange

def get_midpoints(MopsRegion,MopsRange):
    """ 
    Returns midpoints between offshore and backbeach MOP points.

    Keyword Arguments:
    MopsRegion -- The dataframe returned by the mop_data function.

    Returns:
    backmidX -- Longitudes of backbeach midpoints.
    backmidY -- Latitudes of backbeach midpoints.
    offmidX -- Longitudes of offshore midpoints.
    offmidY -- Latitudes of offshore midpoints.

    Midpoint creation is as follows for MOP Area 'x': 
        On backbeach, finds the midpoint between backbeach point 'x-1' and
            backbeach point 'x', which is then assigned to MOP line 'x' and will 
            exist along the south end of MOP Area 'x'
        The same process is performed for offshore points.
    """
    offmidX = []
    offmidY = []
    backmidX = []
    backmidY = []
    # mops_range = range(len(MopsRegion)) 
    for idx in MopsRange:
        geodesic = pyproj.Geod(ellps='WGS84')
        if idx == 0:
            # For first section there is no 'previous midpoint', so just find
            # points 50m south of backbeach & offshore points
            fwd_az,_,_ = geodesic.inv(MopsRegion['BackLon'][idx+1], MopsRegion['BackLat'][idx+1], 
                                    MopsRegion['BackLon'][idx], MopsRegion['BackLat'][idx])
            mid_dist = 50
            startpointX, startpointY,_ = geodesic.fwd(MopsRegion['BackLon'][idx], 
                                                    MopsRegion['BackLat'][idx],fwd_az,mid_dist)
            backmidX.append(startpointX)
            backmidY.append(startpointY)
            startpointX, startpointY,_ = geodesic.fwd(MopsRegion['OffLon'][idx], 
                                                    MopsRegion['OffLat'][idx],fwd_az,mid_dist)
            offmidX.append(startpointX)
            offmidY.append(startpointY)
        elif idx > 0: #and idx < MopsRange[-1]:
            # Previously described process for areas 2 to (End-1)
            back_fwd_az,_,back_dist = geodesic.inv(MopsRegion['BackLon'][idx-1], MopsRegion['BackLat'][idx-1], 
                                                MopsRegion['BackLon'][idx], MopsRegion['BackLat'][idx])
            backtempX,backtempY,_ = geodesic.fwd(MopsRegion['BackLon'][idx-1], MopsRegion['BackLat'][idx-1],
                                                back_fwd_az,back_dist/2)

            off_fwd_az,_,off_dist = geodesic.inv(MopsRegion['OffLon'][idx-1], MopsRegion['OffLat'][idx-1], 
                                                MopsRegion['OffLon'][idx], MopsRegion['OffLat'][idx])
            offtempX,offtempY,_ = geodesic.fwd(MopsRegion['OffLon'][idx-1], MopsRegion['OffLat'][idx-1],
                                            off_fwd_az,off_dist/2)
            backmidX.append(backtempX)
            backmidY.append(backtempY)
            offmidX.append(offtempX)
            offmidY.append(offtempY)

    return backmidX,backmidY,offmidX, offmidY

def get_cliffbacks(backmidX,backmidY,offmidX, offmidY,backdist):
    """ 
    Returns coordinates for onshore vertexes (back of 'Big' tiles),
    given the previously created midpoints.

    Keyword Arguments:
    backmidX -- Longitudes of backbeach midpoints.
    backmidY -- Latitudes of backbeach midpoints.
    offmidX -- Longitudes of offshore midpoints.
    offmidY -- Latitudes of offshore midpoints.
    backdist -- Inshore distance from backbeach midpoints for cliffback vertexes.

    Returns:
    cliffmidX -- Longitudes of cliffback midpoints.
    cliffmidY -- Latitudes of cliffback midpoints.

    Given offshore and backbeach midpoints between MOP transects, finds coordinates of
        'cliff back' points a given distance (backdist) from the backbeach midpoint along the 
        azimuth from the offshore to backbeach midpoint.
    """
    geodesic = pyproj.Geod(ellps='WGS84')
    cliffbackX = []
    cliffbackY = []
    for idx in range(len(backmidX)):
        fwd_az,_,_ = geodesic.inv(offmidX[idx], offmidY[idx],backmidX[idx], backmidY[idx])
        endpointX, endpointY,_ = geodesic.fwd(backmidX[idx], backmidY[idx],fwd_az,backdist)
        cliffbackX.append(endpointX)
        cliffbackY.append(endpointY)

    return cliffbackX,cliffbackY

def validate_geometry(initial_gdf):
    """ 
    Fixes "bowtie" polygons and returns a geodataframe with entirely valid polygons.

    Finds invalid polygons, usually ones with "bowtie" shapes, and runs shapely's make_valid()
        funtion on them. 
    """
    problem_list = list(initial_gdf[~initial_gdf.is_valid]['ID'])
    valid_geom = []
    for _, row in initial_gdf.iterrows():
        if row['ID'] in problem_list:
            t1 = make_valid(row['geometry'])
            t2 = max(t1.geoms, key=lambda a:a.area)
            valid_geom.append(t2)
        else:
            valid_geom.append(row['geometry'])
    initial_gdf['valid_geometry'] = valid_geom
    valid_gdf = initial_gdf.drop(columns=['geometry'])
    valid_gdf = valid_gdf.set_geometry('valid_geometry')

    return valid_gdf

def remove_overlap(valid_gdf):
    """ 
    Takes in a geodataframe containing entirely valid polygons and removes overlaps.

    Overlaps are removed such that among 2+ overlapping polygons, 
        the last polygon in the list will retain the overlapped area.
    """
    fix_geo = []
    for i, row in valid_gdf.iterrows():
        others = valid_gdf[valid_gdf['ID'] != row['ID']]
        diff = row['valid_geometry'].difference(others.unary_union)
        valid_gdf.at[i,'valid_geometry'] = diff
        fix_geo.append(diff)

    return fix_geo

def build_big_tiles(MopsRegion,MopsRange, backdist):
    """ 
    Builds "Big" (Cliff + Beach) MOP area tiles in both WSG 1984 and UTM Zone 11N coordinates.

    Keyword Arguments:
    MopsRegion -- The dataframe returned by the mop_data function.
    backdist -- Inshore distance from backbeach midpoints for cliffback vertexes (see get_cliffbacks function).

    Returns:
    bigtile_df_wsg -- Geodataframe of big tiles in WSG 1984 coordinates.
    bigtile_df_utm -- Geodataframe of big tiles in UTM Zone 11N coordinates.

    First creates backabeach, offshore, and cliffback midpoints 
        (see internal get_midpoints and get_cliffbacks functions).

    Then builds Polygon "x" from vertices in the following order:
        cliffback midpoint (x), 
        backbeach midpoint (x), 
        offshore midpoint (x),
        offshore MOP transect point (x), 
        offshore midpoint (x+1),
        backbeach midpoint (x+1), 
        cliffback midpoint (x+1)

    Then fixes invalid polygons and removes overlaps using internal
        validate_geometry and remove_overlaps functions, respectively (both described below).
    """

    # Get Offshore and Backbeach midpoints
    backmidX,backmidY,offmidX, offmidY = get_midpoints(MopsRegion,MopsRange)

    # Get Cliffback midpoints
    cliffbackX,cliffbackY = get_cliffbacks(backmidX,backmidY,offmidX, offmidY,backdist)

    id = []
    MOP = []
    geometry = []
    # tile_range = range(len(MopsRegion))
    tile_range = MopsRange
    count = 0
    MopNames = MopsRegion['Name'].str.extract('(\d+)', expand=False)
    for idx in tile_range:
        if idx < tile_range[-1]:
            # Builds BIG tiles (Beach + Cliff) to be split later by backbeach line
            vertsX = [cliffbackX[count],backmidX[count],offmidX[count],MopsRegion['OffLon'][idx],
                      offmidX[count+1],backmidX[count+1],cliffbackX[count+1]]#,cliffbackX[idx]]
            vertsY = [cliffbackY[count],backmidY[count],offmidY[count],MopsRegion['OffLat'][idx],
                      offmidY[count+1],backmidY[count+1],cliffbackY[count+1]]#,cliffbackY[idx]]
            r = LinearRing(list(zip(vertsX,vertsY)))
            geo_tmp = Polygon(r)
            
            count = count+1
            # Remove known problem range for Coronado/Point Loma gap
            problem_range = range(221,240)
            if idx not in problem_range:
                geometry.append(geo_tmp)
                MOP.append(MopNames[idx])
                id.append(idx)
    # Create dictionary and initial gdf
    initial_dict = {
        'ID':id, 
        'MOP':MOP, 
        'geometry':geometry
    }
    initial_gdf = gpd.GeoDataFrame(initial_dict,geometry=geometry,crs='epsg:4326')
    # Fix invalid polys and remove overlaps
    valid_gdf = validate_geometry(initial_gdf)
    fix_geo = remove_overlap(valid_gdf)
    fix_dict = {
        'ID':id, 
        'MOP':MOP, 
        'geometry':fix_geo
    }
    # Output both wsg and utm versions of tiles
    bigtile_df_wsg = gpd.GeoDataFrame(fix_dict,geometry=fix_geo,crs='epsg:4326')
    bigtile_df_utm = bigtile_df_wsg.to_crs('epsg:26911')

    return bigtile_df_wsg,bigtile_df_utm

def split_tiles(cliffline_path,bigtile_df_utm,MopsRegion):
    """ 
    Splits "Big" (Cliff + Beach) MOP area tiles into separate cliff and beach tiles 
    using the shapefile of line delineating the cliff base

    Keyword Arguments:
    cliffline_path -- The path to the shapefile of the cliff line
    bigtile_df_utm -- Geodataframe of big tiles in UTM Zone 11N coordinates, created by build_big_tiles
    MopsRegion -- The dataframe returned by the mop_data function.

    Returns:
    MOPtile_df -- Geodataframe of split beach/cliff tiles, including 
        MOP number and Class code (1:beach,0:cliff) for each tile

    The line read from the shapefile is first cleaned, merged, and filled until 
        it is a single shapely LineString object, which is necessary for shapely's
        split() method.
    Next, "check points" are created 10m inshore from the MOP transect offshore point.
        These are used to determine if each post-split tile is a beach or cliff tile.
    """
  
    # Find backbeach line shapefile
    cliffline = gpd.read_file(cliffline_path)

    # If gaps exist, fill them
    if len(cliffline) > 1:
        full_line = fix_line_shp(cliffline)
    else:
        full_line = cliffline.geometry[0]
        
    # Create 'check point' in middle of each beach tile to later assign it a beach class code
    geodesic = pyproj.Geod(ellps='WGS84')
    checkpts = []
    for FID in bigtile_df_utm['ID']:
        fwd_az,_,_ = geodesic.inv(MopsRegion['OffLon'][FID], MopsRegion['OffLat'][FID], 
                                  MopsRegion['BackLon'][FID], MopsRegion['BackLat'][FID])
        dist = 10
        X, Y,_ = geodesic.fwd(MopsRegion['OffLon'][FID], MopsRegion['OffLat'][FID],fwd_az,dist)
        checkpts.append(Point(X,Y))
    checkpts = gpd.GeoSeries(checkpts,crs='epsg:4326')
    checkpts = checkpts.to_crs('epsg:26911')

    # Split big tiles into beach/cliff sections
    splitgeo = []
    for poly in bigtile_df_utm.geometry:
        tempgeo = split(poly,full_line)
        splitgeo.append(tempgeo)
    bigtile_split = bigtile_df_utm.set_geometry(splitgeo)
    # Assign class codes to each tile and build attribute lists
    MOP = []
    FID = []
    MOP_idx = []
    Class = []
    final_poly = []
    beach_tot = 0
    cliff_tot = 0
    for idx in range(len(bigtile_split)):
        # Assign boolean determining if each subtile contains a check point
        check = gpd.GeoSeries(bigtile_split['geometry'][idx].geoms).apply(lambda x: checkpts.within(x).any())
        # Create new entry in each attribute list for each subtile
        if len(bigtile_split['geometry'][idx].geoms) > 1:
            for g_idx in range(len(bigtile_split['geometry'][idx].geoms)):
                MOP.append(bigtile_split['MOP'][idx])
                MOP_idx.append(bigtile_split['ID'][idx])
                final_poly.append(bigtile_split['geometry'][idx].geoms[g_idx])  
                if check[g_idx]:
                    beach_tot = beach_tot+1
                    Class.append(1)
                    FID.append(beach_tot)
                else:
                    cliff_tot = cliff_tot+1
                    Class.append(0)
                    FID.append(len(bigtile_split)+cliff_tot-1)
    # Create Final GeoDataframe from attribute lists
    final_dict = {
        'FID':FID, 
        'MOP_num':MOP,
        'MOP_idx':MOP_idx, 
        'Class':Class, 
        'geometry':final_poly
    }
    MOPtile_df = gpd.GeoDataFrame(final_dict,crs='epsg:26911')

    return MOPtile_df
            
def build_mop_lines(MopsRegion,MopsRange):
    """ 
    Builds MOP Transect Lines in both WSG 1984 and UTM Zone 11N coordinates.

    Keyword Arguments:
    MopsRegion -- The dataframe returned by the mop_data function.

    Returns:
    MOPline_df_wsg -- Geodataframe of big tiles in WSG 1984 coordinates.
    MOPline_df_utm -- Geodataframe of big tiles in UTM Zone 11N coordinates.
    """
    id = []
    MOP = []
    geometry = []
    line_range = MopsRange
    for idx in line_range:
        if idx < line_range[-1]:
            ### builds MOP Lines
            vertsX = [MopsRegion['OffLon'][idx],MopsRegion['BackLon'][idx]]
            vertsY = [MopsRegion['OffLat'][idx],MopsRegion['BackLat'][idx]]
            geo_tmp = LineString(zip(vertsX,vertsY))

            geometry.append(geo_tmp)
            MOP.append(MopsRegion['Name'][idx][-4:])
            id.append(idx)
    intermed_dict = {
        'ID':id, 
        'MOP':MOP, 
        'geometry':geometry
    }
    MOPline_df_wsg = gpd.GeoDataFrame(intermed_dict,geometry=geometry,crs='epsg:4326')
    MOPline_df_utm = MOPline_df_wsg.to_crs('epsg:26911')

    return MOPline_df_wsg,MOPline_df_utm

def writetxt(MOPtile_df,cliffline_path):
    """ 
    Writes FIDs of beach & cliff tiles to respective txt files.

    """

    folder,cliffline_name = os.path.split(cliffline_path)

    # Write beach/cliff txt files
    beach_fid = list(MOPtile_df.loc[MOPtile_df['Class'] == 1]['FID'])
    cliff_fid = list(MOPtile_df.loc[MOPtile_df['Class'] == 0]['FID'])
    path = os.path.join(folder,cliffline_name[:-4])
    # Write beach FIDs to txt file
    with open(path+'_beach.txt',"w+") as fb:
        for item in beach_fid:
            fb.write(f'{item}\n')
    fb.close()
    # Write beach FIDs to txt file
    with open(path+'_cliff.txt',"w+") as fc:
        for item in cliff_fid:
            fc.write(f'{item}\n')
    fc.close()
    return

def build_polys(mop_csv_path,basepoints_path):
    """ 
    Executes previous functions, taking in cliff base line & building MOP polys.

    Full description of process laid out in specific function docstrings.

    Basic Arguments:
    mop_csv_path -- path to csv of MOP data
    basepoints_path -- path to csv of cliff base points, as generated by cliff_delineatool.py

    """

    cliffline_path = basepoints_to_line(basepoints_path)
    folder = os.path.split(cliffline_path)[0]

    # Build Polys based off cliff line
    cliffline_name = os.path.basename(cliffline_path)
    # Get MOP data for area that the cliffline overlaps
    MopsRegion,MopsRange  = mop_data(mop_csv_path, cliffline_path)
    # backdist is how far inland (in meters) from the midpoints the cliff tiles should extend
    backdist = 250
    # Build big tiles containing both beach + cliff areas
    _,bigtile_df_utm = build_big_tiles(MopsRegion, MopsRange, backdist)
    # Split big tiles using cliffline
    MOPtile_df = split_tiles(cliffline_path,bigtile_df_utm,MopsRegion)

    # Save Polys
    poly_name = cliffline_name[:-4]+'_poly.shp'
    save_path = os.path.join(folder,poly_name)
    MOPtile_df.to_file(save_path)
    # Write attribute table to csv file
    header = ['FID','MOP_num','Class']
    csv_path = os.path.join(folder,cliffline_name[:-4]+'_table.txt')
    MOPtile_df.to_csv(csv_path,columns=header)


    return
