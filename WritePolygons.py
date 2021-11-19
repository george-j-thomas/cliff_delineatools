from shapely.geometry import Polygon,LineString, LinearRing, MultiLineString
from shapely.validation import make_valid
import pandas as pd
import pyproj
import geopandas as gpd
from shapely.geometry.polygon import LinearRing


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

def BuildBigTiles(MopsRegion, backdist):
    """ Builds "Big" (Cliff + Beach) MOP area tiles in both WSG 1984 and UTM Zone 11N coordinates.

    Keyword Arguments:
    MopsRegion -- The dataframe returned by the MopData function.
    backdist -- Inshore distance from backbeach midpoints for cliffback vertexes (see GetCliffBacks function).

    Returns:
    bigtile_df_wsg -- Geodataframe of big tiles in WSG 1984 coordinates.
    bigtile_df_utm -- Geodataframe of big tiles in UTM Zone 11N coordinates.

    First creates backabeach, offshore, and cliffback midpoints 
        (see internal GetMidPoints and GetCliffBacks functions).

    Then builds Polygon "x" from vertices in the following order:
        cliffback midpoint (x), 
        backbeach midpoint (x), 
        offshore midpoint (x),
        offshore MOP transect point (x), 
        offshore midpoint (x+1),
        backbeach midpoint (x+1), 
        cliffback midpoint (x+1)

    Then fixes invalid polygons and removes overlaps using internal
        MakeValid and RemoveOverlaps functions, respectively (both described below).
    """

    def GetMidPoints(MopsRegion):
        """ Returns midpoints between offshore and backbeach MOP points.

        Keyword Arguments:
        MopsRegion -- The dataframe returned by the MopData function.

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
        mops_range = range(len(MopsRegion)) 
        for idx in mops_range:
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
            elif idx > 0 and idx < mops_range[-1]:
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
            elif idx == mops_range[-1]:
                # For last area, similar process to first area
                fwd_az,_,_ = geodesic.inv(MopsRegion['BackLon'][idx-1], MopsRegion['BackLat'][idx-1], 
                                        MopsRegion['BackLon'][idx], MopsRegion['BackLat'][idx])
                mid_dist = 50
                endpointX, endpointY,_ = geodesic.fwd(MopsRegion['BackLon'][idx], MopsRegion['BackLat'][idx],fwd_az,mid_dist)
                backmidX.append(endpointX)
                backmidY.append(endpointY)
                endpointX, endpointY,_ = geodesic.fwd(MopsRegion['OffLon'][idx], MopsRegion['OffLat'][idx],fwd_az,mid_dist)
                offmidX.append(endpointX)
                offmidY.append(endpointY)

        return backmidX,backmidY,offmidX, offmidY

    def GetCliffBacks(backmidX,backmidY,offmidX, offmidY,backdist):
        """ Returns coordinates for onshore vertexes (back of 'Big' tiles),
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

    def MakeValid(initial_gdf):
        """ Fixes "bowtie" polygons and returns a geodataframe with entirely valid polygons.

        Finds invalid polygons, usually ones with "bowtie" shapes, and runs shapely's "make_valid"
            funtion on them. 
        """
        problem_list = list(initial_gdf[~initial_gdf.is_valid]['MOP'])
        valid_geom = []
        for _, row in initial_gdf.iterrows():
            if row['MOP'] in problem_list:
                t1 = make_valid(row['geometry'])
                t2 = max(t1.geoms, key=lambda a:a.area)
                valid_geom.append(t2)
            else:
                valid_geom.append(row['geometry'])
        initial_gdf['valid_geometry'] = valid_geom
        valid_gdf = initial_gdf.drop(columns=['geometry'])
        valid_gdf = valid_gdf.set_geometry('valid_geometry')

        return valid_gdf
    
    def RemoveOverlap(valid_gdf):
        """ Takes in a geodataframe containing entirely valid polygons and removes overlaps.

        Overlaps are removed such that among 2+ overlapping polygons, 
            the last polygon in the list will retain the overlapping area.
        """
        fix_geo = []
        for i, row in valid_gdf.iterrows():
            others = valid_gdf[valid_gdf['ID'] != row['ID']]
            diff = row['valid_geometry'].difference(others.unary_union)
            valid_gdf.at[i,'valid_geometry'] = diff
            fix_geo.append(diff)

        return fix_geo

    # Get Offshore and Backbeach midpoints
    backmidX,backmidY,offmidX, offmidY = GetMidPoints(MopsRegion)

    # Get Cliffback midpoints
    cliffbackX,cliffbackY = GetCliffBacks(backmidX,backmidY,offmidX, offmidY,backdist)

    id = []
    MOP = []
    geometry = []
    tile_range = range(len(MopsRegion))
    for idx in tile_range:
        if idx < tile_range[-1]:
            # Builds BIG tiles (Beach + Cliff) to be split later by backbeach line
            vertsX = [cliffbackX[idx],backmidX[idx],offmidX[idx],MopsRegion['OffLon'][idx],
                      offmidX[idx+1],backmidX[idx+1],cliffbackX[idx+1]]#,cliffbackX[idx]]
            vertsY = [cliffbackY[idx],backmidY[idx],offmidY[idx],MopsRegion['OffLat'][idx],
                      offmidY[idx+1],backmidY[idx+1],cliffbackY[idx+1]]#,cliffbackY[idx]]
            r = LinearRing(list(zip(vertsX,vertsY)))
            geo_tmp = Polygon(r)
            
            # Remove known problem range for Coronado/Point Loma gap
            problem_range = range(221,240)
            if int(MopsRegion['Name'][idx][-4:]) not in problem_range:
                geometry.append(geo_tmp)
                MOP.append(MopsRegion['Name'][idx][-4:])
                id.append(idx)
    # Create dictionary and initial gdf
    initial_dict = {
        'ID':id, 
        'MOP':MOP, 
        'geometry':geometry
    }
    initial_gdf = gpd.GeoDataFrame(initial_dict,geometry=geometry,crs='epsg:4326')
    # Fix invalid polys and remove overlaps
    valid_gdf = MakeValid(initial_gdf)
    fix_geo = RemoveOverlap(valid_gdf)
    fix_dict = {
        'ID':id, 
        'MOP':MOP, 
        'geometry':fix_geo
    }
    # Output both wsg and utm versions of tiles
    bigtile_df_wsg = gpd.GeoDataFrame(fix_dict,geometry=fix_geo,crs='epsg:4326')
    bigtile_df_utm = bigtile_df_wsg.to_crs({'init':'epsg:26911'})

    return bigtile_df_wsg,bigtile_df_utm

def SplitTiles(cliffline_path,bigtile_df_utm,MopsRegion):
    """ Splits "Big" (Cliff + Beach) MOP area tiles into separate cliff and beach tiles 
        using the shapefile of line delineating the cliff base

    Keyword Arguments:
    cliffline_path -- The path to the shapefile of the cliff line
    bigtile_df_utm -- Geodataframe of big tiles in UTM Zone 11N coordinates, created by BuildBigTiles
    MopsRegion -- The dataframe returned by the MopData function.

    Returns:
    MOPtile_df -- Geodataframe of split beach/cliff tiles, including 
        MOP number and Class code (1:beach,0:cliff) for each tile

    The line read from the shapefile is first cleaned, merged, and filled until 
        it is a single shapely LineString object, which is necessary for shapely's
        split() method.
    Next, "check points" are created 10m inshore from the MOP transect offshore point.
        These are used to determine if each post-split tile is a beach or cliff tile.
    """
    from shapely.ops import split, linemerge
    from shapely.geometry import Point

    # Find backbeach line shapefile
    cliffline = gpd.read_file(cliffline_path)
    # Merge whole line together, filling gaps along the way
    cliffline_multi = MultiLineString(list(cliffline['geometry']))
    cliffline_merge = linemerge(cliffline_multi)
    # Remove 'bad' cliff linestrings that are only 2 points
    goodline_list = []
    for line in cliffline_merge.geoms:
        if len(line.coords) > 2:
            goodline_list.append(line)
    # Fill Gaps in between 'good' cliff linestrings with small straight lines 
    # (Need continuous linestring for splitting)
    connectors = []
    for line_i in range(len(goodline_list)-1):
        first = Point(goodline_list[line_i].coords[0])
        last = Point(goodline_list[line_i+1].coords[-1])
        connectors.append(LineString([first,last]))
    full_list = list(goodline_list) + list(connectors)
    full_line = linemerge(full_list)
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
    checkpts = checkpts.to_crs({'init':'epsg:26911'})

    # Split big tiles into beach/cliff sections
    splitgeo = []
    for poly in bigtile_df_utm.geometry:
        tempgeo = split(poly,full_line)
        splitgeo.append(tempgeo)
    bigtile_split = bigtile_df_utm.set_geometry(splitgeo)
    # Assign class codes to each tile and build attribute lists
    MOP = []
    FID = []
    Class = []
    final_poly = []
    beach_tot = 0
    cliff_tot = 0
    for idx in range(len(bigtile_split)):
        # Assign boolean determining if each subtile contains a check point
        check = gpd.GeoSeries(bigtile_split['geometry'][idx].geoms).apply(lambda x: checkpts.within(x).any())
        # Create new entry in each attribute list for each subtile
        for g_idx in range(len(bigtile_split['geometry'][idx].geoms)):
            MOP.append(bigtile_split['MOP'][idx])
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
        'MOP':MOP, 
        'Class':Class, 
        'geometry':final_poly
    }
    MOPtile_df = gpd.GeoDataFrame(final_dict,crs='epsg:26911')

    return MOPtile_df
            
def BuildMOPLines(MopsRegion):
    """ Builds MOP Transect Lines in both WSG 1984 and UTM Zone 11N coordinates.

    Keyword Arguments:
    MopsRegion -- The dataframe returned by the MopData function.

    Returns:
    MOPline_df_wsg -- Geodataframe of big tiles in WSG 1984 coordinates.
    MOPline_df_utm -- Geodataframe of big tiles in UTM Zone 11N coordinates.
    """
    id = []
    MOP = []
    geometry = []
    line_range = range(len(MopsRegion))
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
    MOPline_df_utm = MOPline_df_wsg.to_crs({'init':'epsg:26911'})

    return MOPline_df_wsg,MOPline_df_utm

csv_path = r"C:\Users\g4thomas\Documents\CliffDelineation\Mop.csv"
SDMops = range(931)
MopsRange = SDMops
MopsRegion = MopData(csv_path, MopsRange)
backdist = 250
bigtile_df_wsg,bigtile_df_utm = BuildBigTiles(MopsRegion,backdist)
MOPline_df_wsg,MOPline_df_utm = BuildMOPLines(MopsRegion)
cliffline_path = r"C:\Users\g4thomas\Documents\CliffDelineation\BackBeachLine.shp"
MOPtile_df = SplitTiles(cliffline_path,bigtile_df_utm,MopsRegion)

MOPtile_df.to_file("MOPtiles_py.shp")
