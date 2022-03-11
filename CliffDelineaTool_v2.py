# CliffDelineaToolPy v1.0.0
# An algorithm to map coastal cliff base and top from topography
# https://github.com/zswirad/CliffDelineaTool
# Zuzanna M Swirad (zswirad@ucsd.edu), Scripps Institution of Oceanography, UC San Diego
# Originally coded in MATLAB 2019a (v1.2.0; 2021-11-24)
# Help in debugging: George Thomas (2021-12-2)
# Last updated on 2021-12-3 (Python 3.8)

import os
import pandas as pd
import numpy as np
import math
import statsmodels.api as sm
import glob

# Indicate input data parameters:
os.chdir(r"C:\Users\g4thomas\Documents\CliffDelineation\Files\2_To_Delineate_txt")

fileName = '20211103_00518_00568_NoWaves_Blacks_CliffPoints.txt'#'20200220_00982_01170_NoWaves_MASS_Camp_Pendleton_CliffPoints.txt'
top = False

# Set the calibrated input variables:
nVert = 5#20 # How many adjacent points to consider as a local scale?
baseMaxElev = 10 # What is the top limit for cliff base elevation (m)?
baseSea = 20#15; # What is the max seaward slope for cliff base (deg)?
baseLand = 20#25; # What is the min landward slope for cliff base (deg)? - seems to have larger impact than bS
topSea = 20 # What is the min seaward slope for cliff top (deg)?
topLand = 15 # What is the max landward slope for cliff top (deg)?
propConvex = 0.5 # What is the minimal proportion of the distance from trendline #2 to replace modelled cliff top location?
smoothWindow = 10 # What is the alongshore moving window for cross-shore smoothing (points)? # changing to 5 did nothing

def get_slopes(data):
    """ Get slopes of transect lines for later delineation.

    Takes in original cliff points data and returns table including:
        -SeaSlope
        -LandSlope
        -Trendline1
        -Difference1
    """
    table = pd.DataFrame()

    for n in range(min(data['TransectID']), max(data['TransectID'])+1):
        sub = data[data['TransectID'] == n]
        rowCount = sub.shape[0]
        if rowCount > 0:
            sub = sub.sort_values(['Distance'])
            sub = sub.reset_index(drop=True)
            
            # Fill data gaps:
            sub.loc[sub['Elevation'] < -50, ['Elevation']] = np.nan
            sub = sub.interpolate()
            sub = sub.fillna(method = 'ffill')
            sub = sub.fillna(method = 'bfill')
            
            # Calculate local slopes:
            zeros = np.zeros(rowCount + 1)
            sub['SeaSlope'] = pd.Series(zeros) # seaward slope (average slope between the point and nVert consecutive seaward points)
            sub['LandSlope'] = pd.Series(zeros) # landward slope (average slope between the point and nVert consecutive landward points)
            sub['Trendline1'] = pd.Series(zeros) # trendline #1
            sub['Difference1'] = pd.Series(zeros) # elevations - trendline #1
            sub = sub.fillna(0)
            
            for z in range(nVert, rowCount - nVert):
                count = 0
                for s in range(1, nVert + 1):
                    if sub.Elevation[z] != sub.Elevation[z-s]: 
                        angle = math.degrees(math.atan((sub.Elevation[z] - sub.Elevation[z-s])/(sub.Distance[z] - sub.Distance[z-s])))
                        if angle < 0:
                            angle = 0
                            
                        # sub.SeaSlope[z] = sub.SeaSlope[z] + angle
                        sub.loc[z,'SeaSlope'] = sub.loc[z,'SeaSlope'] + angle

                        count += 1
    
                # sub.SeaSlope[z] = sub.SeaSlope[z]/count
                sub.loc[z,'SeaSlope'] = sub.loc[z,'SeaSlope']/count
                
                count = 0
                for s in range(1, nVert + 1):
                    if sub.Elevation[z] != sub.Elevation[z-s]:
                        angle = math.degrees(math.atan((sub.Elevation[z+s] - sub.Elevation[z])/(sub.Distance[z+s] - sub.Distance[z])))
                        if angle < 0:
                            angle = 0
                            
                        # sub.LandSlope[z] = sub.LandSlope[z] + angle
                        sub.loc[z,'LandSlope'] = sub.loc[z,'LandSlope'] + angle

                        count += 1
    
                # sub.LandSlope[z] = sub.LandSlope[z]/count
                sub.loc[z,'LandSlope'] = sub.loc[z,'LandSlope']/count

            # Limit the transect landwards to the highest point + nVert:
            indMax = np.argmax(sub['Elevation'], axis=0)
            # If the highest point is more than 250 m along the transect, limit it to 250 
            if indMax > 250:
                indMax = 250 - nVert
            # Limit the transect landwards to indMax
            if rowCount > indMax + nVert:
                allDrop = rowCount - (indMax + nVert + 1)
                sub.drop(sub.tail(allDrop).index,inplace = True)    
            rowCount = sub.shape[0]
            sub = sub.reset_index(drop=True)
            
            # Draw trendline #1 (straight line between the seaward and landward transect ends):
            sub.loc[0,'Trendline1'] = sub.loc[0,'Elevation']
            sub.iloc[-1,sub.columns.get_loc('Trendline1')] = sub.iloc[-1, sub.columns.get_loc('Elevation')]
            for z in range(1, rowCount):
                sub.loc[z,'Trendline1'] = ((sub.loc[z,'Distance'] - sub.loc[0,'Distance'])
                                    * (sub.iloc[-1, sub.columns.get_loc('Elevation')] - sub.loc[0,'Elevation']) 
                                    / (sub.iloc[-1, sub.columns.get_loc('Distance')] - sub.loc[0,'Distance']) 
                                    + sub.loc[0,'Elevation'])
                
            # Calculate vertical distance between actual elevations and trendline #1:
            sub['Difference1'] = sub['Elevation'] - sub['Trendline1']
                
            table = table.append(sub)
            
    table = table.reset_index(drop=True)
    return table

def model_base(table):
    """ Finds potential cliff base locations using input variables

    """
    # Find potential cliff base locations:
    potential_base = table[(table['Elevation'] < baseMaxElev) & (table['SeaSlope'] < baseSea) & (table['LandSlope'] > baseLand) & (table['Difference1'] < 0)]   
    
    # From the points that satisfy the criteria, for each transect select one with the largest vertical difference between the elevation and trendline #1:
    modelled_base = pd.DataFrame() 
    if potential_base.shape[0] > 0:
        cliffed_profiles = potential_base['TransectID'].unique()
        for n in range(potential_base['TransectID'].min(), potential_base['TransectID'].max()+1):
            for m in range(cliffed_profiles.shape[0]):
                if n == cliffed_profiles[m]:
                    sub = potential_base[potential_base['TransectID'] == n]
                    sub = sub.sort_values(by=['Difference1'])
                    modelled_base = modelled_base.append(sub.iloc[0])
    return modelled_base, cliffed_profiles

def model_top(modelled_base, cliffed_profiles,table):
    """ Finds potential cliff top locations

    """
    modelled_top = pd.DataFrame()
    for n in range(int(modelled_base['TransectID'].min()), int(modelled_base['TransectID'].max()+1)):
        for m in range(cliffed_profiles.shape[0]):
            if n == cliffed_profiles[m]:
                sub = table[table['TransectID'] == n]
                sub = sub.reset_index(drop=True)
                
                # Remove points seawards from the cliff base:
                sub_base = modelled_base[modelled_base['TransectID'] == n]   
                sub_base = sub_base.reset_index(drop=True)
                sub_base_dist = sub_base.Distance[0]
                sub.drop(sub[sub['Distance'] < sub_base_dist].index, inplace = True)
                sub = sub.reset_index(drop=True)
                
                # Draw trendline #2 between cliff base and landward transect end:
                rowCount = sub.shape[0]
                zeros = np.zeros(rowCount + 1)
                sub['Trendline2'] = pd.Series(zeros) # trendline #2
                sub['Difference2'] = pd.Series(zeros) # elevation - trendline #2
                sub = sub.fillna(0)

                sub.loc[0,'Trendline2'] = sub.loc[0,'Elevation']
                sub.iloc[-1,sub.columns.get_loc('Trendline2')] = sub.iloc[-1, sub.columns.get_loc('Elevation')]
                for z in range(1, rowCount):
                    sub.loc[z,'Trendline2'] = ((sub.Distance[z]-sub.Distance[0]) 
                                        * (sub.iloc[-1, sub.columns.get_loc('Elevation')]-sub.Elevation[0]) 
                                        / (sub.iloc[-1, sub.columns.get_loc('Distance')]-sub.Distance[0]) 
                                        + sub.Elevation[0])
                sub['Difference2'] = sub['Elevation'] - sub['Trendline2']
                
                # Find potential cliff top locations:
                potential_top = sub[(sub['SeaSlope'] > topSea) & (sub['LandSlope'] < topLand) & (sub['Difference2'] > 0)]
                
                if potential_top.shape[0] > 0:
                    potential_top = potential_top.sort_values(by=['Difference2'])
                    
                    # From the points that satisfy the criteria, for each transect select one with the largest vertical difference between the elevation and trendline #2:
                    modelled_top0 = potential_top.iloc[-1]   
                    
                    # Check whether the selected point is part of within-cliff flattening:
                    if potential_top['Distance'].max() > modelled_top0.Distance + nVert:
                        subNew = sub.copy()
                        sub_top_dist = potential_top.iloc[-1, sub.columns.get_loc('Distance')]
                        subNew.drop(subNew[subNew['Distance'] < sub_top_dist].index, inplace = True) # remove points seawards from the modelled cliff top
                        rowCountNew = subNew.shape[0]
                        zerosNew = np.zeros(rowCountNew + 1)
                        subNew['Trendline3'] = pd.Series(zerosNew)
                        subNew['Difference3'] = pd.Series(zerosNew)
                        subNew = subNew.fillna(0)
                        subNew = subNew.reset_index(drop=True)
                        
                        subNew.loc[0,'Trendline3'] = subNew.loc[0,'Elevation']
                        subNew.iloc[-1, subNew.columns.get_loc('Trendline3')] = subNew.iloc[-1, subNew.columns.get_loc('Elevation')]                            
                        for z in range(1, rowCountNew):
                            subNew.loc[z,'Trendline3'] = ((subNew.Distance[z]-subNew.Distance[0]) 
                                                        * (subNew.iloc[-1, subNew.columns.get_loc('Elevation')]-subNew.Elevation[0]) 
                                                        / (subNew.iloc[-1, subNew.columns.get_loc('Distance')]-subNew.Distance[0]) 
                                                        + subNew.Elevation[0])
                        subNew['Difference3'] = subNew['Elevation'] - subNew['Trendline3']
                        
                        potential_top2 = potential_top.copy()
                        potential_top2.drop(potential_top2[potential_top2['Distance'] < modelled_top0.Distance].index, inplace = True)
                        rowCountNew = potential_top2.shape[0]
                        zerosNew = np.zeros(rowCountNew + 1)
                        potential_top2['Difference3'] = pd.Series(zerosNew)
                        potential_top2 = potential_top2.fillna(0)
                        potential_top2 = potential_top2.reset_index(drop=True)
                        
                        for p in range(potential_top2.shape[0]):
                            subNewTemp = subNew[subNew['Distance'] == potential_top2.Distance[p]]
                            potential_top2.iloc[p,potential_top2.columns.get_loc('Difference3')] = subNewTemp.Difference3

                        potential_top2 = potential_top2[(potential_top2['Difference3'] > 0) & (potential_top2['Difference2'] >= modelled_top0.Difference2*propConvex) & (potential_top2['Distance'] >= modelled_top0.Distance + nVert)]
                        if potential_top2.shape[0] > 0:
                            potential_top2 = potential_top2.sort_values(by=['Difference2'])
                            potential_top2.drop(['Difference3'], axis=1)
                            modelled_top0 = potential_top2.iloc[-1]
                        
                    modelled_top = modelled_top.append(modelled_top0)
    return modelled_top

def fix_top(fix, modelled_top, modelled_base, table):
    """ Re-models cliff top locations, removing outliers
    
    """
    fix = fix.reset_index(drop=True)
    modelled_top.drop(['StandResidual', 'Outlier'], axis=1) 
    for c in range(fix.shape[0]):
        sub = table[table['TransectID'] == fix.TransectID[c]]
        sub = sub.reset_index(drop=True)
        outlier = modelled_top[modelled_top['TransectID'] == fix.TransectID[c]]
        
        # Remove points seawards from the cliff base:
        sub_base = modelled_base[modelled_base['TransectID'] == c] # THIS ORIGINALLY HAD "== n" which I don't think makes sense   
        sub_base = sub_base.reset_index(drop=True)
        sub_base_dist = sub_base.Distance[0]
        sub.drop(sub[sub['Distance'] < sub_base_dist].index, inplace = True)
        sub = sub.reset_index(drop=True)
        
        # Draw trendline #2 between cliff base and landward transect end:
        rowCount = sub.shape[0]
        zeros = np.zeros(rowCount + 1)
        sub['Trendline2'] = pd.Series(zeros) # trendline #2
        sub['Difference2'] = pd.Series(zeros) # elevation - trendline #2
        sub = sub.fillna(0)

        sub.loc[0,'Trendline2'] = sub.loc[0,'Elevation']
        sub.iloc[-1, sub.columns.get_loc('Trendline2')] = sub.iloc[-1, sub.columns.get_loc('Elevation')]
        for z in range(1, rowCount):
            sub.loc[z,'Trendline2'] = ((sub.Distance[z]-sub.Distance[0]) 
                                        * (sub.iloc[-1, sub.columns.get_loc('Elevation')]-sub.Elevation[0]) 
                                        / (sub.iloc[-1, sub.columns.get_loc('Distance')]-sub.Distance[0]) 
                                        + sub.Elevation[0])
        sub['Difference2'] = sub['Elevation'] - sub['Trendline2']
        
        # Find potential cliff top locations:
        potential_top = sub[(sub['SeaSlope'] > topSea) & (sub['LandSlope'] < topLand) & (sub['Difference2'] > 0)]  
        rowCount = potential_top.shape[0]          
        zeros = np.zeros(rowCount + 1)
        potential_top['SmoothedDistance'] = pd.Series(zeros) # smoothed distance
        potential_top['DistanceFromSmoothed'] = pd.Series(zeros) # distance from smoothed distance
        
        potential_top = potential_top.fillna(0)
        potential_top['SmoothedDistance'] = fix.SmoothedDistance[c]
        potential_top['DistanceFromSmoothed'] = abs(potential_top['Distance'] - potential_top['SmoothedDistance'])
        potential_top = potential_top.sort_values(by=['DistanceFromSmoothed'])
        potential_top = potential_top.iloc[0]
        
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'PointID'] = potential_top['PointID']
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'Elevation'] = potential_top['Elevation']
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'Distance'] = potential_top['Distance']
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'SeaSlope'] = potential_top['SeaSlope']
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'LandSlope'] = potential_top['LandSlope']
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'Trendline1'] = potential_top['Trendline1']
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'Difference1'] = potential_top['Difference1']
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'Trendline2'] = potential_top['Trendline2']
        modelled_top.loc[(modelled_top['TransectID'] == potential_top['TransectID']),'Difference2'] = potential_top['Difference2']

    rowCount = modelled_top.shape[0]          
    zeros = np.zeros(rowCount + 1)
    modelled_top['StandResidual'] = pd.Series(zeros) # standardized residuals
    modelled_top['Outlier'] = pd.Series(zeros) # outliers
    modelled_top = modelled_top.fillna(0)
    
    model = sm.OLS(modelled_top['Distance'],modelled_top['SmoothedDistance'])
    results = model.fit()
    influence = results.get_influence()
    modelled_top['StandResidual'] = influence.resid_studentized_internal               
    modelled_top.loc[abs(modelled_top['StandResidual']) > 2, ['Outlier']] = 1 

    modelled_top.drop(modelled_top[modelled_top['Outlier'] == 1].index, inplace = True) # ignore new cliff top positions if standardized residuals did not improve
    return modelled_top
    
def delineate_cliff(fileName,top,nVert=20,baseMaxElev=10,baseSea=20,baseLand=20,topSea=20,topLand=15,propConvex=0.5,smoothWindow=10):
    with open(fileName, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(fileName, 'w') as fout:
        fout.writelines(data[1:])

    data = pd.read_csv(fileName, header = None)
    data.columns =['PointID', 'TransectID', 'Elevation', 'Distance','UTM_E','UTM_N'] # Give column names; you can add XY coordinates etc if they were exported from GIS software
    data = data[data.Elevation != -9999]

    table = pd.DataFrame()
   
    table = get_slopes(data)
            
    modelled_base, cliffed_profiles = model_base(table)
     
    # Find cliff top locations for transects with cliff base:
    if modelled_base.shape[0] > 0:    

        # Save the data:
        saveName1 = fileName[:9]+fileName[29:-16] + '_nVert{}_bME{}_bS{}_bL{}_tS{}_tL{}_pC{}_sW{}_base.txt'.format(nVert,baseMaxElev,baseSea,baseLand,topSea,topLand,propConvex,smoothWindow)

        folder = r"C:\Users\g4thomas\Documents\CliffDelineation\Files\3_Delineated_csv"
        outpath = os.path.join(folder,saveName1)
        modelled_base_save = modelled_base[['PointID', 'TransectID','UTM_E','UTM_N']] # Select which columns to save; you may want to add XY coordinates if they were present
        modelled_base_save.to_csv(outpath, header=False, index = False) # change to header=True if exporting with header     

        if top:
            modelled_top = model_top(modelled_base, cliffed_profiles,table)
                            
            # Remove alongshore outliers:
            # 1. Find outliers:   
            modelled_base = modelled_base.sort_values(by=['TransectID'])
            modelled_top = modelled_top.sort_values(by=['TransectID'])
            
            rowCount = modelled_top.shape[0]          
            zeros = np.zeros(rowCount + 1)
            modelled_top['SmoothedDistance'] = pd.Series(zeros) # smoothed distance
            modelled_top['StandResidual'] = pd.Series(zeros) # standardized residuals
            modelled_top['Outlier'] = pd.Series(zeros) # outliers (https://online.stat.psu.edu/stat462/node/172/; accessed on 2021/06/04)
            modelled_top = modelled_top.fillna(0)
            
            modelled_top['SmoothedDistance'] = modelled_top['Distance'].rolling(window = smoothWindow).median()
            modelled_top['SmoothedDistance'] = modelled_top['SmoothedDistance'].fillna(method = 'ffill')
            modelled_top['SmoothedDistance'] = modelled_top['SmoothedDistance'].fillna(method = 'bfill')
            
            model = sm.OLS(modelled_top['Distance'],modelled_top['SmoothedDistance'])
            results = model.fit()
            influence = results.get_influence()
            modelled_top['StandResidual'] = influence.resid_studentized_internal
            
            modelled_top.loc[abs(modelled_top['StandResidual']) > 2, ['Outlier']] = 1
    
            fix = modelled_top[modelled_top['Outlier'] == 1]
            
            # 2. Delete or replace outliers with more suitable potential cliff tops:
            # (Repeat cliff top detection for the transects with outliers.)
            if fix.shape[0] > 0:
                modelled_top = fix_top(fix, modelled_top, modelled_base, table)

            if modelled_top.shape[0] > 0:
                saveName2 = fileName[:-4] + '_top.txt'
                modelled_top_save = modelled_top[['PointID', 'TransectID','UTM_E','UTM_N']] # Select which columns to save; you may want to add XY coordinates if they were present
                outpath = os.path.join(folder,saveName2)
                modelled_top_save.to_csv(outpath, header=False, index = False) # change to header=True if exporting with header    

def average_cliffline(fileName,top=False):
    import geopandas as gpd
    from shapely.geometry import LineString 
    folder = r"C:\Users\g4thomas\Documents\CliffDelineation\Files\3_Delineated_csv"
    fulldata = pd.DataFrame(columns =['PointID', 'TransectID','UTM_E','UTM_N'])
    search = os.path.join(folder,fileName[:9]+fileName[29:-16]+r'*base.txt')
    for infile in glob.glob(search):
        data = pd.read_csv(infile, header = None)
        data.columns =['PointID', 'TransectID','UTM_E','UTM_N'] 
        fulldata = fulldata.append(data,ignore_index=True)
    transects = sorted(fulldata.TransectID.unique())
    utm_e = []
    utm_n = []
    for tID in transects:
        utm_e.append(np.mean(fulldata['UTM_E'].loc[fulldata['TransectID'] == tID]))
        utm_n.append(np.mean(fulldata['UTM_N'].loc[fulldata['TransectID'] == tID]))
    coordinates = list(zip(utm_e,utm_n))
    line = LineString(coordinates)
    gdr = gpd.GeoSeries(line,crs='epsg:26911')
    fol1,_ = os.path.split(folder)
    outpath = os.path.join(fol1,'4_Cliff_Lines',fileName[:9]+fileName[29:-16]+r'_base_avg.shp')
    gdr.to_file(outpath)

    if top:
        search = os.path.join(folder,fileName[:9]+fileName[29:-16]+r'*top.txt')
        for infile in glob.glob(search):
            data = pd.read_csv(infile, header = None)
            data.columns =['PointID', 'TransectID','UTM_E','UTM_N'] 
            fulldata = fulldata.append(data,ignore_index=True)
        transects = fulldata.TransectID.unique()
        utm_e = []
        utm_n = []
        for tID in transects:
            utm_e.append(np.mean(fulldata['UTM_E'].loc[fulldata['TransectID'] == tID]))
            utm_n.append(np.mean(fulldata['UTM_N'].loc[fulldata['TransectID'] == tID]))
        coordinates = list(zip(utm_e,utm_n))
        line = LineString(coordinates)
        gdr = gpd.GeoSeries(line,crs='epsg:26911')
        fol1,_ = os.path.split(folder)
        outpath = os.path.join(fol1,'4_Cliff_Lines',fileName[:9]+fileName[29:-16]+r'_top_avg.shp')
        gdr.to_file(outpath)

    return

def nVert_test(fileName,top):
    n_values = np.linspace(5,25,5)
    for nVert in n_values:
        delineate_cliff(fileName,top,nVert)