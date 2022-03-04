import math
import geopandas as gpd
import pandas as pd
import os
from shapely.geometry import Polygon, LineString

'''
calcWaterCourseLength

Calculates river Network length up/downstream of given feature

Added functionality to stop distance calculation if a blockage is reached (fish migration distance blocked by weirs)

Entirely inspired by Ed B's FlowTrace (C) 2019 [https://github.com/boesiii/flowtrace/blob/master/flowtrace.py].
Re-written for geopandas by Sam R.

GNU General Public License
'''


def calcWaterCourseLength(fid_col, starting_fid, network, direction=0, tolerance=.1):
    '''

    :param fid_col: col containing fid
    :param starting_fid: fid of feature to start the analysis on
    :param network: the river network as a geopandas gdf
    :param direction: upstream (0) or downstream (1)
    :param tolerance: gap to search for next feat, default 0.1
    :return: selected (the selected features), length (length of selected features)
    '''
    # todo convert to espg 27700
    if network.crs != 'EPSG:27700':  # 'EPSG:4269'
        network = network.to_crs('EPSG:27700')


    possible_vertex = list(network[network[fid_col] == starting_fid][fid_col]).copy()
    found_vertex = possible_vertex.copy()

    while possible_vertex:
        # get next coord
        feat = network[network[fid_col] == possible_vertex.pop(0)]
        # get first or last cord dependant on direction
        coords = list(zip(list(feat['geometry'].geometry.iloc[0].coords.xy[0]),
                          list(feat['geometry'].geometry.iloc[0].coords.xy[1])))
        cur_coord = coords[0 - direction]

        # select all proximal fids
        buffer = Polygon([
            (cur_coord[0] - tolerance, cur_coord[1] - tolerance),
            (cur_coord[0] - tolerance, cur_coord[1] + tolerance),
            (cur_coord[0] + tolerance, cur_coord[1] + tolerance),
            (cur_coord[0] + tolerance, cur_coord[1] - tolerance)
        ])

        # potential_fids = index.intersects(buffer)
        potential_fids = list(network[network.intersects(buffer)][fid_col])

        # loop connected features
        for fid in potential_fids:
            if fid not in found_vertex:
                next_feat = network[network[fid_col] == fid].copy()
                coords = list(zip(list(next_feat.geometry.iloc[0].coords.xy[0]),
                                  list(next_feat.geometry.iloc[0].coords.xy[1])))
                next_coord = coords[direction - 1]

                # dist between coords
                dist = math.sqrt((next_coord[0] - cur_coord[0]) ** 2
                                 + (next_coord[1] - cur_coord[1]) ** 2)

                if dist <= tolerance:
                    # if within tolerance adds fid to continue exploration & to the final vertex
                    found_vertex.append(fid)
                    possible_vertex.append(fid)

    selected = network[network[fid_col].isin(found_vertex)].copy()
    return selected, sum(selected.geometry.length)


def calcWaterCourseLengthBlockage(fid_col, starting_fid, network, blockage, direction=0, tolerance=.1):
    '''

    :param fid_col: col containing fid
    :param starting_fid: fid of feature to start the analysis on
    :param network: the river network as a geopandas gdf
    :param direction: upstream (0) or downstream (1)
    :param tolerance: gap to search for next feat, default 0.1
    :param blockage: the blockages that stop river length accumulation i.e. weirs

    :return: selected (the selected features), length (length of selected features)
    '''

    if network.crs != 'EPSG:27700':
        network = network.to_crs('EPSG:27700')
    if blockage.crs != 'EPSG:27700':
        blockage = blockage.to_crs('EPSG:27700')

    possible_vertex = list(network[network[fid_col] == starting_fid][fid_col]).copy()
    found_vertex = possible_vertex.copy()

    while possible_vertex:
        # get next coord
        feat = network[network[fid_col] == possible_vertex.pop(0)]
        # get first or last cord dependant on direction
        coords = list(zip(list(feat['geometry'].geometry.iloc[0].coords.xy[0]),
                          list(feat['geometry'].geometry.iloc[0].coords.xy[1])))
        cur_coord = coords[0 - direction]

        # select all proximal fids
        buffer = Polygon([
            (cur_coord[0] - tolerance, cur_coord[1] - tolerance),
            (cur_coord[0] - tolerance, cur_coord[1] + tolerance),
            (cur_coord[0] + tolerance, cur_coord[1] + tolerance),
            (cur_coord[0] + tolerance, cur_coord[1] - tolerance)
        ])

        # potential_fids = index.intersects(buffer)
        potential_fids = list(network[network.intersects(buffer)][fid_col])

        # loop connected features
        for fid in potential_fids:
            if fid in found_vertex:
                continue
            next_feat = network[network[fid_col] == fid].copy()
            coords = list(zip(list(next_feat.geometry.iloc[0].coords.xy[0]),
                              list(next_feat.geometry.iloc[0].coords.xy[1])))
            next_coord = coords[direction - 1]

            # dist between coords
            dist = math.sqrt((next_coord[0] - cur_coord[0]) ** 2
                             + (next_coord[1] - cur_coord[1]) ** 2)

            if dist <= tolerance:
                the_blockage = blockage[blockage.within(next_feat.buffer(1).unary_union)]
                if the_blockage.empty:
                    # if within tolerance adds fid to continue exploration & to the final vertex
                    found_vertex.append(fid)
                    possible_vertex.append(fid)
                else:
                    # if within tolerance but blockage only add section up til blockage in network to found
                    # & don't add to possible vertex to stop exploration

                    blockage_x, blockage_y = the_blockage.geometry.x.values[0], the_blockage.geometry.y.values[0]
                    dists = [math.sqrt((coord[0] - blockage_x) ** 2
                                       + (coord[1] - blockage_y) ** 2) for coord in coords]
                    line_til_weir = LineString(coords[:dists.index(min(dists))])
                    network.loc[network[fid_col] == fid, 'geometry'] = line_til_weir

                    found_vertex.append(fid)

    selected = network[network[fid_col].isin(found_vertex)].copy()
    return selected, sum(selected.geometry.length)


if __name__ == "__main__":
    dir = os.path.dirname(os.path.realpath(__file__))

    rivers = gpd.read_file(os.path.join(dir, r"test data\RiverWorfeBarriersDRN.shp"))
    weirs = gpd.read_file(os.path.join(dir,r"test data\RiverWorfeBarriers.shp"))

    # test features:
    # eaew1001000000258269 # test id (30 upstream & 26 downstream)
    # eaew1001000000225135 # blockage in next feature 444 ish total len
    # eaew1001000000222782 # blockage about 2500m away, unblocked length 14000m.

    selected, length = calcWaterCourseLength(fid_col='DRN_ID', starting_fid="eaew1001000000258269", network=rivers)
    print(selected, length)
    selected, length = calcWaterCourseLengthBlockage(fid_col='DRN_ID', starting_fid="eaew1001000000222782", network=rivers, blockage=weirs)
    print(selected, length)
