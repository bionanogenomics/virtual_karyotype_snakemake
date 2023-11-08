import sys
import os
sys.path.append(os.getcwd())
global centro
centro = {}
import pandas as pd
from OMKar.main import (
    Vertices, Graph, find_bp_in_segment, merge_list, detect_sv_directions,
    find_nodes, find_start_end, next_prev_label, detect_overlap_map, Plot_graph,
    find_in_smap, dfs, find_connected_components, return_all_edges_in_cc,
    calculate_seg_length, estimating_edge_multiplicities_in_CC, remove_edge,
    dfs_count, isValidNextEdge, check_traverse_segment, printEulerUtil,
    detect_segment_vertices, detect_segment_odd_degree, scoring_paths,
    printEulerTour, detect_del_dup_cn, detect_duplicatioon_inversion_cn,
    detect_receprical_translocation, check_non_centromeric_path,
    convert_path_to_segment, check_exiest_call, extend_segments_cn,
    average_values_greater_than_a, cn_in_mask_N_region, is_overlapping,
    merge_segments_all_seg_smap, reverse_path, find_indices_with_sum_of_2,
    share_same_segments, convert_segment_to_path, swap_segment, fix_dicentric
)
from collections import defaultdict
from os.path import exists
from parsers import *
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import pulp as p
import argparse
import itertools
from copy import deepcopy
import math

def map_cyto_coords(coords, cytoband_filtered, orientation):
    """
    Map cytogenetic coordinates based on the provided cytoband data.
    
    Args:
        coords (str): Coordinate string in the format 'chrom:start-end'.
        cytoband_filtered (pd.DataFrame): Filtered cytoband data.
        orientation (str): Orientation of the coordinates.
    
    Returns:
        str: A formatted string with chromosome, band, start, end, and orientation info.
    """
    chrom = coords.split(':')[0].replace('chr','')
    start,end = map(float,coords.split(':')[1].split('-'))
    start = int(start)
    end = int(end)
    band = cytoband_filtered[(cytoband_filtered['chr'] == chrom) & ((cytoband_filtered['start'] >= start) | (cytoband_filtered['end'] >= start)) & ((cytoband_filtered['end'] <= end)|(cytoband_filtered['start'] <= end))]['band'].values
    if len(band) > 1:
        band = '{}{}'.format(band[0],band[-1])
    if len(band) == 1:
        band = band[0]
    out_handle = "{chrom}{band}({start}-{end})x1({orientation})".format(chrom=chrom,band=band,start=start,end=end,orientation=orientation)
    return out_handle

def read_in_cyto(cytobands = 'resources/cytoBand.txt'):
    """
    Read and process cytoband data from a given file.
    
    Args:
        cytobands (str): Path to the cytoband file. Default is 'resources/cytoBand.txt'.
    
    Returns:
        pd.DataFrame: Filtered and processed cytoband data.
    """
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
    cytoband = pd.read_table(cytobands,header=None)
    cytoband.columns = ['chr','start','end','band','stain']
    cytoband_filtered = cytoband[cytoband['chr'].isin(chroms)]
    cytoband_filtered['chr'] = cytoband_filtered['chr'].str.replace("chr", '').replace('X','23').replace('Y','24')
    cytoband_filtered['band_base'] = cytoband_filtered['band'].str.split('.',expand=True)[0]
    cytoband_filtered['start'] = cytoband_filtered['start'].astype(int)
    cytoband_filtered['end'] = cytoband_filtered['end'].astype(int)
    return cytoband_filtered

def node_to_map(svs, xmap, g):
    """
    Create mapping dictionaries for nodes based on SVs.
    
    Args:
        svs (list): List of structural variations.
        xmap (xmap.parser): Information about the xmap.
        g (Graph): Graph structure.
    
    Returns:
        tuple: Two dictionaries - mapping of nodes to map and nodes to smap.
    """
    node_to_map_dict = {}
    node_to_smap_dict = {}
    for sv in svs:
        q_id = sv.q_id
        smap_id = sv.smap_id
        n_type1, n_type2 = detect_sv_directions(sv, xmap)
        #a and b are two nodes that ara connected by sv edge
        a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type1)
        b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type2)
        a = str(a)
        b = str(b)
        if a not in node_to_map_dict:
            node_to_map_dict[a] = []
            node_to_smap_dict[a] = []
        node_to_map_dict[a].append(q_id)
        node_to_smap_dict[a].append(smap_id)
        if b not in node_to_map_dict:
            node_to_map_dict[b] = []
            node_to_smap_dict[b] = []
        node_to_map_dict[b].append(q_id)
        node_to_smap_dict[b].append(smap_id)
    return node_to_map_dict, node_to_smap_dict

def return_mapids(v,u,node_to_map_dict):
    """
    Return map IDs for the given nodes v and u.
    
    Args:
        v (int): Node v.
        u (int): Node u.
        node_to_map_dict (dict): Mapping dictionary of nodes to map IDs.
    
    Returns:
        str: String representation of intersecting nodes and differences between nodes v and u.
    """
    v = str(v)
    u = str(u)
    if v in node_to_map_dict:
        start_node_maps = node_to_map_dict[v]
    else:
        start_node_maps = ''
    if u in node_to_map_dict:
        end_node_maps = node_to_map_dict[u]
    else:
        end_node_maps = ''
    intersect_nodes = set(start_node_maps).union(set(end_node_maps))
    start_dif = set(start_node_maps) - set(end_node_maps)
    end_dif = set(end_node_maps) - set(start_node_maps)
    map_str = '{}|{}|{}'.format(','.join(map(str,intersect_nodes)),','.join(map(str,start_dif)),','.join(map(str,end_dif)))
    return map_str

def convert_path(structure, path_map, cytoband_filtered):
    """
    Convert the provided structure into merged coordinates and ISCN coordinates.
    
    Args:
        structure (str): String representation of the structure.
        path_map (dict): Mapping of paths to coordinates.
        cytoband_filtered (pd.DataFrame): Filtered cytoband data.
    
    Returns:
        tuple: Merged coordinates and ISCN coordinates as strings.
    """
    structure_split = structure.split()
    coord_list = []
    iscn_list = []
    for s in structure_split:
        orientation = s[-1]
        p = s[:-1]
        coords = path_map[p]
        update_coords = '{}({})'.format(coords,orientation)
        coord_list.append(update_coords)
        iscn_mapped = map_cyto_coords(coords, cytoband_filtered, orientation)
        iscn_list.append(iscn_mapped)
    merged_coords = ','.join(coord_list)
    iscn_coords = ','.join(iscn_list)
    return merged_coords, iscn_coords

def smap_to_segment(pathnumber, smapids, smap):
    """ Associates smap ids with pathnumber, pulls smap entries

    Args:
        pathnumber (int): path 
        smapids (_type_): _description_
        smap (list): list of parsed smap entries from parsers.SmapEntry
    """
    all_ids = smapids.split('|')[0].split(',')
    frame_list = []
    if (len(all_ids) >=1) & (all_ids[0] != ''):
        for smapid in all_ids:
            smap_entry = find_in_smap(int(smapid),smap)
            smap_frame = pd.DataFrame(columns=['smap_id','q_id','ref_c_id1','ref_c_id2','ref_start','ref_end','confidence','sv_type','size','VAF'],data=[[smapid, smap_entry.q_id, smap_entry.ref_c_id1, smap_entry.ref_c_id2, smap_entry.ref_start, smap_entry.ref_end, smap_entry.confidence, smap_entry.sv_type, smap_entry.size, smap_entry.VAF]],index=['Segment {}'.format(pathnumber)])
            frame_list.append(smap_frame)
    if len(frame_list)>0:
        joined_frame_list = pd.concat(frame_list)
    else:
        joined_frame_list = None
    return joined_frame_list

def find_sv_node_edges(svs, xmap, g):
    """
    Identify nodes in a graph that are connected by a structural variation (SV) edge.

    Parameters:
    - svs (list): A list of structural variations (SVs) where each SV has attributes:
        * q_id: Query identifier of the SV.
        * smap_id: Smap identifier of the SV.
        * ref_c_id1: Reference chromosome ID for the start position of the SV.
        * ref_start: Start position of the SV on the reference chromosome.
        * ref_c_id2: Reference chromosome ID for the end position of the SV.
        * ref_end: End position of the SV on the reference chromosome.
    - xmap (object): An object representing the cross-map for the SVs, used by the `detect_sv_directions` function.
    - g (object): A graph object with an attribute `vertices` representing all the vertices (or nodes) in the graph.

    Returns:
    - set: A set of tuple pairs, where each tuple (a, b) represents two nodes (a and b) in the graph 
           that are connected by an SV edge. The set includes both (a, b) and (b, a) for bidirectionality.

    Note:
    The function makes use of two additional functions:
    1. `detect_sv_directions(sv, xmap)` - Determines the directions/types of the nodes based on the SV and xmap.
    2. `find_nodes(ref_c_id, position, vertices, n_type)` - Finds nodes in the graph based on chromosome ID, 
       position, list of vertices, and node type.

    Example:
    Consider svs as a list where each element has attributes like q_id, smap_id, ref_c_id1, etc.
    Given svs, xmap, and g, the function will return a set of node pairs connected by SV edges.
    """
    tuples_list = []
    for sv in svs:
        q_id = sv.q_id
        smap_id = sv.smap_id
        n_type1, n_type2 = detect_sv_directions(sv, xmap)
        #a and b are two nodes that ara connected by sv edge
        a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type1)
        b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type2)
        tuples_list.append((a,b))
        tuples_list.append((b,a))
    sv_tuples_set = set(tuples_list)
    return sv_tuples_set

def associate_segments_to_svs(paths,g,sv_tuples_set,node_to_smap_dict,centro):
    """
    Associate segments to structural variations (SVs) by converting Eulerian paths with vertex IDs to segment paths.

    Parameters:
    - paths (list of lists): A list of Eulerian paths where each path is represented as a list of vertices.
    - g (object): A graph object representing the underlying graph.
    - sv_tuples_set (set): A set of tuple pairs, where each tuple (a, b) represents two nodes in the graph 
                           connected by an SV edge.
    - node_to_smap_dict (dict): A dictionary mapping nodes to their corresponding SMAP IDs.

    Returns:
    - list: A list of segments associated to SVs, where each segment is represented as a list 
            [Segment_Name, Path_Name, SMAP_IDs].

    Internal Functions and Steps:
    1. return_all_edges_in_cc: Returns all edges in a connected component of the graph.
    2. detect_segment_vertices: Identifies the vertices that belong to segments in the graph.
    3. check_non_centromeric_path: Checks if a path is non-centromeric.
    4. return_mapids: Returns the SMAP IDs for a pair of nodes.

    Example:
    Given paths, g, sv_tuples_set, and node_to_smap_dict, the function will return a list of segments associated to SVs.

    Notes:
    - This function prints intermediate paths.
    - Assumes some functions (like return_all_edges_in_cc, detect_segment_vertices, check_non_centromeric_path, 
      return_mapids) exist in the surrounding context.
    """
    c = 1
    segs_list = []
    for p in paths:
        component = list(set(p))
        component_edges = return_all_edges_in_cc(component, g)
        segment_vertices = detect_segment_vertices(component, component_edges)
        ans = []
        temp = [p[0]]
        for i in range(1,len(p)-1):
            if p[i] in segment_vertices and p[i-1]==p[i+1]:
                temp.append(p[i])
                if check_non_centromeric_path(temp,g,centro):
                    ans.append(temp)
                temp = [p[i]]
            else:
                temp.append(p[i])
        temp.append(p[-1])
        if check_non_centromeric_path(temp,g,centro):
            ans.append(temp)
        ans2 = []
        for subpath in ans:
            temp = ''
            path = f"Path {c}"
            print(path)
            for i in range(0,len(subpath)-1,2):
                seg_number = int((max(subpath[i],subpath[i+1])+1)/2)
                try:
                    print((subpath[i+1],subpath[i+2]))
                    if (subpath[i],subpath[i+1]) in sv_tuples_set:
                        smap_ids = return_mapids(str(subpath[i]),str(subpath[i+1]),node_to_smap_dict)
                        segs_list.append([f"Segment {seg_number}", path, smap_ids])
                    if (subpath[i+1],subpath[i+2]) in sv_tuples_set:
                        smap_ids = return_mapids(str(subpath[i+1]),str(subpath[i+2]),node_to_smap_dict)
                        segs_list.append([f"Segment {seg_number}", path, smap_ids])
                except:
                    continue
            c += 1
    return segs_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-cnv", "--cnv", help="path to cnv call (cnv_call_exp.txt)", required=True)
    parser.add_argument("-smap", "--smap", help="path to smap file", required=True)
    parser.add_argument("-rcmap", "--rcmap", help="path to CNV rcmap file (cnv_rcmap_exp.txt)", required=True)
    parser.add_argument("-xmap", "--xmap", help="path to contig alignments file xmap", required=True)
    parser.add_argument("-centro", "--centro", help="path to file contains centromere coordinates", required=False)
    parser.add_argument("-cyto", "--cyto", help="path to file contains cytoband coordinates", required=False)
    parser.add_argument("-n", "--name", help="output name", required=True)
    parser.add_argument("-o", "--output", help="path to output dir", required=True)
    parser.add_argument("-sv", "--sv_output", help="path to output dir", required=True)
    args = parser.parse_args()
    if args.centro is not None: # this will parse centromere region. It can be hard coded.
        centro = parse_centro(args.centro)
    else:
        centro = None
    segments, all_seg = parse_cnvcall(args.cnv)
    sv_output = args.sv_output
    smap = parse_smap(args.smap)
    segments = merge_segments_all_seg_smap(segments, all_seg, smap, centro) # Need to debug this function
    segments.sort(key=lambda x: (int(x.chromosome), x.start))
    rcov, rcop = parse_rcmap(args.rcmap)
    chrY_cn = int(np.average(list(rcop['24'].values())) + 0.5)
    # chrX_cn = 2
    chrX_cn = round(np.average(list(rcop['23'].values())))
    if chrY_cn > 0:
        chrX_cn = 1
    xmap = parse_xmap(args.xmap)
    output = args.output+'/'+ args.name + '.txt'
    iscn_output = args.output+'/'+ args.name + '_ISCN' + '.txt'
    file = args.output+'/'+ args.name + '.png'
    name = args.name
    cytoband_filtered = read_in_cyto(args.cyto)
    svs = []
    segments = extend_segments_cn(segments, all_seg) #fill the gap between calls. 
    for k in rcop.keys():
        seg_list = []
        label_list = list(rcop[k].keys())
        for s in segments:
            if s.chromosome == k:
                if s.width > 200000: #if call has length greater than 200Kbp assume a segment
                    seg_list.append(s)
        prev_point = list(rcop[k].keys())[0]
        if len(seg_list) == 0: #create segment for start of chromosme
            new_seg = Segments()
            new_seg.start = 0
            new_seg.end = list(rcop[k].keys())[-1]
            new_seg.width = list(rcop[k].keys())[-1]
            new_seg.chromosome = k
            new_seg.fractional_cn = cn_in_mask_N_region(k,new_seg.start,new_seg.end, rcop[k],centro)
            new_seg.int_cn = cn_in_mask_N_region(k,new_seg.start,new_seg.end, rcop[k],centro) #assumption that sample is diploide. default CN = 2
            if int(k) == 23:
                new_seg.fractional_cn = chrX_cn
                new_seg.int_cn = chrX_cn
            if int(k) == 24:
                new_seg.fractional_cn = chrY_cn
                new_seg.int_cn = chrY_cn
            new_seg.bp = [0, list(rcop[k].keys())[-1]]
            segments.append(new_seg)
    for k in rcop.keys():
        seg_list = []
        label_list = list(rcop[k].keys())
        for s in segments:
            if s.chromosome == k:
                if s.width > 200000: #if call has length greater than 200Kbp assume a segment
                    seg_list.append(s)
        prev_point = list(rcop[k].keys())[0]
        if len(seg_list) == 0: #create segment for start of chromosme
            new_seg = Segments()
            new_seg.start = 0
            new_seg.end = list(rcop[k].keys())[-1]
            new_seg.width = list(rcop[k].keys())[-1]
            new_seg.chromosome = k
            new_seg.fractional_cn = cn_in_mask_N_region(k,new_seg.start,new_seg.end, rcop[k], centro)
            new_seg.int_cn = cn_in_mask_N_region(k,new_seg.start,new_seg.end, rcop[k], centro) #assumption that sample is diploide. default CN = 2
            if int(k) == 23:
                new_seg.fractional_cn = chrX_cn
                new_seg.int_cn = chrX_cn
            if int(k) == 24:
                new_seg.fractional_cn = chrY_cn
                new_seg.int_cn = chrY_cn
            new_seg.bp = [0, list(rcop[k].keys())[-1]]
            segments.append(new_seg)
        else:
            seg_list.sort(key=lambda x: x.start)
            for s in seg_list:
                start, end = find_start_end(prev_point, s.start, label_list) # there are labels between two segments then create segment with CN =2 between them
                if start != 0 and end != 0: #
                    new_seg = Segments()
                    new_seg.start = start
                    new_seg.end = end
                    new_seg.width = end - start
                    new_seg.chromosome = k
                    new_seg.fractional_cn = cn_in_mask_N_region(k,new_seg.start,new_seg.end, rcop[k], centro)
                    new_seg.int_cn = cn_in_mask_N_region(k,new_seg.start,new_seg.end, rcop[k], centro) #assumption that sample is diploide. default CN = 2
                    if int(k) == 23:
                        new_seg.fractional_cn = chrX_cn
                        new_seg.int_cn = chrX_cn
                    if int(k) == 24:
                        new_seg.fractional_cn = chrY_cn
                        new_seg.int_cn = chrY_cn
                    new_seg.bp = [start, end]
                    segments.append(new_seg)
                prev_point = s.end
            start, end = find_start_end(s.end, label_list[-1], label_list)
            if start != 0 and end != 0:
                new_seg = Segments()
                new_seg.start = start
                new_seg.end = end
                new_seg.width = end - start
                new_seg.chromosome = k
                new_seg.fractional_cn = cn_in_mask_N_region(k,new_seg.start,new_seg.end, rcop[k], centro)
                new_seg.int_cn = cn_in_mask_N_region(k,new_seg.start,new_seg.end, rcop[k], centro)
                if int(k) == 23:
                    new_seg.fractional_cn = chrX_cn
                    new_seg.int_cn = chrX_cn
                if int(k) == 24:
                    new_seg.fractional_cn = chrY_cn
                    new_seg.int_cn = chrY_cn
                new_seg.bp = [start, end]
                segments.append(new_seg)

    segments.sort(key=lambda x: (int(x.chromosome), x.start))
    # for s in segments:
    #     print('asli', s.chromosome, s.start, s.end, s.int_cn, sorted(s.bp))
    for i in smap:
        # translocation applied filters.
        if i.sv_type.startswith('trans') and i.confidence >= 0.05 and not i.sv_type.endswith(
                'segdupe') and not i.sv_type.endswith('common') and not i.sv_type.endswith('oveerlap') and (i.ref_c_id1!= i.ref_c_id2 or abs(i.ref_end - i.ref_start) > 300000):
            svs.append(i)
            exist,s = detect_receprical_translocation(i, xmap, smap)
            if exist:
                svs.append(s)
                print(s.line.strip())
            print(i.line.strip())
        # indels
        elif i.sv_type.startswith('inse') or i.sv_type.startswith('delet'):
        # elif i.sv_type.startswith('delet'):
            if not i.sv_type.endswith('nbase') and not i.sv_type.endswith('tiny') and i.confidence >= 0:
                if i.sv_type.startswith('delet') and i.size > 200000:
                    print(i.line.strip())
                    if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)[0]:
                        _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)
                        svs.append(i)
                        print(i.line.strip())
                elif i.size > 500000 and abs(i.ref_end - i.ref_start) > 500000: # this would be for insertion length more than 500Kbp
                    svs.append(i)
                    print(i.line.strip())
        #if we have inversion SV
        elif i.sv_type == 'inversion' and i.confidence >= 0.7: # filter low confidance
            start, end = 0, 0
            dir = ''
            dir1, dir2 = detect_sv_directions(i, xmap)
            s = find_in_smap(i.linkID, smap) #inversion has two rwo in smap file. we find them with Link ID
            if dir1 == 'H': #update inversion call. it is to complicated but baisically calculate inversion start and endpoint
                dir = 'left'
                start = s.ref_start
                end = i.ref_end
            else:
                dir = 'right'
                start = i.ref_start
                end = s.ref_start
            start, end = min(start, end), max(start, end)
            i.ref_start = start
            i.ref_end = end
            if abs(end - start) > 800000: #apply filter on size of inversion
                svs.append(i)
                print(i.line.strip(), start, end,dir)
                print(s.line.strip())
        elif i.sv_type == 'inversion_paired' and i.confidence >= 0.7: #if it is full inversion
            s = find_in_smap(i.linkID, smap)
            if abs(i.ref_start - s.ref_end) > 500000:
                i.ref_end = s.ref_end
                print(i.line.strip())
                svs.append(i)
        elif i.sv_type == 'duplication_inverted':
            if detect_duplicatioon_inversion_cn(i, xmap, segments)[0]:
                _, fold_point = detect_duplicatioon_inversion_cn(i, xmap, segments) #because CN is changed
                i.ref_start = fold_point
                i.ref_end = fold_point
                svs.append(i)
                print(i.line.strip())
        elif i.sv_type == 'duplication' or i.sv_type == 'duplication_split':
            if detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)[0]:
                _, i.ref_start, i.ref_end = detect_del_dup_cn(i.ref_c_id1, i.ref_start, i.ref_end, segments)
                svs.append(i)
                print(i.line.strip())
        elif i.size > 500000 and not i.sv_type.startswith('inversion'): #Other type of SV
            svs.append(i)
            print(i.line.strip())

    for sv in svs:#integrate BPs and Segments
        find_bp_in_segment(sv.ref_c_id1, sv.ref_start, segments) #
        find_bp_in_segment(sv.ref_c_id2, sv.ref_end, segments)
    #merging Bps for spiliting segments 
    for s in segments:
        s.bp = merge_list(s.bp)

    #Graph Creation
    aa = 0
    g = Graph()
    counter = 0
    prev_chr = 0
    for index, s in enumerate(segments): #Forach segment creates two vertices
        for inedx_bp, i in enumerate(s.bp):
            if inedx_bp == 0: #if it is first BP in a segment
                aa += 1
                v = Vertices() #Create Vetrices Object
                v.chromosome = s.chromosome
                v.id = counter
                v.pos = i
                v.cn = s.int_cn
                v.type = 'H' # the start node is Head node
                g.vertices.append(v)
                if prev_chr == s.chromosome:
                    g.edges.append((counter - 1, counter, 0, 'R'))  #addign reference edge between two continues segments
                    g.return_node(counter - 1).append_edges(counter) #update adjancency matrix
                    g.return_node(counter).append_edges(counter - 1)
                prev_chr = s.chromosome
                counter += 1
            elif inedx_bp == len(s.bp) - 1: #if it is last BP in a segment
                aa += 1
                v = Vertices()
                v.chromosome = s.chromosome
                v.id = counter
                v.pos = i
                v.cn = s.int_cn
                v.type = 'T' #last node should be Tail node
                g.vertices.append(v)
                g.edges.append((counter - 1, counter, s.int_cn, 'S'))
                g.return_node(counter).append_edges(counter - 1)
                g.return_node(counter - 1).append_edges(counter)
                # g.edges.append((counter,counter-1, s.int_cn,'S'))
                counter += 1
                prev_chr = s.chromosome
            else: #Create both Tail and Head node
                aa += 1
                v = Vertices()
                v.chromosome = s.chromosome
                v.id = counter
                v.pos = i
                v.cn = s.int_cn
                v.type = 'T'
                g.vertices.append(v)
                g.edges.append((counter, counter - 1, s.int_cn, 'S'))
                g.return_node(counter - 1).append_edges(counter)
                g.return_node(counter).append_edges(counter - 1)
                # g.edges.append((counter-1,counter, s.int_cn,'S'))
                counter += 1
                ####
                aa += 1
                v = Vertices()
                v.chromosome = s.chromosome
                v.id = counter
                v.pos = i + 1
                v.type = 'H'
                v.cn = s.int_cn
                g.vertices.append(v)
                g.edges.append((counter - 1, counter, 0, 'R'))
                g.return_node(counter - 1).append_edges(counter)
                g.return_node(counter).append_edges(counter - 1)
                counter += 1
    for sv in svs:
        n_type1, n_type2 = detect_sv_directions(sv, xmap)
        #a and b are two nodes that ara connected by sv edge
        a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type1)
        b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type2)
        if int(sv.ref_c_id1) > int(sv.ref_c_id2) or (int(sv.ref_c_id1) == int(sv.ref_c_id2) and sv.ref_end < sv.ref_start):
            a = find_nodes(sv.ref_c_id1, sv.ref_start, g.vertices, n_type2)
            b = find_nodes(sv.ref_c_id2, sv.ref_end, g.vertices, n_type1)
        # if sv.sv_type == 'inversion' or sv.sv_type == 'inversion_paired':
        if sv.sv_type == 'inversion_paired': #Lets complete the inversion
            if g.return_node(a).type == 'H':
                new_edge = (a-1,b-1,0,'SV')
            else:
                new_edge = (a+1,b+1,0,'SV')
            if new_edge not in g.edges:
                    g.edges.append(new_edge)
                    g.return_node(new_edge[0]).append_edges(new_edge[1])
                    g.return_node(new_edge[1]).append_edges(new_edge[0])
        if b == a and (a, b, 0, 'SV') not in g.edges: #it can be happend in duplication inversion
            g.edges.append((a, b, 0, 'SV'))
            g.return_node(a).append_edges(b)
        #     print(sv.line)
        #     print(a,b)
        elif sv.sv_type.startswith('delet'): #telomere cite deletion prevent
            if (a, b, 0, 'SV') not in g.edges and abs(a-b)!=1:
                g.edges.append((a, b, 0, 'SV'))
                g.return_node(a).append_edges(b)
                g.return_node(b).append_edges(a)
        elif (a, b, 0, 'SV') not in g.edges:
            g.edges.append((a, b, 0, 'SV'))
            g.return_node(a).append_edges(b)
            g.return_node(b).append_edges(a)
        # g.edges.append((b,a,0,'SV'))
    g.print_node()
    print(g.edges)
    Plot_graph(g,file,name,centro)
    connected_components = find_connected_components(g)
    for component in connected_components:
    # if 6 in component:
        component_edges = estimating_edge_multiplicities_in_CC(component, g, xmap)
    connected_components = find_connected_components(g)
    paths = []
    for component in connected_components:
        component_edges = return_all_edges_in_cc(component, g)
        print(component)
        print(component_edges)
        paths.append(printEulerTour(component, component_edges, g, centro))
    node_to_map_dict, node_to_smap_dict = node_to_map(svs, xmap, g)
    path_map = {}
    smap_frames = []

    with open(output , 'w') as f :
        f.write('Segment\tNumber\tChromosome\tStart\tEnd\tStartNode\tEndNode\tMapIDs\tSmapIDs\n')
        number = 1
        for i in range(0,len(g.vertices),2):
            v = g.vertices[i]
            u = g.vertices[i+1]
            mapids = return_mapids(v.id,u.id,node_to_map_dict)
            smapids = return_mapids(v.id,u.id,node_to_smap_dict)
            f.write('Segment\t{id}\t{chrom}\t{start}\t{end}\t{snode}\t{enode}\t{mapids}\t{smapids}\n'.format(id = number, chrom = v.chromosome, start= v.pos, end= u.pos, snode = v.id, enode = u.id, mapids=mapids, smapids=smapids))
            p_n = str(number)
            smap_frames.append(smap_to_segment(p_n, smapids, smap))
            if p_n not in path_map:
                path_map[p_n] = 'chr{}:{}-{}'.format(v.chromosome,v.pos,u.pos)
            number += 1
        c = 1
        for p in paths:
            print(p)
            # structures, scores in convert_path_to_segment(p,g,centro)
            structures, scores = convert_path_to_segment(p,g,centro)
            print(structures)
            for jj in range(len(structures)):
            # for structure,scores in convert_path_to_segment(p,g,centro):
                structure = structures[jj]
                merged_coords,iscn_coords = convert_path(structure, path_map, cytoband_filtered)
                f.write('Path'+str(c)+ ' = '+structure+'\n')
                f.write('Path'+str(c)+ ' = '+merged_coords+'\n')
                f.write('Path'+str(c)+ ' = '+iscn_coords+'\n')
                c+=1
    sv_tuples_set = find_sv_node_edges(svs, xmap, g)
    segs_list = associate_segments_to_svs(paths, g ,sv_tuples_set, node_to_smap_dict,centro)
    subset = [x for x in smap_frames if isinstance(x,pd.DataFrame)]
    subset_smap_frame = pd.concat(subset)
    subset_smap_frame['Paths']= subset_smap_frame.apply(lambda x: [], axis=1)
    subset_smap_frame.index.name = 'Segments'
    subset_smap_frame.reset_index(inplace=True)
    for segment in segs_list:
        matching_rows = subset_smap_frame[subset_smap_frame['Segments'] == segment[0]]
        if not matching_rows.empty:
            idx = matching_rows.index[0]
            subset_smap_frame.at[idx,'Paths'].append(segment[1])
    node_to_map_dict,_ = node_to_map(svs, xmap, g)
    path_map = {}
    with open(iscn_output, 'w') as f :
        number = 1
        for i in range(0,len(g.vertices),2):
            v = g.vertices[i]
            u = g.vertices[i+1]
            mapids = return_mapids(v.id,u.id,node_to_map_dict)
            p_n = str(number)
            if p_n not in path_map:
                path_map[p_n] = 'chr{}:{}-{}'.format(v.chromosome,v.pos,u.pos)
            number += 1
        c = 1
        for p in paths:
            structures,scores = convert_path_to_segment(p,g,centro)
            # for structure,scores in convert_path_to_segment(p,g,centro):
            for jj in range(len(structures)):
                structure = structures[jj]
                split_structure = structure.split()
                segments_list = ["Segment {}".format(x.replace('-','').replace('+','')) for x in split_structure]
                path = f"Path {c}"
                merged_coords,iscn_coords = convert_path(structure, path_map, cytoband_filtered)
                f.write('Path'+str(c)+ ' = '+iscn_coords+'\n')
                c+=1
    subset_smap_frame.to_csv(sv_output,index=False)


if __name__ == "__main__":
    main()