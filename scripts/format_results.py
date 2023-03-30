import numpy as np
import pandas as pd
import argparse
import re
pd.set_option('mode.chained_assignment', None)


def format_iscn_results(iscn_format, cytobands, cytobands_out, genome_out, contig_orientation_out, kprect_out):
    """
    Formats ISCN Virtual Karyotype results into KaryoploteR format

    Parameters
    ----------
    iscn_format: str, required
        {sample}_results_ISCN.txt generated by rule: run_vk
    cytobands: str, required
        Cytoband file located in resources folder (cytoBand.txt)
        columns: chr, start, end, band, gieStain
    cytobands_out: str, required
        relative path to sample specific KaryoploteR formatted cytoband file
    genome_out: str, required
        relative path to sample specific KaryoploteR formatted genome file
    contig_orientation_out: str, required
        relative path to sample specific contig orientation file 
    kprect_out: str, required
        relative path to sample specific kprect parameters file       
    """
    cytoband_filtered = read_in_custom_cyto(cytobands = cytobands)
    iscn_results = pd.read_table(iscn_format,skiprows=1,sep=' = ',engine='python',header=None)
    iscn_results.columns = ['paths','iscn']
    iscn_split_coords = iscn_results['iscn'].str.split(',')
    iscn_paths = iscn_results['paths']
    cytoband_paths = []
    for path,coord in zip(iscn_paths,iscn_split_coords):
        for seg in coord:
            chrom_seg = seg.split('(')[0]
            chrom = return_chrom(chrom_seg)
            start,end = seg.split('(')[1].split(')')[0].split('-')
            strand = seg.split('(')[2].split(')')[0]
            sub_cyto = map_cyto_coords(chrom, start, end, strand, cytoband_filtered, path)
            cytoband_paths.append(sub_cyto)
    iscn_cyto = pd.concat(cytoband_paths)
    iscn_cyto_grouped = iscn_cyto.groupby('name')
    intervals_out_list = []
    for path,sub_cyto_frame in iscn_cyto_grouped:
        intervals = check_paths(sub_cyto_frame, path)
        intervals_out_list.append(intervals)
    interval_out = pd.concat(intervals_out_list)
    iscn_cyto['chrom'] = iscn_cyto['chr']
    iscn_cyto['chr'] = iscn_cyto['name']
    path_map_dict = interval_out.loc[:,['path','contig_id']].drop_duplicates().set_index('path').to_dict()['contig_id']
    update_mapper = {k:'{} : {}'.format(k,v) for k,v in path_map_dict.items()}
    iscn_cyto['chr'] = iscn_cyto['chr'].map(update_mapper)
    iscn_cyto.rename(columns={'stain':'gieStain','name':'Path','band':'name'},inplace=True)
    resolved_cyto_non_orient, resolved_stain = resolve_iscn_cyto(iscn_cyto)
    resolved_cyto = orient_p_to_q(resolved_cyto_non_orient)
    resolved_cyto.rename(columns={'gieStain':'original_stain','Resolved_stain':'gieStain','name':'original_name','Resolved_band':'name'},inplace=True)
    resolved_cyto.to_csv(cytobands_out, sep='\t',index=False)
    resolved_cyto['delta'] = resolved_cyto['end'] - resolved_cyto['start']
    path_order = resolved_cyto['chr'].unique().tolist()
    custom_genome = resolved_cyto.loc[:,['chr','delta']].groupby('chr').sum().sort_values('chr').reset_index()
    custom_genome.chr = custom_genome.chr.astype('category')
    custom_genome.chr = custom_genome.chr.cat.set_categories(path_order)
    custom_genome['start'] = 0
    custom_genome['end'] = custom_genome['delta']
    custom_genome_out = custom_genome.reindex(['chr','start','end'],axis=1).sort_values('chr')
    custom_genome_out.to_csv(genome_out,sep='\t',index=False)
    generate_karyoploter_plotting_params(resolved_cyto, contig_orientation_out, kprect_out)

def generate_karyoploter_plotting_params(resolved_cyto, contig_orientation_out, kprect_out):
    """
    Generates both kprect and contig_orientation plotting parameters to be used during figure generation step

    Parameters
    ----------
    resolved_cyto: DataFrame, required
        DataFrame of resolved ISCN cytoband embedded paths generated by function: resolve_iscn_cyto()
    contig_orientation_out : str, required
        handle for output contig orientation file 
    kprect_out: str, required
        handle for output contig orientation file
    """
    contig_orientation_frame = generate_orientation_arrow_coordinates(resolved_cyto)
    kprect_frame = generate_kprect_coordinates(resolved_cyto)
    contig_orientation_frame.to_csv(contig_orientation_out,sep='\t',index=False)
    kprect_frame.to_csv(kprect_out,sep='\t',index=False)

def swap_negative_orientation(collapsed_segments):
    """
    Swaps start and end positions for segments of orientation file with negative orientation

    Parameters
    ----------
    collapsed_segments: DataFrame, required
        DataFrame of collapsed orientation segments generated by function: generate_orientation_arrow_coordinates()
    Returns
    -------
    collapsed_segments: DataFrame
        Orientation resolved collapsed_segments Dataframe
    """
    idx = collapsed_segments['strand'] == '-'
    collapsed_segments.loc[idx,['start','end']] = collapsed_segments.loc[idx,['end','start']].values
    return collapsed_segments

def generate_orientation_arrow_coordinates(resolved_cyto):
    """
    Generates orienation arrow plotting parameter file

    Parameters
    ----------
    resolved_cyto: DataFrame, required
        DataFrame of collapsed orientation segments generated by function: resolved_cyto()
    Returns
    -------
    collapsed_segments: DataFrame
        Orientation resolved collapsed_segments Dataframe
    """
    collapsed_segment_list = []
    grouped_iscn = resolved_cyto.groupby('Path',sort=False)
    for g,frame in grouped_iscn:
        frame = frame.reset_index(drop=True)
        indx_final = frame.index.max()
        for indx,row in frame.iterrows():
            if indx == 0:
                segment_start = row
            if (segment_start.strand != row.strand) | (segment_start.chrom != row.chrom):
                append_row = frame.iloc[indx-1]
                append_row.start = segment_start.start
                collapsed_segment_list.append(append_row.to_frame().T)
                segment_start = row
            if (indx == indx_final) & (segment_start == row).all():
                collapsed_segment_list.append(row.to_frame().T)
            if (indx == indx_final) & (segment_start != row).any():
                append_row = row
                append_row.start = segment_start.start
                collapsed_segment_list.append(append_row.to_frame().T)
    collapsed_segments= pd.concat(collapsed_segment_list).sort_values(['Path','start','end'])
    map_dict = {'+':'blue','-':'red'}
    collapsed_segments['color'] = collapsed_segments['strand'].map(map_dict)
    collapsed_segments = swap_negative_orientation(collapsed_segments)
    return collapsed_segments

def generate_kprect_coordinates(resolved_cyto):
    """
    Generates kprect - contig specific plotting parameter file

    Parameters
    ----------
    resolved_cyto: DataFrame, required
        DataFrame of collapsed orientation segments generated by function: resolved_cyto()
    Returns
    -------
    collapsed_segments: DataFrame
        Orientation resolved collapsed_segments Dataframe
    """
    collapsed_segment_list = []
    grouped_iscn = resolved_cyto.groupby('Path',sort=False)
    for g,frame in grouped_iscn:
        frame = frame.reset_index(drop=True)
        indx_final = frame.index.max()
        for indx,row in frame.iterrows():
            if indx == 0:
                segment_start = row
            if (segment_start.chrom != row.chrom):
                append_row = frame.iloc[indx-1]
                append_row.start = segment_start.start
                collapsed_segment_list.append(append_row.to_frame().T)
                segment_start = row
            if (indx == indx_final) & (segment_start == row).all():
                collapsed_segment_list.append(row.to_frame().T)
            if (indx == indx_final) & (segment_start != row).any():
                append_row = row
                append_row.start = segment_start.start
                collapsed_segment_list.append(append_row.to_frame().T)
    collapsed_segments= pd.concat(collapsed_segment_list).sort_values(['Path','start','end'])
    collapsed_segments['midpoint'] = collapsed_segments.loc[:,['start','end']].mean(axis=1)
    return collapsed_segments

def resolve_iscn_cyto(iscn_cyto):
    """
    Resolves non-contiguous genomic segments and embeds major cytobands in paths detected by VK algorithm

    Parameters
    ----------
    iscn_cyto: DataFrame, required
        DataFrame of ISCN coordinates generated within the function: format_iscn_results()
    Returns
    -------
    resolved_stain_collapsed_cyto: DataFrame
        Orientation resolved collapsed_segments Dataframe generated by collapse_bands()
    resolved_cyto: DataFrame
        Orientation resolved collapsed_segments Dataframe
    resolved_stain_cyto: DataFrame
        Orientation resolved collapsed_segments Dataframe generated by resolve_iscn_bands()
    """
    iscn_cyto_copy = iscn_cyto.copy()
    cols = iscn_cyto_copy.columns.tolist()
    cols = cols + ['bp','original_start', 'original_end']
    iscn_cyto_copy['bp'] = False
    iscn_cyto_copy['original_start'] = iscn_cyto_copy['start']
    iscn_cyto_copy['original_end'] = iscn_cyto_copy['end']
    iscn_cyto_copy = iscn_cyto_copy.reindex(cols,axis=1)
    cyto_list = []
    cols = iscn_cyto_copy.columns.tolist()
    grouped_iscn = iscn_cyto_copy.groupby('Path')
    for g,frame in grouped_iscn:
        frame = frame.reset_index(drop=True)
        for i,row in frame.iterrows():
            if i == 0:
                if row.start != 0:
                    new_end = (row.end - row.start) 
                    frame.iloc[i] = [row.chr, 0, new_end, row['name'], row.gieStain, row.band_base, row.Updated_stain, row.Resolved_band, row.Resolved_stain, row.Path, row.strand, row.chrom, False, row.original_start, row.original_end]
                else:
                    frame.iloc[i] = row.values.tolist()
            else:
                previous_frame = frame.iloc[i-1]
                if frame.index.max == i:
                    new_end = (row.end - row.start) + previous_frame.end
                    frame.iloc[i] = [row.chr,previous_frame.end, new_end, row['name'], row.gieStain, row.band_base, row.Updated_stain, row.Resolved_band, row.Resolved_stain, row.Path, row.strand, row.chrom, False, row.original_start, row.original_end]
                else:
                    if row.start == previous_frame.end:
                        frame.iloc[i] = row.values.tolist()
                    else:
                        new_end = (row.end - row.start) + previous_frame.end
                        frame.iloc[i] = [row.chr,previous_frame.end,new_end,row['name'],row.gieStain,row.band_base, row.Updated_stain, row.Resolved_band, row.Resolved_stain, row.Path, row.strand, row.chrom, True, row.original_start, row.original_end]
        cyto_list.append(frame)               
    resolved_cyto = pd.concat(cyto_list)
    resolved_stain_cyto = resolve_iscn_bands(resolved_cyto)
    resolved_stain_collapsed_cyto = collapse_bands(resolved_stain_cyto)
    return resolved_stain_collapsed_cyto, resolved_cyto

def resolve_iscn_bands(resolved_cyto):
    """
    Selects gieStain based on prevelance in segments and retains centromere stain when present

    Parameters
    ----------
    resolved_cyto: DataFrame, required
        DataFrame of ISCN coordinates generated within the function: format_iscn_results()
    Returns
    -------
    resolved_stain_cyto: DataFrame
        Orientation resolved, minor cytoband coordinate collapsed_segments Dataframe
    """
    resolved_stain_list = []
    cols = resolved_cyto.columns.tolist()
    grouped_iscn = resolved_cyto.groupby(['Path','chrom','strand'],sort=False)
    for (g,c,s),frame in grouped_iscn:
        frame = frame.reset_index(drop=True)
        band_base_group = frame.groupby('band_base')
        for b,band_frame in band_base_group:
            stain_counts = band_frame['gieStain'].value_counts()
            if 'acen' in stain_counts:
                band_frame['Updated_stain'] = 'acen'
            else:
                band_frame['Updated_stain'] = stain_counts.index[0]
            resolved_stain_list.append(band_frame)
    cols.append('Updated_stain')
    resolved_stain_cyto = pd.concat(resolved_stain_list).sort_values(['Path','start','end'])
    return resolved_stain_cyto

def collapse_bands(resolved_stain_cyto):
    """
    Collapses minor cytobands into major cytobands

    Parameters
    ----------
    resolved_stain_cyto: DataFrame, required
        DataFrame of ISCN coordinates generated within the function: resolve_iscn_bands()
    Returns
    -------
    collapsed_bands: DataFrame
        Orientation resolved, minor cytoband coordinate collapsed_segments Dataframe
    """
    collapsed_band_list = []
    grouped_iscn = resolved_stain_cyto.groupby('Path',sort=False)
    for g,frame in grouped_iscn:
        frame = frame.reset_index(drop=True)
        indx_final = frame.index.max()
        for indx,row in frame.iterrows():
            if indx == 0:
                segment_start = row
            if (segment_start.chrom != row.chrom) or (segment_start.strand != row.strand) or (segment_start.Resolved_band!= row.Resolved_band):
                append_row = frame.iloc[indx-1]
                append_row.start = segment_start.start
                collapsed_band_list.append(append_row.to_frame().T)
                segment_start = row
            if (indx == indx_final) & (segment_start == row).all():
                collapsed_band_list.append(row.to_frame().T)
            if (indx == indx_final) & (segment_start != row).any():
                append_row = row
                append_row.start = segment_start.start
                collapsed_band_list.append(append_row.to_frame().T)
    collapsed_bands = pd.concat(collapsed_band_list).sort_values(['Path','start','end'])
    return collapsed_bands

def return_chrom(chrom):
    """
    Iterates through string ISCN chromosome and returns chromsome with ISCN nomenclature removed

    Parameters
    ----------
    chrom: str, required
        string ISCN chromosome
    Returns
    -------
    chrom_num: str
        chromsome with ISCN nomenclature removed
    """
    i=1
    for i in range(len(chrom)):
        if chrom[:i].isdigit():
            chrom_num = str(chrom[:i])
        else:
            continue
        i+=1
    return chrom_num

def find_intervals(sub_cyto_frame):
    """
    Finds gaps between segments and reports intervals

    Parameters
    ----------
    sub_cyto_frame: DataFrame, required
        Grouped by Path subsetted cytoband DataFrame
    Returns
    -------
    intervals: DataFrame
        DataFrame with gap filled intervals
    """
    # Perform this on chrom specific segments
    endpoints = sub_cyto_frame.loc[:,['start','end']].stack().sort_values().reset_index(drop=True)
    intervals = pd.DataFrame({'start':endpoints.shift().fillna(0), 
                            'end':endpoints}).astype(int)
    # construct the list of intervals from the endpoints
    intervals['intv'] = [pd.Interval(a,b) for a,b in zip(intervals.start, intervals.end)]
    # these are the original intervals
    orig_invt = pd.arrays.IntervalArray([pd.Interval(a,b) for a,b in zip(sub_cyto_frame.start, sub_cyto_frame.end)])
    # walk through the intervals and compute the intersections
    intervals['total'] = intervals.intv.apply(lambda x: orig_invt.overlaps(x).sum())
    intervals = intervals[~(intervals['start']==intervals['end'])]
    return intervals

def check_paths(sub_cyto_frame, path):
    """
    Function checks paths and identifies potential derivative chromosomes

    Parameters
    ----------
    sub_cyto_frame: DataFrame, required
        Grouped by Path subsetted cytoband DataFrame
    path: str, required
        Path identified by Virtual Karyotype algorithm
    Returns
    -------
    path_intervals: DataFrame
        DataFrame with gap filled intervals
    """
    sub_cyto_sorted = sub_cyto_frame.sort_values(['chr','start','end'])
    if len(sub_cyto_sorted['chr'].unique()) == 1:
        path_intervals = find_intervals(sub_cyto_sorted)
        path_intervals['chr'] = sub_cyto_sorted['chr'].unique()[0]
        contig_id = "{chrom}".format(chrom=sub_cyto_sorted['chr'].unique()[0])
    else:
        intervals_list = []
        sub_cyto_contigs = sub_cyto_sorted.groupby('chr')
        contig_id = "({})".format(';'.join(sub_cyto_sorted['chr'].value_counts().index))
        for chrom,sub_cyto in sub_cyto_contigs:
            sub_path_intervals = find_intervals(sub_cyto)
            sub_path_intervals['chr'] = chrom
            intervals_list.append(sub_path_intervals)
        path_intervals = pd.concat(intervals_list)
    path_intervals['path'] = path
    path_intervals['contig_id'] = contig_id
    return path_intervals

def map_cyto_coords(chrom, start, end, strand, cytoband_filtered, path):
    """
    Maps cytoband coordinates to genomic segments identified by VK algorithm
    
    Parameters
    ----------
    chrom: str, required
        chromosome to subset cytoband_filtered by
    start: str, required
        starting coordinates to subset cytoband_filtered by
    end: str, required
        ending coordinates to subset cytoband_filtered by
    strand: str, required
        if the segment is negative, reverse the ordering of the cytoband segment
    cytoband_filtered: Dataframe, required
        cytoband_filtered Dataframe load by function: read_in_cyto()
    path: str, required
        path identified by VK algorithm
    Returns
    -------
    sub_cyto: DataFrame
        subsetted cytoband Dataframe
    """
    start = int(start)
    end = int(end)
    sub_cyto = cytoband_filtered[(cytoband_filtered['chr'] == chrom) & ((cytoband_filtered['start'] >= start) | (cytoband_filtered['end'] >= start)) & ((cytoband_filtered['end'] <= end)|(cytoband_filtered['start'] <= end))]
    sub_cyto['name'] = path
    sub_cyto.iloc[0,1] = start
    sub_cyto.iloc[-1,2] = end
    sub_cyto['strand'] = strand
    if strand == '-':
        sub_cyto = sub_cyto.iloc[::-1]
    return sub_cyto

def read_in_cyto(cytobands = 'resources/cytoBand.txt'):
    """
    Loads cytoband file as Dataframe
    
    Parameters
    ----------
    cytobands: str, required
        relative path to cytoBand.txt
    Returns
    -------
    cytoband_filtered: Dataframe
        filtered cytoband Dataframe
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

def read_in_custom_cyto(cytobands = 'resources/hg38_400_level_cytoband.tsv'):
    """
    Loads cytoband file as Dataframe
    
    Parameters
    ----------
    cytobands: str, required
        relative path to cytoBand.txt
    Returns
    -------
    cytoband_filtered: Dataframe
        filtered cytoband Dataframe
    """
    cytoband = pd.read_table(cytobands)
    cytoband['chr'] = cytoband['chr'].astype(str)
    cytoband['start'] = cytoband['start'].astype(int)
    cytoband['end'] = cytoband['end'].astype(int)
    return cytoband

def check_strand(frame):
    if (len(frame['strand'].unique()) == 1) and (frame['strand'].unique()[0] == '-'):
        return True
 
def orient_p_to_q(resolved_cyto):
    """
    Resolves non-contiguous genomic segments and embeds major cytobands in paths detected by VK algorithm

    Parameters
    ----------
    resolved_cyto: DataFrame, required
        DataFrame of ISCN coordinates generated within the function: resolve_iscn_cyto()
    Returns
    -------
    resolved_strand_iscn_sorted: DataFrame
        sorted by presence of cetromere and ordered by autosome, iscn relabled cytoband Dataframe from relabel_paths()
    """
    iscn_cyto_copy = resolved_cyto.copy()
    cyto_list = []
    grouped_iscn = iscn_cyto_copy.groupby('Path')
    for g,frame in grouped_iscn:
        frame = frame.reset_index(drop=True)
        if check_strand(frame):
            frame = frame.iloc[::-1].reset_index(drop=True)
            for i,row in frame.iterrows():
                if i == 0:
                    if row.start != 0:
                        new_end = (row.end - row.start) 
                        frame.iloc[i] = [row.chr,0,new_end,row['name'],row.gieStain,row.band_base, row.Updated_stain, row.Resolved_band, row.Resolved_stain, row.Path, row.strand, row.chrom, row.bp, row.original_start, row.original_end]
                    else:
                        frame.iloc[i] = row.values.tolist()
                else:
                    previous_frame = frame.iloc[i-1]
                    if frame.index.max == i:
                        new_end = (row.end - row.start) + previous_frame.end
                        frame.iloc[i] = [row.chr,previous_frame.end,new_end,row['name'],row.gieStain,row.band_base, row.Updated_stain, row.Resolved_band, row.Resolved_stain, row.Path, row.strand, row.chrom, row.bp, row.original_start, row.original_end]
                    else:
                        if row.start == previous_frame.end:
                            frame.iloc[i] = row.values.tolist()
                        else:
                            new_end = (row.end - row.start) + previous_frame.end
                            frame.iloc[i] = [row.chr,previous_frame.end,new_end,row['name'],row.gieStain,row.band_base, row.Updated_stain, row.Resolved_band, row.Resolved_stain, row.Path, row.strand, row.chrom, row.bp, row.original_start, row.original_end]
            cyto_list.append(frame)
        else:
            cyto_list.append(frame)              
    resolved_strand_cyto = pd.concat(cyto_list)
    resolved_strand_iscn_sorted = relabel_paths(resolved_strand_cyto)
    return resolved_strand_iscn_sorted

def relabel_paths(resolved_strand_cyto):
    """
    Resolves non-contiguous genomic segments and embeds major cytobands in paths detected by VK algorithm
    Parameters
    ----------
    resolved_cyto: DataFrame, required
        DataFrame of ISCN coordinates generated within the function: format_iscn_results()
    Returns
    -------
    resolved_strand_iscn_sorted: DataFrame
        sorted by presence of cetromere and ordered by autosome, iscn relabled cytoband Dataframe
    """
    iscn_cyto_copy = resolved_strand_cyto.copy()
    cyto_list = []
    grouped_iscn = iscn_cyto_copy.groupby('Path')
    for g,frame in grouped_iscn:
        frame = frame.reset_index(drop=True)
        if frame['Updated_stain'].isin(['acen']).any():
            centromere = frame[frame['Updated_stain'] == 'acen']['chrom'].unique()
            if len(centromere) == 1:
                frame['centromere'] = 0
                frame['sort_chrom'] = centromere[0]
            else:
                frame['centromere'] = 1
                frame['sort_chrom'] = frame[frame['chrom'].isin(centromere)]['chrom'].value_counts().index[0]
            contigs = frame['chrom'].unique()
            if len(contigs) == 1:
                frame['iscn'] = "{path} : {contig}".format(path=g,contig=contigs[0])
            else:
                if len(centromere) == 1:
                    non_centromere = list(set(contigs) - set(centromere))
                    frame['iscn'] = "{path} : der({centromere})t({contig})".format(path=g,centromere=centromere[0],contig=';'.join(map(str,non_centromere)))
                else:
                    non_centromere = list(set(contigs) - set(centromere))
                    frame['iscn'] = "{path} : der({centromere})t({contig})".format(path=g,centromere=';'.join(map(str,centromere)),contig=';'.join(map(str,non_centromere)))
            cyto_list.append(frame)
        else:
            frame['centromere'] = 2
            contigs = frame['chrom'].unique()
            if len(contigs) == 1:
                frame['iscn'] = "{path} : {contig}".format(path=g,contig=contigs[0])
                frame['sort_chrom'] = contigs[0]
            else:
                frame['iscn'] = "{path} : {contig}".format(path=g,contig=';'.join(map(str,contigs)))
                frame['sort_chrom'] = frame[frame['chrom'].isin(contigs)]['chrom'].value_counts().index[0]
            cyto_list.append(frame)
    resolved_strand_iscn_cyto = pd.concat(cyto_list)
    index_cols = ['iscn', 'start', 'end', 'name', 'gieStain', 'band_base', 'Path', 'strand', 'chrom', 'bp', 'original_start', 'original_end', 'Updated_stain', 'centromere', 'chr','sort_chrom','Resolved_band','Resolved_stain']
    resolved_strand_iscn_cyto['sort_chrom'] = resolved_strand_iscn_cyto['sort_chrom'].astype(int)
    resolved_strand_iscn_sorted = resolved_strand_iscn_cyto.sort_values(['centromere','sort_chrom','chr','start','end']).reindex(index_cols,axis=1).rename(columns={'iscn':'chr','chr':'iscn'})
    return resolved_strand_iscn_sorted

def main():
    parser = argparse.ArgumentParser(
        """Function accepts generates formats Virtual Karyotype results to be visualized using KaryoploteR"""
    )
    parser.add_argument('--iscn_format', type=str, help="relative path to VK ISCN results file")
    parser.add_argument('--cytoband', type=str, help="relative path to cytoband file")
    parser.add_argument('--cytobands_out', type=str, help="relative path to granges formatted output file")
    parser.add_argument('--genome_out', type=str, help="relative path to granges formatted output file")
    parser.add_argument('--contig_orientation_out', type=str, help="relative path to granges formatted output file")
    parser.add_argument('--kprect_out', type=str, help="relative path to granges formatted output file")

    # parse command line arguments
    args = parser.parse_args()
    print(args)
    iscn_format = args.iscn_format
    cytoband = args.cytoband
    cytobands_out = args.cytobands_out
    genome_out = args.genome_out
    contig_orientation_out = args.contig_orientation_out
    kprect_out = args.kprect_out

    format_iscn_results(iscn_format=iscn_format, cytobands=cytoband, cytobands_out=cytobands_out, genome_out=genome_out, contig_orientation_out=contig_orientation_out, kprect_out=kprect_out)

if __name__ == "__main__":
    main()