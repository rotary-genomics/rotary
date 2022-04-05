#!/usr/bin/env python
# Get end regions of contigs

import pandas as pd

# Variables
assembly_info_filepath = '../../assembly/assembly_info.txt' # -i
circular_list_outfile = 'circular.list' # -a
end_regions_outfile = 'ends.bed' # -b
length_threshold = 100000 # -l
contig_id = 'contig_2' # -n

assembly_info = pd.read_csv('../../assembly/assembly_info.txt', sep='\t')
circular_contigs = assembly_info[assembly_info['circ.'] == 'Y'][['#seq_name','length','circ.']]
circular_contigs

# TODO - add multi-contig support
if contig_id is not None:
    circular_contigs = circular_contigs[circular_contigs['#seq_name'] == contig_id]

if circular_list_outfile is not None:
    circular_contigs['#seq_name'].to_csv(circular_list_outfile, index=False, header=None)

if end_regions_outfile is not None:
    contigs = []
    starts = []
    stops = []

    for index,row in circular_contigs.iterrows():
        contig,length,circ = row

        if length < length_threshold:
            contigs.append(contig)
            starts.append(0)
            stops.append(length)

        else:
            half_threshold = int(round(length_threshold / 2,0))

            contigs.append(contig)
            starts.append(0)
            stops.append(half_threshold)

            contigs.append(contig)
            starts.append(length - half_threshold)
            stops.append(length)

    end_regions = pd.DataFrame({'contig':contigs, 'start':starts, 'stop':stops})
    end_regions.to_csv(end_regions_outfile, sep='\t', header=None, index=False)
