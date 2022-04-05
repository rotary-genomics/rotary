#!/usr/bin/env python
# check circlator status
# Jackson M. Tsuji, 2022

import pandas as pd
import sys

input_circlator_logfile = 'merge.circularise.log' # -i

circlator_info = pd.read_csv(input_circlator_logfile, sep='\t')[['#Contig','circularised']]

if circlator_info['circularised'].sum() == 0:
    print('true')

elif circlator_info['circularised'].sum() > 0:
    print('false')

else:
    sys.exit('File processing error. # of non-circularized contigs is not 0 but is not > 0 either. Exiting...')
