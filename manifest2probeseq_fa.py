#! /usr/bin/env python

import os
import sys

if len( sys.argv ) != 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    sys.stderr.write( "Usage: %s <illumina manifest file>\n" % ( os.path.basename( sys.argv[0] ) ) )
    sys.exit()


start = False
fh = open( sys.argv[1] )
for line in fh:
    if 'IlmnID' in line:
        start = True
    elif 'Controls' in line:
        start = False
    elif start:
        line = line.rstrip( '\n' )
        if not line: continue

        fields = line.split( ',' )
        name = fields[1]
        snp  = fields[3].replace( '[', '' ).replace( ']', '' )
        probe_seq = fields[5]
        chrom = fields[9]
        pos   = fields[10]
        
        seq_id = '%s;%s;%s;%s' % ( name, snp, chrom, pos )
        print '>%s\n%s' % ( seq_id, probe_seq )
fh.close()
