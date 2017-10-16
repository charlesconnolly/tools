#! /usr/bin/env python

import os
import sys

import subprocess
import argparse

from collections import defaultdict 

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( '--hwe_cutoff', type=float, help='default: 0.0001',
                         default=0.001 )
    parser.add_argument( '--max_mcr', help='max missing call rate: default 0.02',
                         type=float, default=0.02 )
    parser.add_argument( '--keep_indels', action='store_true', help='default: exclude indels' )
    parser.add_argument( '--keep_monomorphic', action='store_true', help='exclude monomorphic' )
    parser.add_argument( '--plink_path', default='plink', 
                         help='default: version in current path' )
    parser.add_argument( 'bfile', help='plink bfile stem' )
    args = parser.parse_args()

    return args



args = parse_args()
bfile = os.path.basename( args.bfile )
log_file = bfile + '_qc.log'
log_fh = open( log_file, 'w' )

# - variant qc (do in the following order)
#   - remove hardy-weinberg < hwe_cutoff

log_fh.write( '\n##### HWE filtering #####\n' )
log_fh.flush()
sys.stderr.write( 'HWE filtering\n' )
hwe_bfile = bfile + '_hwe%f' % ( args.hwe_cutoff )
cmd = '%s --bfile %s ' % ( args.plink_path, args.bfile )
cmd += '--hwe %f ' % ( args.hwe_cutoff )
cmd += '--make-bed --out %s' % ( hwe_bfile )
subprocess.check_call( cmd, shell=True, stderr=log_fh, stdout=log_fh )

#  - remove variants with missing call rate > cutoff
log_fh.write( '\n##### filtering on MCR %f #####\n' % ( args.max_mcr ) )
log_fh.flush()
sys.stderr.write( 'filtering on MCR %f\n' % ( args.max_mcr ) )
mcr_bfile = hwe_bfile + '_mcr%f' % ( args.max_mcr )
cmd =  '%s --bfile %s ' % ( args.plink_path, hwe_bfile )
cmd += '--geno %f ' % ( args.max_mcr )
cmd += '--make-bed --out %s ' % ( mcr_bfile )
subprocess.check_call( cmd, shell=True, stderr=log_fh, stdout=log_fh )

#  - remove monomorphic
monomorphic_bfile = None
if not args.keep_monomorphic:
    log_fh.write( '\n##### removing momomorphic #####\n' )
    log_fh.flush()
    sys.stderr.write( 'removing momomorphic\n' )

    monomorphic_bfile = mcr_bfile + '_no_monomorphic'
    cmd = '%s --bfile %s ' % ( args.plink_path, mcr_bfile )
    cmd += '--mac 1 --make-bed --out %s ' % ( monomorphic_bfile )
    subprocess.check_call( cmd, shell=True, stderr=log_fh, stdout=log_fh )


# remove indels
# need to remove indels before removing dups, beacuse they can overlap

indels_bfile = None
if not args.keep_indels:
    bim_file_stem = monomorphic_bfile or mcr_bfile
    bim_file = bim_file_stem + '.bim'
    indels = []
    indel_alleles = set( [ 'D', 'I' ] )
    fh = open( bim_file )
    for line in fh:
        line = line.rstrip( '\n' )
        if not line: continue

        chrom, var_name, cm, pos, a1, a2 = line.split()
        if a1 in indel_alleles or a2 in indel_alleles:
            indels.append( var_name )
    fh.close()

    indels_file = bim_file_stem + '_indels.txt'
    indel_fh = open( indels_file, 'w' )
    indel_fh.write( '\n'.join( indels ) + '\n' )
    indel_fh.close()

    log_fh.write( '\n##### removing indels in %s #####\n' % ( indels_file ) )
    log_fh.flush()

    sys.stderr.write( 'removing indels in %s\n' % ( indels_file) )
    indels_bfile = bim_file_stem + '_no_indels'
    cmd = '%s --bfile %s ' % ( args.plink_path, bim_file_stem )
    cmd += '--exclude %s ' % ( indels_file )
    cmd += '--make-bed --out %s ' % ( indels_bfile )
    subprocess.check_call( cmd, shell=True, stderr=log_fh, stdout=log_fh )




# remove duplicates
bim_file_stem = indels_bfile or monomorphic_bfile or mcr_bfile
bim_file_stem + '.bim'
fh = open( bim_file )

chrom_pos2var_names = defaultdict( list )
for line in fh:
    line = line.rstrip( '\n' )
    if not line: continue

    chrom, var_name, cm, pos, a1, a2 = line.split()

    chrom_pos = '%s_%s' % ( chrom, pos )
    chrom_pos2var_names[chrom_pos].append( var_name )
fh.close()

dups = []
for chrom_pos, var_names in chrom_pos2var_names.iteritems():
    if len( var_names ) > 1:
        dups += var_names

if dups:
    dupsfile = bim_file_stem + '_dups.txt'
    out_fh = open( dupsfile, 'w' )
    out_fh.write( '\n'.join( dups ) + '\n' )
    out_fh.close()

    log_fh.write( '\n##### removing duplicate variants in %s #####\n' % ( dupsfile ) )
    log_fh.flush()
    sys.stderr.write( 'removing duplicate variants in %s\n' % ( dupsfile) )

    dedup_bfile = bim_file_stem + '_dedup'
    cmd = '%s --bfile %s ' % ( args.plink_path, bim_file_stem )
    cmd += '--exclude %s --make-bed --out %s' % ( dupsfile, dedup_bfile )
    subprocess.check_call( cmd, shell=True, stderr=log_fh, stdout=log_fh )




