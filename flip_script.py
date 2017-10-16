#! /usr/bin/env python

import os
import sys

import subprocess
import argparse
import pysam

from collections import defaultdict


# get ref allele for genomic coordinate
# param: string, chromosome
# param: string, position
# param: string, variant name
# param: argparse object
def get_ref_allele( chrom, pos, var_name, args ):
    
    # TBX: global vaiable representing pysam tabix object
    if TBX:
        pos = int( pos )
        start = pos - 1

        if args.chrom_as_int:
            if chrom == 'X': 
                chrom = '23'
            elif chrom == 'Y': 
                chrom = '24'
            elif chrom == 'M' or chrom == 'MT':
                chrom = '26'

        if args.chrom_as_N:
            chrom = chrom.replace( 'chr', '' )

        for row in TBX.fetch( chrom, start, pos ):
            chrom, pos, ref = row.split()

    # VAR2REF: global variable, dict where keys are variant names, values are ref alleles
    elif VAR2REF:
        ref = VAR2REF.get( var_name, None )
        
    return ref


# check out .bim, see if ambiguous alleles (A/T, C/G) are present
# param: string, filepath to bim file
# return: tuple of set of ambiguou alleles, set of non-ambiguous alleles, dict where keys are alleles, values are the number of each allele
def check_alleles( bim_file ):
    alleles2count = defaultdict( int )
    ambiguous = set()
    non_ambiguous = set()

    fh = open( bim_file )
    for line in fh:
        line = line.rstrip( '\n' )
        if not line: continue

        fields = line.split()
        var = fields[1]
        a1  = fields[4]
        a2  = fields[5]
        alleles = '%s\t%s' % ( a1, a2 )
        alleles2count[alleles] += 1

        if a1 in NTS and a2 in NTS and a1 == COMP[a2]:
            ambiguous.add( var )
        else:
            non_ambiguous.add( var )

    fh.close()

    return ( ambiguous, non_ambiguous, alleles2count )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( '--chrom_as_int', help='chrX/X->23, chrY/Y->24, chrM/chrMT/M/MT->26',
                         action='store_true' )
    parser.add_argument( '--chrom_as_N', action='store_true', help="chrom names lack 'chr'" )
    parser.add_argument( '--exclude_ambiguous', action='store_true' )
    parser.add_argument( '--flip_non_ambiguous', action='store_true' )
    parser.add_argument( '--flip_from_strand_file', action='store_true' )
    parser.add_argument( '--force', help='overwrite existing intermediary files', action='store_true' )
    parser.add_argument( '--strand_file', required=True )
    parser.add_argument( '--plink_path' )
    parser.add_argument( '--ref_genome_tabix_file' )
    parser.add_argument( '-v', '--var_to_ref_file', 
                         help='tab-separated var, ref_allele file' )
    parser.add_argument( 'bfile' )
    args = parser.parse_args()

    return args



# check if file exists and whether to overwrite it
# param: string, filepath
# param: boolena, whether to overwrite
# return: boolean, whether to make a file
def make_file( fname, force ):
    if os.path.isfile(fname) and not force:
        return False
    elif not os.path.isfile( fname ):
        return True


# make files for plink input mapping variant names to position, separately to chromosome
# param: argparse object
# actions: creates files
# return: tuple of strings representing file paths
def make_chr_pos_update_files( args ):
    var2chrom = {}
    var2pos   = {}

    chrom_file = os.path.basename( args.strand_file + '.var2chrom.tsv' )
    pos_file   = os.path.basename( args.strand_file + '.var2pos.tsv' )

    fh = open( args.strand_file )
    for line in fh:
        line = line.rstrip( '\n' )
        if not line: continue
        if line.startswith( '#' ): continue
        
        fields = line.split()
        var = fields[0]
        chrom = fields[1]
        pos   = fields[2]

        var2chrom[var] = chrom
        var2pos[var] = pos
    fh.close()

    log_fh.write( 'making update chromosome file: %s\n' % ( chrom_file ) )
    chr_fh = open( chrom_file, 'w' )
    for var, chrom in var2chrom.iteritems():
        chr_fh.write( '\t'.join( [ var, chrom ] ) + '\n' )
    chr_fh.close()

    log_fh.write( 'making update position file: %s\n' % ( chrom_file ) )
    pos_fh = open( pos_file, 'w' )
    for var, pos in var2pos.iteritems():
        pos_fh.write( '\t'.join( [ var, pos ] ) + '\n' )
    pos_fh.close()

    return ( chrom_file, pos_file )
    
# wrapper to plink command to update variant chromosome
# param: string, path to plink bfile stem
# param: string, path to file mapping variant name to chromosome
# param: string, path to plink
# param: string, log filehandle
# actions: calls plink --update_chr, creates output plink bedfile
# return: string, updated bfile stem
def update_chrom( bfile, chr_file, plink_path, log_fh ):
    sys.stderr.write( 'updating chromosome\n')
    log_fh.write( '\n##### updating chromosome #####\n' )
    log_fh.flush()
    plink_path = args.plink_path or 'plink'
    updated_chrom_bfile = bfile + '_updated_chrom'

    cmd = '%s --bfile %s ' % ( plink_path, args.bfile )
    cmd += '--update-chr %s  ' % ( chr_file )
    cmd += '--make-bed --out %s' % ( updated_chrom_bfile )
    plink_output = subprocess.check_call( cmd, shell=True, stdout=log_fh, stderr=log_fh )
    log_fh.flush()

    return updated_chrom_bfile

# wrapper to plink command to update variant position
# param: string, path to plink bfile stem
# param: string, path to file mapping variant name to position
# param: string, path to plink command
# param: log filehandle
# actions: calls plink --update_map, creates output plink bedfile
# return: string, updated bfile stem
def update_pos( bfile, pos_file, plink_path, log_fh ):
    sys.stderr.write( 'updating position \n' )
    log_fh.write( '\n##### updating position #####\n' )
    updated_pos_bfile = bfile + '_updated_pos'
    cmd = '%s --bfile %s ' % ( plink_path, bfile )
    cmd += '--update-map %s ' % ( pos_file )
    cmd += '--make-bed --out %s' % ( updated_pos_bfile )
    plink_output = subprocess.check_call( cmd, shell=True, stdout=log_fh, stderr=log_fh )
    log_fh.flush()

    return updated_pos_bfile

# wrapper to plink command to extract variants represented by well-aligned probesets
# param: string, path to plink bfile stem
# param: string, path to file with list of good variants
# param: string, path to plink command
# param: log filehandle
# actions: calls plink --extract, creates output plink bedfile
# return: string, updated bfile stem
def extract_aligned( bfile, good_var_file, plink_path, log_fh ):
    sys.stderr.write( 'extracting variants with good alignment\n' )
    log_fh.write( '\n##### extracting variants with good alignment #####\n')
    filtered_bfile = bfile + '_filtered'
    cmd = '%s --bfile %s ' % ( plink_path, bfile )
    cmd += '--extract %s ' % ( good_var_file )
    cmd += '--make-bed --out %s' % ( filtered_bfile )
    plink_output = subprocess.check_call( cmd, shell=True, stdout=log_fh, stderr=log_fh )
    log_fh.flush()

    return filtered_bfile


# check status of ref in alleles
# alleles could be G,A,T,C,.,+,-,D,I,0
# ref could be G,A,T,C,N
# param: string, path to plink bim file
# return: tuple of dict where keys are status of ref, values are variant names, dict where keys are variant names, values are ref alleles

def check_ref( bimfile ):
    log_fh.write( '\n##### checking ref alleles in %s #####\n' % ( bimfile ) )
    sys.stderr.write( 'checking ref alleles in %s\n' % ( bimfile ) )

    indel_alleles = set( [ 'D', 'I', '+', '-' ] )
    missing_alleles = set( [ '.', '0' ] )
    var_group = defaultdict( list )
    var2ref = {}
    fh = open( bimfile )

    for line in fh:
        line = line.rstrip( '\n' )
        if not line: continue

        chrom, var_name, cm, pos, a1, a2 = line.split( '\t' )
        ref = get_ref_allele( chrom, pos, var_name, args )
        ref = ref.upper()
        a1 = a1.upper()
        a2 = a2.upper()

        var2ref[var_name] = ref

        try:
            if ref == 'N':
                var_group['no_ref'].append( var_name )
            elif a1 in NTS and a2 in NTS:
                if a1 == ref or a2 == ref:
                    var_group['ref_present'].append( var_name )
                elif a1 == COMP[ref] or a2 == COMP[ref]:
                    var_group['flipped'].append( var_name )
                else:
                    var_group['impossible'].append( var_name )
            elif a1 in indel_alleles and a2 in indel_alleles:
                var_group['indels'].append( var_name )
            elif a1 in missing_alleles and a2 in missing_alleles:
                var_group['missing_variant'].append( var_name )
            else:
                var_group['illegal_alleles'].append( var_name )
        except Exception, e:
            print str( e )
            print line
            exit()


    fh.close()
    return ( var_group, var2ref )


# wrapper to plink command to flip variant alleles
# param: string, path to plink bfile stem
# param: string, path to list of variants to flip
# param: string, path to plink command
# param: log filehandle
# actions: calls plink --flip, creates output plink bedfile
# return: string, updated bfile stem
def flip( bfile, var_list_file, plink_path, log_fh ):
    sys.stderr.write( 'flipping vars in %s\n' % ( bfile ) )
    log_fh.write( '\n##### flipping vars in %s #####\n' % ( bfile ) )

    flipped_bfile = bfile + '_flipped'
    cmd = '%s --bfile %s ' % ( plink_path, bfile )
    cmd += '--flip %s ' % ( var_list_file )
    cmd += '--make-bed --out %s ' % ( flipped_bfile )

    subprocess.check_call( cmd, shell=True, stdout=log_fh, stderr=log_fh )
    log_fh.flush()

    return flipped_bfile


# wrapper to plink command to set a1 allele to ref allele
# param: string, path to plink bfile stem
# param: string, path to file with mapping of variant to ref allele
# param: string, path to plink command
# param: log filehandle
# actions: calls plink --a1-allele, creates output plink bedfile
# return: string, updated bfile stem
def set_real_ref( bfile, var2ref_file, plink_path, log_fh ):
    sys.stderr.write( 'setting a1 to ref allele in %s\n' % ( bfile )  )
    log_fh.write( '\n#####setting a1 to ref allele in %s #####\n' % ( bfile )  )

    real_ref_bfile = bfile + '_real_ref'
    try:
        cmd = '%s --bfile %s ' % ( plink_path, bfile )
        cmd += '--a1-allele %s ' % ( var2ref_file )
        cmd += '--make-bed --out %s ' % ( real_ref_bfile )
        subprocess.check_call( cmd, shell=True, stdout=log_fh, stderr=log_fh )
    except subprocess.CalledProcessError, e:
        sys.stderr.write( str( e ) )
        log_fh.write( str( e ) )
        sys.exit()
    log_fh.flush()

    return real_ref_bfile


# write file mapping variant name to ref allele
# param: string, path to output file
# param: var2ref, dict mapping variant name to ref allele
# param: log filehandle
# return: nothing
def write_var2ref_file( var2ref_file, var2ref, log_fh ):
    log_fh.write( '\n##### writing var2ref file: %s #####\n' % ( var2ref_file ) )
    sys.stderr.write( 'writing var2ref file: %s\n' % ( var2ref_file ) )

    var2ref_fh = open( var2ref_file, 'w' )
    for var in var2ref:
        ref = var2ref[var]
        var2ref_fh.write( '\t'.join( [ var, ref ] ) + '\n' )
    var2ref_fh.close()


# write file with list of variants to flip
# param: string, path to output file
# param: list, list of variants
# param: log filehandle
# return: nothing
def write_flip_file( outfile, variants, log_fh ):
    log_fh.write( '\n##### writing flip file %s #####\n' % ( outfile ) )
    sys.stderr.write( 'writing flip file %s\n' % ( outfile ) )
    flip_fh = open( outfile, 'w' )
    flip_fh.write( '\n'.join( variants ) + '\n' )
    flip_fh.close()


# Extract variants to flip from strand file in the format of W Rayner
# param: string, path to strand file
# return: list, variants to flip
def get_flips_from_strand_file( strand_file ):
    flips = []
    fh = open( strand_file )
    for line in fh:
        line = line.rstrip( 'n' )
        if not line: continue
        fields = line.split( '\t' )
        var = fields[0]
        strand = fields[4]
        if strand == '-':
            flips.append( var )
    fh.close()
    return flips


# wrapper to plink command to exclude ambiguous (A/T, C/G) variants from plink bfile
# param: string, path to plink bfile
# param: string, path to file with list of ambiguous variants
# param: log filehandle
# return: string, updated bfile stem
def exclude_ambiguous( bfile, ambiguous_file, log_fh ):
    log_fh.write( '\n##### excluding ambiguous variants' )
    sys.stderr.write( 'excluding ambiguous variants\n' )
    
    no_ambiguous_bfile = bfile + '_ambiguous_excluded'
    cmd = '%s --bfile %s ' % ( plink_path, bfile )
    cmd += '--exclude %s --make-bed --out %s' % ( ambiguous_file, no_ambiguous_bfile )
    subprocess.check_call( cmd, shell=True, stdout=log_fh, stderr= log_fh )
    
    return no_ambiguous_bfile


################################################################################
################################################################################
##
## main
##
################################################################################
################################################################################
args = parse_args()
bfile = os.path.basename( args.bfile )
log_fh = open( bfile + '_flip.log', 'w' )
var2ref_file = ''


TBX = None
VAR2REF = None
if args.ref_genome_tabix_file:
    TBX = pysam.TabixFile( args.ref_genome_tabix_file )
elif args.var_to_ref_file:
    var2ref_file = args.var_to_ref_file
    VAR2REF = {}
    fh = open( args.var_to_ref_file )
    for line in fh:
        line = line.rstrip( '\n' )
        if line.startswith( '#' ): continue
        if not line: continue

        var, ref = line.split()
        VAR2REF[var] = ref
    fh.close()
else:
    sys.stderr.write( 'need to specify --ref_genome_tabix_file or var_to_ref_file\n' )
    log_fh.write( 'need to specify --ref_genome_tabix_file or var_to_ref_file\n' )
    sys.exit()

NTS = set( [ 'G', 'A', 'T', 'C' ] )
COMP = { 'A':'T', 'C':'G', 'G':'C', 'T':'A' }

plink_path = args.plink_path or 'plink'

sys.stderr.write( 'checking for ambiguous alleles\n' )
log_fh.write( 'checking for ambiguous alleles\n' )
bimfile = args.bfile + '.bim'
ambiguous_vars, non_ambiguous_vars, allele2count = check_alleles( bimfile )

if ambiguous_vars:
    sys.stderr.write( 'ambiguous alleles present in %s\n' % ( args.bfile + '.bim' ) )
    log_fh.write( 'ambiguous alleles present in %s\n' % ( args.bfile + '.bim' ) )
    log_fh.flush()

ambiguous_file = None
if args.exclude_ambiguous:
    ambiguous_file = bfile + '_ambiguous_variants.txt'
    ambiguous_fh = open( ambiguous_file, 'w' )
    ambiguous_fh.write( '\n'.join( list( ambiguous_vars ) )  + '\n' )
    ambiguous_fh.close()
    

sys.stderr.write( 'making chrom, pos update files\n' )
chr_file, pos_file = make_chr_pos_update_files( args )

updated_chrom_bfile = update_chrom( bfile, plink_path, log_fh )
updated_pos_bfile = update_pos( updated_chrom_bfile, pos_file, plink_path, log_fh )
extracted_bfile = extract_aligned( updated_pos_bfile, pos_file, plink_path, log_fh )

# are alleles aligned to plus strand, minus strand, are missing, illegal, or impossible
# write out var2ref file for assigning ref alleles later
# var2ref_status = { var_name : { 'ref': ref, 
#                                 'status': ref|comp|impossible|illegal|missing
# impossible_allele: eg ref is G , allele is A/T

var_group, var2ref = check_ref( extracted_bfile + '.bim' )

if not var2ref_file:
    var2ref_file = extracted_bfile + '_var2ref.tsv'
    write_var2ref_file( var2ref_file, var2ref, log_fh )

flipped_bfile = None
if args.flip_non_ambiguous:
    if non_ambiguous_vars:
        vars_to_flip = var_group['flipped']
        flips_outfile = bfile + '_non-ambiguous_alleles_to_flip'
        write_flip_file( flips_outfile, vars_to_flip, log_fh )
        flipped_bfile = flip( extracted_bfile, flips_outfile, plink_path, log_fh )
    else:
        sys.stderr.write( "no non_ambiguous vars to flip\n" )
        log_fh.write( "\n##### no non_ambiguous vars to flip #####\n" )

if args.flip_from_strand_file:
    strand_flips = get_flips_from_strand_file( args.strand_file )
    flips_outfile = args.strand_file + '_vars_to_flip.txt'
    write_flip_file( flips_outfile, strand_flips, log_fh )
    flipped_bfile = flip( extracted_bfile, flips_outfile, log_fh )


if flipped_bfile:
    real_ref_bfile = set_real_ref( flipped_bfile, var2ref_file, plink_path, log_fh )
else:
    real_ref_bfile = set_real_ref( extracted_bfile, var2ref_file, plink_path, log_fh )

if args.exclude_ambiguous:
    non_ambiguous_bfile = exclude_ambiguous( real_ref_bfile, ambiguous_file, log_fh )


