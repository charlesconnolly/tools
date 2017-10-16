#! /usr/bin/env python

import os
import sys

import argparse
import math

from collections import defaultdict



# calculate percent identity as blat does
# ported from UCSC pslScore.pl
# param: int, query start
# param: int, query end
# param: int, number query gaps
# param: int, target start
# param: int, target end
# param: int, number target gaps
# param: int, number matches
# param: int, number mismatch
# param: int, number repeat matches
def get_pct_id( q_start, q_end, q_gaps, t_start, t_end, 
                t_gaps, n_match, n_mismatch, n_rep_match ):
    
    q_aln_size = q_end - q_start
    t_aln_size = t_end - t_start
    aln_size = min( [ q_aln_size, t_aln_size ] )

    if aln_size <= 0:
        return 0
    
    # ignore introns
    size_diff = q_aln_size - t_aln_size
    if size_diff < 0: size_diff = 0

    n_gaps = q_gaps + t_gaps

    size = n_match + n_rep_match + n_mismatch
    
    n_total = n_match + n_rep_match + n_mismatch
    if n_total:
        round_away_from_zero = 3*math.log( 1 + size_diff )

        if round_away_from_zero < 0:
            round_away_from_zero = int( round_away_from_zero - 0.5 )
        else:
            round_away_from_zero = int( round_away_from_zero + 0.5 )

        milli_bad = (1000*(n_mismatch + n_gaps + round_away_from_zero ))  / n_total
        
    pct_identity = 100 - ( milli_bad * 0.1 )
    
    return pct_identity


# sort by chromosome
# param: string, chromosome
# return: integer, chromosome
def by_chrom( q_name ):
    chrom2int = { 'X':23, 'Y':24, 'XY':25,'MT':26,'M':26 }

    fields = q_name.split( ';' )
    chrom = fields[2]
    chrom = chrom.replace( 'chr', '' )

    try:
        chrom = int( chrom )
    except ValueError:
        chrom = chrom2int[chrom]

    return chrom
    



parser = argparse.ArgumentParser()
parser.add_argument( '-o', '--outfile_stem', help='output file name stem, default: input file stem', )
parser.add_argument( '-p', '--print_strandfile', action='store_true', 
                     help='send strand file output to STDOUT' )
parser.add_argument( 'blat_file', help='blat .psl output file' )
args = parser.parse_args()


name2lines = defaultdict( list )
name2alns = {}
start = False
fh = open( args.blat_file )
for line in fh:
    line = line.rstrip( '\n' )
    if not line: continue

    if line.startswith( '----------' ):
        start = True
    elif start:
        fields = line.split()

        n_match      = int( fields[0] )
        n_mismatch   = int( fields[1] )
        n_rep_match  = int( fields[2] )
        q_gap        = int( fields[4] )
        
        q_size       = int( fields[10] )
        q_start      = int( fields[11] )
        q_end        = int( fields[12] )
        q_name       = fields[9]
        strand       = fields[8]
        
        t_gap   = int( fields[6] )
        t_chrom = fields[13]
        t_start = int( fields[15] )
        t_end   = int( fields[16] )
                 
        q_name_fields = q_name.split( ';' )
        name = q_name_fields[0]
#        name, snp, chrom, pos = q_name.split( ';' )

        if int( q_end ) == int( q_size ):
            aln3_prime = True
        else:
            aln3_prime = False


        # score defn from UCSC pslScore-1.pl
        psl_score = n_match + ( n_rep_match / 2 ) - n_mismatch - q_gap - t_gap

        # pct_identity from UCSC pslScore-1.pl
        pct_id = get_pct_id( q_start, q_end, q_gap, t_start, t_end, t_gap,
                             n_match, n_mismatch, n_rep_match )


        # need to keep track of variants without matches
        if not name in name2alns:
            name2alns[name] = []

        # if no matches, note alignments in .missing file
        name2lines[name].append( line )
                             
        # require 90% match (as per W Rayner)
        if pct_id > 0.9:
            if strand == '+':
                pos = str( t_end + 1 )
            elif strand == '-':
                pos = str( t_start )
            else:
                sys.stderr.write( "Bad strand for %s: %s\n" % ( q_name, strand ) )
                sys.exit()
        
            name2alns[name].append( { 'q_name' : q_name,  'pct_id'     : pct_id, 
                                      'n_match': n_match, 'psl_score'  : psl_score, 
                                      't_chrom': t_chrom, 'pos'        : pos, 
                                      'strand' : strand,  'aln3_prime' : aln3_prime,
                                      't_gap'  : t_gap,   'q_gap'      : q_gap,
                                      'q_start': q_start, 'q_end'      : q_end,
                                      'q_size' : q_size,  'n_mismatch' : n_mismatch
                                  } )

fh.close()


if not args.outfile_stem:
    args.outfile_stem = os.path.splitext( os.path.basename( args.blat_file ) )[0]

if args.print_strandfile:
    out_fh = sys.stdout
else:
    out_fh = open( args.outfile_stem + '.strand', 'w' )
    missing_fh = open( args.outfile_stem + '.missing', 'w' )
    missing_fh.write( '# entries without >90% match on illumina specified chromosome\n' )
    
    multi_fh   = open( args.outfile_stem + '.multi',   'w' )
    multi_fh.write( '# entries with more than one >90% match on illumina specified chromosome\n' )
    
    disc_fh    = open( args.outfile_stem + '.disc_strand',   'w' )
    disc_fh.write( '# entries with more than one good match on the illumina specified chromosome with different strands\n' )

    bad_3prime_fh = open( args.outfile_stem + '.bad_3-prime_alignment', 'w' )
    bad_3prime_fh.write("# 3'-end not aligned as detected by query_length != query_end\n")


out_fh.write( '\t'.join( [ '#name', 'chrom', 'position', 'pct_match', 'strand', 'illumina_alleles' ] ) + "\n" )

for name in sorted( name2alns ):
    n_matches = len( name2alns[name] )
    best_match = None

    # if no matches
    if n_matches == 0:
        missing.fh.write( '\n'.join( name2lines[name] ) + '\n' )

    elif n_matches == 1:
        best_match = name2alns[name][0]

    else:
        full_length_matches = []
        good_matches = []
        for match in name2alns[name]:

            # if 3' end does not align, exclude match, noting in .bad_aln file
            if not match['aln3_prime']:

                q_name = match['q_name']
                chrom  = match['t_chrom']
                pos    = match['pos']
                bad_3prime_fh.write( '\t'.join( [ q_name, chrom, pos ] ) + "\n" )
                
            # is it a full length match?
            elif ( match['n_mismatch'] == 0 and 
                   match['q_gap'] == 0      and 
                   match['t_gap'] == 0      and
                   (match['q_end'] - match['q_start'] == match['q_size'] )
                   ):

                full_length_matches.append( match )
                best_match = match

            elif best_match is None:
                best_match = match
                good_matches.append( match )
            elif best_match['psl_score'] < match['psl_score']:
                best_match = match
                good_matches.append( match )

        if len( full_length_matches ) > 1:
            multi_fh.write( '%s: multiple full length matches\n' % ( name ) )
            best_match = None
        elif len( good_matches ) > 1:    
            multi_fh.write( '%s: multiple pct_id > 90%% matches\n' % ( name ) )

    if best_match is not None:
        q_name = best_match['q_name']
        name_fields = q_name.split( ';' )
        ilmn_alleles = name_fields[1]
        
        chrom = best_match['t_chrom']
        pos   = best_match['pos']
        strand = best_match['strand']
        pct_id = best_match['pct_id']

        out_fh.write( '\t'.join( [ name, chrom, str( pos ), str( pct_id ), strand, ilmn_alleles ] ) + '\n' )





