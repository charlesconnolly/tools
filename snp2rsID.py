#! /usr/bin/env python

import os
import sys
import sqlite3
import argparse
import select
import textwrap

from collections import defaultdict


# queries sqlite3 db to find rsID
# input: chromosome
# input: position
# output: rsID or 'not found'


################################################################################
#
# parse_args()
# param: none
# return: argparse Namespace object
#
################################################################################
def parse_args():
    parser = argparse.ArgumentParser( formatter_class=argparse.RawDescriptionHelpFormatter,     
                                      description=textwrap.dedent('''\
This program takes a list of snp_names, rsids or chromosome/position pairs, and prints
the original input, current rsid, chromosome, and start position, end position, and alleles.
The start position is 0-based, and the end position is 1-based, as in a bedfile.
When the query is rsID or snp_name, it can be taken from STDIN or a whitespace-separated file. 
When the query is chromosome and position, it must be read from a file. The relevant columns 
can be specified if different from the default.''' ))
                                      
                                      
    parser.add_argument( '-i', '--input_type', required=True,
                         choices=[ 'snp_name', 'chrom_pos', 'rsid' ] ) 
    parser.add_argument( '--snp_db', help='path to dbSNP sqlite3 db', required=True )

    parser.add_argument( '--merge_db', help='path to dbSNP sqlite3 merge db', required=True )

    parser.add_argument( '--snp_name_column', type=int, default=1,
                         help='index of snp_name column, defaults to 1' )
    parser.add_argument( '--rsid_column', type=int, default=1,
                         help='index of rsid column, defaults to 1' )
    parser.add_argument( '--chromosome_column', default=1, type=int,
                         help='index of chromosome column, defaults to 1' )
    parser.add_argument( '--position_column', default=2, type=int,
                         help='index of columns with position, defaults to 2' )
    
    parser.add_argument( '--get_class', action='store_true', 
                         help='when input is rsid, get snp class eg single, deletion' )
    parser.add_argument( 'input', nargs='*' )
    args = parser.parse_args()

    args.snp_name_column -= 1
    args.rsid_column -= 1
    args.chromosome_column -= 1
    args.position_column -= 1

    return args




################################################################################
#
# get_inputs()
# param: args
# return: list either chrom:pos, snp_name, or rsid
#
################################################################################
def get_inputs( args ):
    inputs = []

    if args.input:
        for idx in range( len( args.input ) ):
            elem = args.input[idx]
    
            if os.path.isfile( elem ):
                fh = open( elem )
                for line in fh:
                    line = line.strip()
                    if not line: continue
                    if line.startswith( '#' ): continue
                
                    if args.input_type == 'chrom_pos':
                        fields = line.split()
                        chrom = fields[args.chromosome_column]
                        pos   = fields[args.position_column]
                        chrom = chrom.replace( 'chr', '' )
                        inputs.append( '%s:%s' % ( chrom, pos ) )
                    elif args.input_type == 'snp_name':
                        fields = line.split()
                        snp_name = fields[args.snp_name_column]
                        inputs.append( snp_name )

                    elif args.input_type == 'rsid':
                        fields = line.split()
                        #line = line.replace( 'chr', '' )
                        rsid = fields[args.rsid_column]
                        inputs.append( rsid )

            else:
                if args.input_type == 'snp_name' or args.input_type == 'rsid':
                    inputs.append( elem )
                else:
                    sys.stderr.write( 'only rsIDs or snp_names can be stdin inputs' )
                    sys.exit()
    elif select.select([sys.stdin,],[],[],0.0)[0]:
        inputs = [ line.rstrip( '\n' ) for line in sys.stdin.readlines() ]

    if not inputs:
        raise ValueError( "No inputs" )

    return inputs    




################################################################################
#
# get_cursor()
# params: args
# returns: sqilte3 cursor object
#
################################################################################
def get_cursor( db ):
    if not os.path.isfile( db ):
        sys.stderr.write( "sqlite3 database '%s' not found\n" % ( db ) )
        sys.stderr.write( 'please specify a sqlite3 db\n' )
        sys.exit()

    conn = sqlite3.connect( db )
    curs = conn.cursor()

    return curs


def get_chrom_start_end( snp_name ):
    # snp_name could look like 
    #  chr1:1234
    #  chr1:1234_C/A
    #  chr1:1234_C/AAA
    #  chr1:1234_CCCC/A
    # need to find correctly formatted start and end
    # dbSNP uses 0-based start, 1-based end

    fields = snp_name.split( '_' )
    chrom, pos = fields[0].split( ':' )
    pos = int( pos )
    
    # chrom may look like 'chr1', or '1'
    # we need 'chr1'
    chrom = chrom.replace( 'chr', '' )
    chrom = 'chr%s' % ( chrom )
    
    # if entry is of the form chr1:1234
    if len( fields ) == 1:
        chromStart = int( pos ) - 1
        chromEnd = pos
    else:
        ref, alt = fields[1].split( '/' )

        # deletion eg  10:8025137_ACAGT/A
        if len( ref ) > 1 and not ref.isdigit():
            chromStart = pos 
            chromEnd   = pos + len( ref ) - 1

#        # SV eg  10:7932533_7942689/<DEL>
#        elif alt == '<DEL>':
#            end = ref
#            ref = '-'
#            alt = '-'
            
        # insertion eg 10:8012618_A/AC
        elif len( alt ) > 1 and not alt == '<DEL>':
            chromStart = pos
            chromEnd   = pos

        # SNP eg 11:128446_G/A
        else:
            chromStart = pos -1 
            chromEnd   = pos


    return ( chrom, chromStart, chromEnd )





################################################################################
#
# prep_input() 
# param: list of inputs
# param: string; input_type
# return: tuple; sql statement, list of list of sql params
#
################################################################################
def prep_input( input_elem, input_type ):
    sql = ''
    params = []
    if input_type == 'snp_name':
        sql = '''SELECT name, chrom, chromStart, chromEnd, alleles FROM dbsnp 
                 WHERE chrom = ? AND chromStart = ? AND chromEnd = ?'''

        # elem could be snp_name like 1:1234_C/T or 1:1234 or 1:1234_CCC/A or 1:1234_C/AAAA
        chrom, start, end = get_chrom_start_end( input_elem )
        params = [ chrom, start, end ] 
            
    elif input_type == 'chrom_pos':
        sql = '''SELECT name, chrom, chromStart, chromEnd, alleles FROM dbsnp 
                 WHERE chrom = ? and chromEnd = ?'''
        
        chrom, pos = input_elem.split( ':' )
        chrom = chrom.replace( 'chr', '' )
        chrom = 'chr%s' % ( chrom )
        params = [ chrom, pos ]

    elif input_type == 'rsid':
        sql = 'SELECT name, chrom, chromStart, chromEnd, alleles'
        if args.get_class:
            sql += ', class'
        sql += ' FROM dbsnp WHERE name = ?'
        params = [ input_elem ]

    else:
        sys.stderr.write( 'Bad input type: "%s"\n' % ( input_type ) )
        sys.exit()

    return ( sql, params )





def get_merge_output( rsid ):
    rsid = int( rsid.replace( 'rs', '' ) )

    rsid_sql = '''SELECT DISTINCT rsCurrent from rs_merge 
                  WHERE rsLow = ? or rsHigh = ?'''
    rsid_curs.execute( rsid_sql, [ rsid, rsid ] )
    rows = rsid_curs.fetchall()

    return rows

    

# returns name, chrom, chromEnd if successful
# returns None if unsuccessful
def get_snp_output( snp_curs, sql, param ):
    snp_curs.execute( sql, param )
    rows = snp_curs.fetchall()

    if rows:
        return rows
    else:
        return None




################################################################################
################################################################################
##
## main
##
################################################################################
################################################################################

args = parse_args()
try:
    inputs = get_inputs( args )
except ValueError, e:
    sys.stderr.write( "No inputs specified\n" )
    sys.exit()

snp_curs = get_cursor( args.snp_db )

if args.input_type == 'rsid':
    rsid_curs = get_cursor( args.merge_db )

found = []
not_found = []
for elem in inputs:
    sql, param = prep_input( elem, args.input_type )

    rows = get_snp_output( snp_curs, sql, param )
    if rows:
        for row in rows:
            found.append( [elem] + list( row ) )

    else:
        # sometimes get a retired rsID, check for new
        if args.input_type == 'rsid':
            rows = get_merge_output( elem )

            if rows:
                for row in rows:
                    rs_current = 'rs%d' % ( row[0] )
                    snp_rows = get_snp_output( snp_curs, sql, [rs_current] )

                    if snp_rows:
                        for snp_row in snp_rows:
                            found.append( [elem] + list( snp_row ) )
                    
                    else:
                        not_found.append( elem )
            else:
                not_found.append( elem )

        else:
            not_found.append( elem )


print '# db: %s' % ( args.snp_db )
if not_found:
    print '# these identifiers were not found in the db'
    for elem in not_found:
        print '# %s' % ( elem )
    print


if found:
    labels = [ '# input', 'current_rsID', 'chrom', 'start_position', 'end_position', 'alleles' ]
    if args.input_type == 'rsid' and args.get_class:
        labels.append( 'class' )

    print '\t'.join( labels )
                       
    for entry in found:
        print '\t'.join( map( str, entry ) )

        

