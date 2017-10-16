#! /usr/bin/env python

import os
import sys

import argparse
import psycopg2
import select
import re
import time 


# check if illegal chars are present in string
# if yes, exit
# legal strings = a-z, A-Z, dot, underscore, 0-9
# param: string to sanitize
# param: string descriptor, eg tablename
def sanitize( s, d='' ):
    mo = re.search( r'([^\a-zA-Z0-9\._])', s )
    if mo:
        illegal_char = mo.group(1)
        sys.stderr.write( "illegal character %s in %s %s\n" % ( illegal_char, d, s ) )
        sys.exit()


# get db credentials from tab-separated file
# there should only be one set of credentials
# param: string, filepath
# return: tuple of username, password, db_port, host
def get_creds( creds_file ):
    fh = open( creds_file )
    for line in fh:
        line = line.rstrip( '\n' )
        line = line.lstrip()
        if not line: continue
        if line.startswith( '#' ): continue

        username, password, port, host = line.split( '\t' )
        return ( username, password, port, host )


# connect to postgres db
# param: string, db name
# param: string, path to tab-separated db credentials file
# return: dn connection object, db cursor object
def postgres_connect( db, creds_file ):
    user, pwd, port, host = get_creds( creds_file )
    dsn = 'host=%s port=%s dbname=%s user=%s password=%s' % ( host, port, db, user, pwd )
    conn = psycopg2.connect( dsn )
    curs = conn.cursor()

    return( conn, curs )



# process each line of annotation file
# param: string representing line of file
# param: argparse object
# return: tuple of variant name, annotation value
def process_line( line, args ):
    line = line.strip()
    if not line: return snp2value
    
    if args.header_starter:
        if line.startswith( args.header_starter ): return None

    fields = line.split( args.delim )

    snp_name = fields[args.snp_name_column]

    if args.db_column_datatype.lower() == 'boolean':
        if not args.present_snp_value:
            sys.stderr.write( "specify --present_snp_value\n" )
            sys.exit()
        else:
            value = args.present_snp_value
    else:
        # if value is missing, substitute default_value
        if len( fields ) < ( args.value_column + 1 ):
            value = args.default_value
        else:
            value = fields[args.value_column]
            
    return [ snp_name, value ]


# for summary table, get range of values in annotation
# param: argparse object
# return: string representing range
def get_range( args ):
    if args.db_column_datatype == 'boolean':
        f_range = 'True, False'
    elif args.db_column_datatype == 'varchar':
        f_range = 'NA'
    elif args.db_column_datatype == 'float' or args.db_column_datatype == 'integer':

        sql = 'SELECT MAX( %s ), MIN( %s ) FROM %s' % ( args.db_column_name,
                                                        args.db_column_name,
                                                        args.table )
        curs.execute( sql )
        col_max, col_min = curs.fetchone()
        f_range = '%s, $s' % ( col_min, col_max )
    else:
        msg = "Can't determine feature range for datatype %s\n" % ( args.db_column_datatype )
        raise ValueError( msg )


    return f_range


def parse_args():
    parser = argparse.ArgumentParser( description='''Add to or update snp annotations in a db table.''' )

    parser.add_argument( '--create_column', action='store_true', 
                         help='drop, create column; DELETES EXISTING DATA IN COLUMN' )
    parser.add_argument( '--create_index', action='store_true', help='add index on annotation column' )
    parser.add_argument( '--creds_file', 
                         help='''file with db credentials; tab-separated user, password, port, and host
                                 default: /home/cconnoll/chuckworking/annotation_db/gecco_db_creds''',
                         default='/home/cconnoll/chuckworking/annotation_db/gecco_db_creds' )
    parser.add_argument( '--default_value', default=None, help='default: None; used for missing values' )
    parser.add_argument( '--delim', default='\t', help='field delimiter, defaults to tab' )
    parser.add_argument( '--db', help="name of db; defaults to 'functional_annotation'",
                         default='functional_annotation' )
    parser.add_argument( '--db_column_name', required=True, help='name of feature' )
    parser.add_argument( '--db_column_datatype', required=True,
                         help='valid SQL datatype, needed if creating column here' )
    parser.add_argument( '--header_starter', help='start of header line if present' )
    parser.add_argument( '-p', '--present_snp_value', help='value to insert if snp is present in list' ) 
    parser.add_argument( '--snp_files', nargs='*',
                         help='file with snp names and optionally value' )
    parser.add_argument( '--snp_name_column', type=int, default=1,
                         help='column in snp_file with snp_name; defaults to 1')
    parser.add_argument( '--summary_table', 
                         help='table with summary stats; default:variant_annotation_feature_summary',
                         default='variant_annotation_feature_summary' )
    parser.add_argument( '--table', help="name of db table; defaults to 'variant_annotation'",
                         default='variant_annotation' )
    parser.add_argument( '--update_column', action='store_true', 
                         help="update column, requires that column exists already" )
    parser.add_argument( '--value_column', type=int, default=2,
                         help='column in snp_file with value; defaults to 2' )
    
    args = parser.parse_args()

    for arg in [ args.db, args.db_column_name, args.db_column_datatype, args.table ]:
        sanitize( arg )


    args.snp_name_column -= 1
    args.value_column -= 1
    return args



################################################################################
################################################################################
##
## main
##
################################################################################
################################################################################

# in some fashion, obtain list of snps to be added and associated feature value
# put these snps into a dict where snp_name is key, feature value is value
# loop over the summary files, appending the appropriate
# feature value for each snp
# need to append a header field as well
# a two column input (file or STDIN) is a key-value pair
# a one column input requires a feature value to be 
# appended for the snps in the list, and a value for
# absent snps


if __name__ == '__main__':
    args = parse_args()

    conn, curs = postgres_connect( args.db, args.creds_file )

    if args.create_column:
        drop_column_sql = '''ALTER TABLE {0} DROP COLUMN IF EXISTS {1}'''.format( args.table, args.db_column_name )
        curs.execute( drop_column_sql )
        conn.commit()

        add_column_sql = '''ALTER TABLE {0} ADD COLUMN {1} {2}'''.format( args.table, 
                                                                          args.db_column_name, 
                                                                          args.db_column_datatype )
        curs.execute( add_column_sql )
        conn.commit()
    elif args.update_column:
        sql = """SELECT column_name
                 FROM information_schema.columns 
                 WHERE table_name='{0}' and column_name='{1}'""".format( args.table, args.db_column_name )
        curs.execute( sql )
        row = curs.fetchone()
        if not row:
            sys.stderr.write( "column %s does not exist in table %s\n" % ( args.db_column_name, 
                                                                           args.table ) )
            sys.exit()
    else:
        sys.stderr.write( "please specify '--update_column' or '--create_column'\n" )
        sys.exit()

    update_sql = 'UPDATE {0} SET {1} = %s WHERE variant_name = %s'.format( args.table, 
                                                                           args.db_column_name )

    # read snps from stdin
    while sys.stdin in select.select([sys.stdin,],[],[],0.0)[0]:
        line = sys.stdin.readline()
        if line.startswith( '#' ): continue
        if line:
            snp_name, value = process_line( line, args )
            curs.execute( update_sql, [ value, snp_name] )
        
    # read snps from files
    values = []
    for filename in args.snp_files:

        # might be nice to fork and sbatch here, committing at the end, 
        # would need to set up intermachine communication. db? ipc over tcp/ip?
        # work ok as is, if resources become strained consider implementation

        sys.stderr.write( 'processing %s\n' % ( os.path.basename( filename ) ) )
        sys.stderr.flush()
        fh = open( filename )
        for line in fh:
            if args.header_starter:
                if line.startswith( args.header_starter ): continue
            snp_name, value = process_line( line, args )
            curs.execute( update_sql, [ value, snp_name ] )
            values.append( value )
        fh.close()
    sys.stderr.write( "committing\n" )
    conn.commit()

    if args.create_index:
        sys.stderr.write( "creating index\n" )
        index_sql = 'create index on {0} ( {1} )'.format( args.table, args.db_column_name ) 
        curs.execute( index_sql )
        conn.commit()

    # add to summary
    summary_sql = '''INSERT INTO {0} (feature_name, feature_datatype, feature_count, feature_range )
                     VALUES ( %s, %s, %s, %s )
                     ON CONFLICT DO UPDATE 
                       SET feature_name = %s, feature_datatype = %s, 
                           feature_count = %s, feature_range = %s'''.format( args.summary_table )

    if 'char' in args.db_column_datatype:
        args.db_column_datatype = 'string'

    f_count = len( values )
    f_range = get_range( args )

    curs.execute( summary_sql, [ args.db_column_name, args.db_column_datatype, f_count, f_range ] )
    conn.commit()
                                
