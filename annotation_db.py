#! /app/easybuild/software/Python/2.7.13-foss-2016b/bin/python

import os
import sys

# or be sure psycopg2 can be loaded
sys.path.append( '/app/easybuild/software/Python/2.7.13-foss-2016b/lib/python2.7/site-packages' )

import psycopg2
import argparse
import re
import select
import textwrap
import traceback


# An iterator that uses fetchmany to keep memory usage down
# param: cursor
# param: (optional) ineger array size, default= 1000
# yields: 1 result
# throws: nothing
def ResultIter(cursor, arraysize=1000):
    while True:
        results = cursor.fetchmany(arraysize)
        if not results:
            break
        for result in results:
            yield result




# param: none
# returns: list of feature names
# throws: nothing
def get_feature_names():
    sql = '''SELECT feature_name FROM variant_annotation_feature_summary'''
    CURS.execute( sql )
    features = []
    for row in CURS.fetchall():
        features.append( row[0] )

    return features

def get_column_names():
    sql = '''SELECT column_name FROM information_schema.columns
             WHERE table_schema = %s AND table_name = %s'''
    CURS.execute( sql, [ SCHEMA, VAR_TABLE ] )

    column_names = []
    for row in ResultIter( CURS ):
        column_names.append( row[0] )

    return column_names

    



# param: list of features
# returns: dictionary where keys are features, values are counts
# throws: ValueError if param is not a list
def get_feature_count( features ):

    if not isinstance( features, list ):
        raise ValueError( "param is not a list" )

    validate_feature_names( features )

    # sql_tmpl = 'SELECT COUNT( %s ) FROM {0}.{1}'.format( SCHEMA, TABLE )

    sql = '''SELECT feature_count 
             FROM variant_annotation_feature_summary
             WHERE feature_name = %s'''
    feature2count = {}
    for feature in features:
        CURS.execute( sql, [ feature ] )
        count = CURS.fetchone()[0]
        feature2count[feature] = count

    return feature2count


# param: list of features
# returns: dictionary where keys are features, values are strings representing value range
# throws: ValueError if param is not a list
def get_feature_range( features ):
    if not isinstance( features, list ):
        raise ValueError( "param is not a list" )

    validate_feature_names( features )

    sql = '''SELECT feature_range 
             FROM variant_annotation_feature_summary 
             WHERE feature_name = %s'''

    feature2range = {}
    for feature in features:
        CURS.execute( sql, [ feature ] )
        feature2range[feature] = CURS.fetchone()[0]

    return feature2range

# param: list of features
# returns: dictionary where keys are features, values are datatypes
# throws: ValueError if param is not a list
def get_feature2datatype( features ):
    if not isinstance( features, list ):
        raise ValueError( "param is not a list" )

    validate_feature_names( features )


    sql = '''SELECT data_type FROM information_schema.columns
             WHERE table_schema = 'public' AND table_name = 'variant_annotation' AND column_name = %s'''

    feature2datatype = {}
    for feature in features:
        # CURS.execute( sql, [ SCHEMA, TABLE, feature ] )
        CURS.execute( sql, [ feature ] )
        datatype = CURS.fetchone()[0]

        if datatype.startswith( 'char' ):
            datatype = str
        elif datatype == 'integer':
            datatype = int
        elif datatype == 'boolean':
            datatype = bool
        elif datatype.startswith( 'double' ):
            datatype = float
        else:
            raise ValueError(  "Unknown data type '%s' for feature '%s'" % ( datatype, feature ) )
        
        feature2datatype[feature] = datatype

    return feature2datatype


def get_feature_datatype( feature ):
    f2dtype = get_feature2datatype( [ feature ] )
    return f2dtype[feature]



def get_dd():
    return dict


def get_rsids():
    # incorporate sqlite3 db?
    # incorporate merge db for up-to-date ids?
    return list

def get_variant_names():
    return list


def validate_feature_names( features ):

    if not isinstance( features, list ):
        raise ValueError( "parameter is not a list" )

    real_features = get_feature_names()
    for feature in features:
        if not feature in real_features:
            msg = "feature not found: %s" % ( feature )
            raise ValueError( msg )


'''
main actions:
- get_features
- get_feature_count, requires --feature
- get_feature_range, requires --feature
- get_feature_datatype, requires --feature
- get_variant_names, requires --rsid_list
- get_rsids, requires --snp_name_list
- get_variants, requires --where
'''


desc = '''This program can be used to 
  - describe the annotations in the database
  - annotate a list of variants
  - filter the variants in the database on annotation criteria
'''

def parse_args():
    parser = argparse.ArgumentParser( formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description=textwrap.dedent( '''\
Access the annotation database. Here, "feature" is 
synonymous with "annotation name". 
Output is written to STDOUT                                      
'''))

    ## db args
    parser.add_argument( '--creds_file', 
                         help='''tab-separated credentials file with username, password, 
                                 port number, and db host, in that order. default:
                                 /fh/fast/peters_u/GECCO_Working/chuckworking/annotation_db/gecco_db_creds''',
                         default='/home/cconnoll/chuckworking/annotation_db/gecco_db_creds' )
    
    parser.add_argument( '--debug', action='store_true' )
    # metadata
    parser.add_argument( '--describe_features', nargs='*',
                         help='''print feature name, datatype, number of variants, and 
                                 range for specified features or all features if none are specified''' )

    # annotate a list of variants
    parser.add_argument( '-a', '--annotate', action='store_true', 
                         help='annotate a list of variants, requires --feature' )
    parser.add_argument( '--features', help='valid annotation name', nargs='+' )

    # get list of variants based on feature criteria
    parser.add_argument( '--filter', action='store_true', 
                         help='get variants by specified criteria. Requires "--where".' )

    # filter criteria
    parser.add_argument( '--where', action='append', help='''"feature op val"; needs to be enclosed in quotes. 
                                                        feature is valid feature name;
                                                        op is either "gt", "lt", or "eq" meaning greater 
                                                        than, less than, or equals, respectively; 
                                                        value is valid datatype for feature; 
                                                        example: "cadd gt 5", "di_are eq True"''' )

    parser.add_argument( '--variants_file', 
                         help='''file with list of variants, 1 per line, as variant_name or  rsID; 
                                 use "stdin" to indicate reading from stdin''' )
    
    parser.add_argument( '-o', '--outfile', help='file to which to write output' )
    args = parser.parse_args()

    if args.filter and not args.where:
        raise ValueError( """specify variant selection criteria with '--where'.\n""" )

    if args.annotate and not args.variants_file:
        raise ValueError( """specify variant names with '--variants_file'.\n""" )
                             

    return args


def is_unclean( s ):
    if re.search( r'[^a-zA-Z0-9_.]', s ):
        return True
    else:
        return False

def get_creds( creds_file ):
    fh = open( creds_file )
    for line in fh:
        line = line.rstrip( '\n' )
        line = line.lstrip()
        if not line: continue
        if line.startswith( '#' ): continue

        # there should only be one set of credentials
        username, password, port, host = line.split( '\t' )
        return ( username, password, port, host )


def postgres_connect( db, creds_file ):
    user, pwd, port, host = get_creds( creds_file )
    dsn = 'host=%s port=%s dbname=%s user=%s password=%s' % ( host, port, db, user, pwd )
    conn = psycopg2.connect( dsn )
    curs = conn.cursor()

    return( conn, curs )


# param: string of the form 'feature op value'
# returns: True
# throws: ValueError if string is incorrectly formatted
def validate_where( where ):
    where = where.strip()
    where = where.lower()


    fields = []
    if 'or' in where:
        try:
            fields = where.split( None, 2 )
        except ValueError:
            pass
    else:
        fields  = where.split()

    if not len( fields ) == 3:
        msg = "invalid format of 'where' constraint\nshould be 'feature op value [or value or value]'\n"
        msg += "where feature is valid feature name, op is 'lt', 'gt', or 'eq'\n"
        msg += "and value is a correct datatype for the feature eg float, boolean\n"
        msg += "Example --where 'cadd gt 5.0' or --where 'coding eg splicing or stopgain'\n"
        raise ValueError( msg )

    feature, op, val = fields
    features = get_feature_names()
    datatype = get_feature_datatype( feature )

    if not feature in features:
        msg = "feature '%s' is not a valid feature name\n" % ( feature )
        raise ValueError( msg )

    if op != 'gt' and op != 'lt' and op != 'eq':
        msg = "invalid op '%s'\nvalid operators are:\n" % ( op )
        msg += "'lt' (less than), 'gt' (greater than), 'eq' (equals)" 
        raise ValueError( msg )

    try:
        datatype( val )
    except ValueError:
        msg = "invalid datatype for '%s': '%s'\n" % ( feature, val )
        msg += "datatype for '%s' is '%s'" % ( feature, datatype.__name__ ) 
        raise ValueError( msg )

    return True

def variant_list( filename ):
    variants = []
    if filename:
        fh = open( filename )
        for line in fh:
            line = line.rstrip( '\n' )
            if not line: continue
            if line.startswith( '#' ): continue
            
            variants.append( line )

    while sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
        line = sys.stdin.readline()
        line = line.rstrip( '\n' )
        if line:
            variants.append( line )

    return variants
            

# param: list of features
# returns: list of lists
# throws: ValueError if feature is not a column in the db table
def describe_features( features ):
    if not isinstance( features, list ):
        raise ValueError( "internal error: param is not a list" )

    sql = '''SELECT feature_name, feature_datatype, Feature_count, feature_range
             FROM variant_annotation_feature_summary 
             WHERE feature_name = %s'''

    rows = [ [ '# name', 'datatype', 'count', 'range' ] ]
    all_features = get_feature_names()
    for feature in features:
        if not feature in all_features:
            msg = "feature not found: %s" % ( feature )
            raise ValueError( msg )
        else:
            CURS.execute( sql, [ feature ] )
            rows.append( list( CURS.fetchone() ) )

    return rows



# param: where constraint in the form 'feature op value'
# param: (optional) boolean to indicate whether to include rsid in output
def get_variants( wheres ):
    if not isinstance( wheres, list ):
        raise ValueError( "internal error: param is not a list" )

    for where in wheres:
        validate_where( where )


    op2symbol = { 'lt':'<', 'eq':'=', 'gt':'>' }
    features = []
    vals = []
    all_clauses = []
    for where in wheres:
        where = where.lower()

        if 'or' in where:
            feature, op, values = where.split( None, 2 )
            values = [ s.strip() for s in values.split( 'or' ) ]
        else:
            feature, op, value = where.split()
            values = [value]

        validate_feature_names( [ feature ] )

        features.append( feature )
        op = op2symbol[op]
        vals += values

        clauses = []
        for value in values:
            clause = '%s %s %%s' % ( feature, op )
            clauses.append( clause )

        clauses = '( ' + ' OR '.join( clauses ) + ' )'
        all_clauses.append( clauses )



    feature_names = [ 'variant_name', 'rsid' ]
    feature_names += features 

    sql = 'SELECT {0} FROM variant_annotation '.format( ', '.join( feature_names ) )
    sql += 'WHERE %s ' % ( ' AND '.join( all_clauses ) )

    CURS.execute( sql, vals )
    rows = [ feature_names ]
    for row in ResultIter( CURS ):
        rows.append( row )
    
    return rows



# param: variant name in the form chrom:pos or chrom:pos_ref/alt
# return: tuple of chrom, pos
# throws: ValueError if variant name is badly formed
def validate_variant_name( variant_name ):
    if not ':' in variant_name:
        msg = "invalid variant name format: '%s'" % ( variant_name ) 
        raise ValueError( msg )
    fields = variant_name.split( '_' )
    chrom, pos = fields[0].split( ':' )
    chrom = chrom.replace( 'chr', '' )
    if chrom == 'X': 
        chrom = 23
    elif chrom == 'Y': 
        chrom = 24
    elif chrom == 'XY': 
        chrom = 26
    try:
        int( chrom )
        int( pos )
    except ValueError:
        msg = "invalid variant name format: '%s'" % ( variant_name ) 
        raise ValueError( msg )

    return ( chrom, pos )



def get_meta():
    sql = '''SELECT feature_name, feature_datatype, feature_count, feature_range 
             FROM variant_annotation_feature_summary'''
    CURS.execute( sql )

    feature2meta = {}
    for f_name, f_dtype, f_count, f_range in CURS.fetchall():
        feature2meta[f_name] = { 'count': f_count, 'range':f_range, 'dtype':f_dtype }

    return feature2meta


def annotate( variants_file, features ):
    if features is None:
        features = get_column_names()

    for f in [ 'rsid', 'variant_name', 'id' ]:
        try:
            features.remove( f )
        except ValueError, e:
            pass

    features = [ 'variant_name', 'rsid' ] + features

    sql_tmpl = '''SELECT {0}
                  FROM variant_annotation 
                  WHERE %s = %%s'''.format( ','.join( features ) )

    rows = [ features ]
    fh = open( variants_file )
    for line in fh:
        line = line.rstrip( '\n' )
        if not line: continue
        if line.startswith( '#' ): continue
        
        input_type = ''
        if line.startswith( 'rs' ):
            input_type = 'rsid'
        elif re.search( r'\d{1,2}:\d+', line ):
            input_type = 'variant_name'
        else:
            raise ValueError( 'unknown input type' )

        CURS.execute( sql_tmpl % ( input_type ), [ line ] )
        for row in ResultIter( CURS ):
            rows.append( list( row ) )
    fh.close()

    return rows
        
# formats None as 'NA', stringify other types
def string_format( s ):
    if s is None:
        return 'NA'
    return str( s )
                      


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#
# main
# 
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

try:
    ARGS = parse_args()
except ValueError, e:
    sys.stderr.write( str( e ) )
    sys.exit()
    
    

db = 'functional_annotation'
CONN, CURS = postgres_connect( db, ARGS.creds_file )
SCHEMA = 'public'
VAR_TABLE  = 'variant_annotation'
SUMMARY_TABLE = 'variant_annotation_feature_summary'

'''
main actions:
- get_features
- describe_feature
- get_variant_names, requires --rsid_list
- get_rsids, requires --snp_name_list
- filter_variants, requires --where
'''
# operator, how may I direct your call?
try:
    feature2meta = get_meta()

    if ARGS.features:
        for f in ARGS.features:
            if not f in feature2meta:
                raise ValueError( 'unknown feature: %s' % ( f ) )


    rows = []
    # these can stand alone
    if ARGS.describe_features is not None:

        if len( ARGS.describe_features ) > 0:
            names = ARGS.describe_features
        else:
            names = feature2meta.keys()

        rows.append( [ 'feature_name', 'feature_datatype', 'feature_count', 'feature_range' ] )
        for feature in sorted( names):            
            try:
                f_count = feature2meta[feature]['count']
                f_range = feature2meta[feature]['range']
                f_dtype = feature2meta[feature]['dtype']
            except KeyError:
                raise Exception( "Unknown feature '%s'\n" % ( feature ) )

            rows.append( [ feature, f_dtype, f_count, f_range ] )


    elif ARGS.annotate:
        rows = annotate( ARGS.variants_file, ARGS.features )

    elif ARGS.filter:
        rows = get_variants( ARGS.where )

    if ARGS.outfile:
        out_fh = open( ARGS.outfile, 'w' )
    else:
        out_fh = sys.stdout

    for row in rows:
        out_fh.write( '\t'.join( map( string_format, row ) ) + '\n' )

        

# error catcher
except Exception, e:
    if ARGS.debug:
        traceback.print_exc()
    sys.stderr.write( "Error: %s\n" % ( str( e ) ) )
    sys.exit( 1 )

