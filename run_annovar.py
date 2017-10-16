#! /usr/bin/env python

import os
import sys
import argparse

import subprocess




parser = argparse.ArgumentParser()
parser.add_argument( '--annovar_path', help='path to annovar', default='annovar' )
parser.add_argument( '--annovar_db_path', help='path to directory with annovar dbs', 
                     required=True )
parser.add_argument( '--db', help='name of annovar db', required=True,
                     choices=[ 'refGene', 
                               'ensGene', 
                               'ljb26_all',
                               'ccdsGene'
                           ],
                     action='append')
parser.add_argument( '-p', '--partition', default='restart', 
                     choices=[ 'restart', 'campus', 'peters_u', 'local' ] )
parser.add_argument( '--sbatch_arg', action='append' )
parser.add_argument( 'input_files', nargs='+', help='annovar-formatted input files' )
args = parser.parse_args()

# how do we know what type of operation each protocol will be?
# should we use a config file to map these?
db2op = { 'refGene' : 'g',
          'ensGene' : 'g',
          'ljb26_all' : 'f',
          'ccdsGene' : 'g'
          }
          
operations = []
for db in args.db:
    if not db in db2op:
        sys.stderr.write( "Unknown db: '%s'\n" % ( db ) )
        sys.exit()
    else:
        operations.append( db2op[db] )

for filename in args.input_files:
    # table_annovar.pl infile annovar_db -out outfile -remove -nastring . -buildver hg19 -protocol refGene,ljb26_all  -operation g,f --otherinfo
    
    outfile = os.path.splitext( os.path.basename( filename ) )[0]
    cmd = '%s %s %s -out %s -remove -nastring . -buildver hg19 ' % ( args.annovar_path, filename, args.annovar_db_path, outfile )
    cmd += '-protocol %s ' % ( ','.join( args.db ) )
    cmd += '-operation %s ' % ( ','.join( operations ) )
    cmd += '--otherinfo '

    

    if args.partition == 'local':
        subprocess.check_call( cmd, shell=True )
    else:
        sbatch_cmd = 'sbatch -p %s ' % ( args.partition )

        if args.sbatch_arg is not None:
            for sbatch_arg in args.sbatch_arg:
                sbatch_cmd += sbatch_arg + ' ' 

        sbatch_cmd += '--wrap "%s"' % ( cmd )
            
#        print sbatch_cmd
        subprocess.check_call( sbatch_cmd, shell=True )
    
