#!/opt/python-3.5/bin/python

#____________________________________________________

__author__  = 'Y.Nakada <nakada@km.phys.sci.osaka-u.ac.jp>'
__version__ = '2.0'
__date__    = '15 April 2018'

#____________________________________________________

import sys

#____________________________________________________

def ExitFailure( message ):
    sys.stderr.write( 'ERROR: ' + message + '\n' )
    sys.exit( 1 )
    return None


#____________________________________________________

def updateJobStat( joblist, path ) :

    f = open( path, 'w' )
    for job, stat, time in joblist :
        info = map( str, job.getInfo() )
        for item in info :
            f.write( item + '\t' )
        f.write( stat + '\t' )
        f.write( '{:}:{:02d}:{:02d}'.format( *time ) + '\n' )
    f.close()

#____________________________________________________

class pycolor :

    black  = '\033[30m'
    red    = '\033[31m'
    green  = '\033[32m'
    yellow = '\033[33m'
    blue   = '\033[34m'
    purple = '\033[35m'
    cyan   = '\033[36m'
    white  = '\033[37m'

    end = '\033[0m'

    bold = '\033[1m'
    underline = '\033[4m'
    invisible = '\033[08m'
    reverce   = '\033[07m'
