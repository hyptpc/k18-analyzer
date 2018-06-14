#!/usr/bin/env python3.5

#____________________________________________________

__author__  = 'Y.Nakada <nakada@km.phys.sci.osaka-u.ac.jp>'
__version__ = '2.1'
__date__    = '12 May 2018'

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

def decodeTime( second ) :

    second = int( second )

    hour    = second // 3600
    second -= hour    * 3600
    minute  = second // 60
    second -= minute  * 60

    return hour, minute, second

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
