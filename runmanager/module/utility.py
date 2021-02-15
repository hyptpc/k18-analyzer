#!/usr/bin/env python3

#____________________________________________________

__author__  = 'Y.Nakada <nakada@km.phys.sci.osaka-u.ac.jp>'
__version__ = '4.0'
__date__    = '2 April 2019'

#____________________________________________________

import sys

#____________________________________________________

def ExitFailure( message ):
    sys.stderr.write( 'ERROR: ' + message + '\n' )
    sys.exit( 1 )
    return None


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
