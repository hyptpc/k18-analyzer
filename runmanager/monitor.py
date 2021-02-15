#!/usr/bin/env python3.5

#____________________________________________________

__author__  = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.0'
__date__    = '2 April 2019'

#____________________________________________________

import os
import sys
import time
import fcntl

import json

sys.path.append( os.path.dirname( os.path.abspath( sys.argv[0] ) )
                 + '/module' )

import RunManager
import utility
from utility import pycolor as cl

#____________________________________________________

DISPLAY_PERIOD = 5  # unit: second

#____________________________________________________

def display( filename ) :

    buff = str()
    with open( filename, 'r' ) as f :
        try :
            fcntl.flock( f.fileno(), fcntl.LOCK_SH )
        except IOError :
            # sys.stderr.write( 'ERROR: I/O error was detected.\n' )
            return
        else :
            buff = f.read()
        finally :
            fcntl.flock( f.fileno(), fcntl.LOCK_UN )

    try :
        info = json.loads( buff )
    except ValueError :
        return

    os.system( 'clear' )

    buff =  cl.reverce + cl.bold\
            + 'KEY'.ljust(8)   + '  '\
            + 'STAT'.ljust(16) + '  '\
            + 'BIN'.ljust(16)  + '  '\
            + 'CONF'.ljust(16) + '  '\
            + 'DATA(#EVENT)'.ljust(24) + '  '\
            + 'ROOT'.ljust(16) + '  '\
            + 'TIME'.ljust(8)\
            + cl.end
    print( buff )


    for key, item in info.items() :

        status = RunManager.SingleRunManager.decodeStatus( item )
        ptime  = RunManager.SingleRunManager.decodeTime( item )

        buff = cl.bold + key[:8].ljust(8) + cl.end + '  ' \
               + '{}{}{}{}{}'.format( cl.reverce, cl.bold, cl.red, status, cl.end ).ljust(16 + 18) + '  ' \
               + os.path.basename( item['bin'] )[-16:].ljust(16) + '  ' \
               + os.path.basename( item['conf'] )[-16:].ljust(16) + '  ' \
               + '{} ({})'.format( os.path.basename( item['data'] ), str( item['nev'] ) )[-24:].ljust(24) + '  ' \
               + os.path.basename( item['root'] )[-16:].ljust(16) + '  '\
               + ptime.rjust(8)
        print( buff )

    buff = cl.reverce + cl.bold + 'Press \'Ctrl-C\' to exit' + cl.end
    print( buff )
    
#____________________________________________________

def main( path ) :

    ptime = time.time()
    while True :
    
        try :
            display( path )
            dtime = DISPLAY_PERIOD - ( time.time() - ptime )
            if dtime > 0 :
                time.sleep( dtime )
            ptime = time.time()
    
        except KeyboardInterrupt :
    
            print( 'KeyboardInterrupt' )
            print( 'Exiting the process...' )
            break
    
        except FileNotFoundError :
    
            sys.stderr.write( 'Cannot find file > ' + fJobInfo + '\n' )
            sys.exit( 1 )

#____________________________________________________

if __name__ == "__main__" :

    argvs = sys.argv
    argc = len( argvs )
    
    if argc != 2 :
        print( 'USAGE: %s [ file ]' % argvs[0] )
        sys.exit( 0 )
    
    if not os.path.exists( argvs[1] ) :
        utility.ExitFailure( 'No such file > ' + argvs[1] )
    
    main( argvs[1] )

