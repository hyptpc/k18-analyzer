#!/opt/python-3.5/bin/python

#____________________________________________________

__author__  = 'Y.Nakada <nakada@km.phys.sci.osaka-u.ac.jp>'
__version__ = '2.0'
__date__    = '15 April 2018'

#____________________________________________________

import os
import sys
import time

sys.path.append( os.path.dirname( os.path.abspath( sys.argv[0] ) )
                 + '/module' )

import utility
from utility import pycolor as cl

#____________________________________________________

SLEEP_TIME = 5

#____________________________________________________

def display( filename ) :

    with open( fstatlist, 'r' ) as f :
        lines = f.readlines()

    data = list()
    for line in lines :
        words = line.split( '\t' )
        tmp = { 'key'     : cl.bold + words[0][:8] + cl.end,
                'nevents' : words[1],
                'unit'    : words[2],
                'binary'  : words[3][:16],
                'conf'    : words[4][:16],
                'data'    : words[5][:16],
                'root'    : words[6][:16],
                'stat'    : cl.reverce + cl.bold + cl.red + words[7] + cl.end,
                'time'    : words[8][:-1] }
        data.append( tmp )

    os.system( 'clear' )

    buff =  cl.reverce + cl.bold\
            + 'KEY'.ljust(8)   + '  '\
            + 'STAT'.ljust(14) + '  '\
            + 'BIN'.ljust(16)  + '  '\
            + 'CONF'.ljust(16) + '  '\
            + 'DATA(#EVENT)'.ljust(26) + '  '\
            + 'ROOT'.ljust(16) + '  '\
            + 'TIME'.ljust(8)\
            + cl.end
    print( buff )

    for item in data :
        buff = item['key'].ljust(8+8) + '  '\
               + item['stat'].ljust(14+18) + '  '\
               + item['binary'].ljust(16)  + '  '\
               + item['conf'].ljust(16)    + '  '\
               + item['data'].ljust(16)    + '(' + item['nevents'].rjust(8) + ')' + '  '\
               + item['root'].ljust(16)    + '  '\
               + item['time'].rjust(8)
        print( buff )

    buff = cl.reverce + cl.bold + 'Press \'Ctrl-C\' to exit' + cl.end
    print( buff )
    
#____________________________________________________

argvs = sys.argv
argc = len( argvs )

if argc != 2 :
    print( 'USAGE: %s [ file ]' % argvs[0] )
    sys.exit( 0 )

if not os.path.exists( argvs[1] ) :
    utility.ExitFailure( 'No such file: %s' % argvs[1] )

fstatlist = argvs[1]

#____________________________________________________


while True :

    try :
        display( fstatlist )
        time.sleep( SLEEP_TIME )

    except KeyboardInterrupt :

        print( 'KeyboardInterrupt' )
        print( 'Exiting the process...' )
        break
