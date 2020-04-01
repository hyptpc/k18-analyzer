#!/usr/bin/env python3.5

#____________________________________________________

__author__  = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.0'
__date__    = '2 April 2019'

#____________________________________________________

import os
import sys
import time
import signal
import fcntl

import json

sys.path.append( os.path.dirname( os.path.abspath( sys.argv[0] ) )
                 + '/module' )

import utility
import RunlistManager
import RunManager


#____________________________________________________

def handler( num, frame ) :

    print( 'KeyboardInterrupt' )

    print( 'Terminating processes...' )

    runMan = RunManager.RunManager()
    runMan.kill()
    runMan.dumpStatus()

    time.sleep( 1 )             # waiting until bsub log files are generated

    fl = False
    tmp = input( 'Keep log files? [y/-] >> ' )
    if 'y' == tmp :
        print( 'Deleting intermediate files...' )
        fl = True
    else :
        print( 'Deleting files...' )

    runMan.finalize( fl )

    sys.exit( 0 )


#____________________________________________________

def main( runlist_path, status_path ) :

    runlistMan = RunlistManager.RunlistManager()
    runlistMan.setRunlistPath( runlist_path )

    runMan = RunManager.RunManager()
    runMan.setRunlistManager( runlistMan )
    runMan.setStatusOutputPath( status_path )
    runMan.initialize()

    runMan.dumpStatus()

    signal.signal( signal.SIGINT, handler )

    runMan.registerStaging()
    runMan.run()


#____________________________________________________

if __name__ == "__main__" :

    argvs = sys.argv
    argc = len( argvs )

    if argc != 2 :
        print( 'USAGE: {} [ file ]'.format( argvs[0] ) )
        sys.exit( 0 )

    if not os.path.exists( argvs[1] ) :
        utility.ExitFailure( 'No such file: {}'.format( argvs[1] ) )

    fRunlist = argvs[1]

    print( 'Press \'Ctrl-C\' to terminate processes.' )

    dJobInfo = os.path.dirname( os.path.abspath( sys.argv[0] ) ) + '/stat'
    if not os.path.exists( dJobInfo ) :
        os.mkdir( dJobInfo )

    fJobInfo = dJobInfo + '/'\
           + os.path.splitext( os.path.basename( fRunlist ) )[0]\
           + '.json'
    if not os.path.exists( fJobInfo ) :
        with open( fJobInfo, 'w' ) as f : pass

    main( fRunlist, fJobInfo )
