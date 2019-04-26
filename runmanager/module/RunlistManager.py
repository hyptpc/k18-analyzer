#!/usr/bin/env python3.5

#____________________________________________________

__author__  = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.0'
__date__    = '2 April 2019'

#____________________________________________________

import os
import copy

import yaml

#____________________________________________________

import utility
import Singleton

#____________________________________________________

#____________________________________________________
#
# Runlist Manager Class
#____________________________________________________

class RunlistManager( metaclass = Singleton.Singleton ) :

    #____________________________________________________
    def __init__( self ) :

        self.__baseName = None
        self.__workdir  = None
        self.__keys     = list()
        self.__runlist  = list()
        self.__flReady  = None


    #____________________________________________________
    def setRunlistPath( self, path ) :

        self.__baseName = os.path.splitext( os.path.basename( path ) )[0]
        self.__workdir  = self.getWorkDir( path )
        self.__keys     = list()
        self.__runlist  = list()
        self.makeRunlist( path )


    #____________________________________________________
    def getTag( self ) :

        return self.__baseName


    #____________________________________________________
    def getNRuns( self ) :

        return len( self.__keys )


    #____________________________________________________
    def getKeys( self ) :

        return list( self.__keys )

    #____________________________________________________
    def getRunInfo( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return dict( self.__runlist[index] )
        else :
            return None


    #____________________________________________________
    def getKey( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return str( self.__runlist[index]['key'] )
        else :
            return None


    #____________________________________________________
    def getExecPath( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return str( self.__runlist[index]['bin'] )
        else :
            return None


    #____________________________________________________
    def getDataPath( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return str( self.__runlist[index]['data'] )
        else :
            return None


    #____________________________________________________
    def getConfPath( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return str( self.__runlist[index]['conf'] )
        else :
            return None


    #____________________________________________________
    def getOutPath( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return str( self.__runlist[index]['root'] )
        else :
            return None


    #____________________________________________________
    def getNProc( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return int( self.__runlist[index]['nproc'] )
        else :
            return None



    #____________________________________________________
    def getBuffPath( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return str( self.__runlist[index]['buff'] )
        else :
            return None


    #____________________________________________________
    def getQueue( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return str( self.__runlist[index]['queue'] )
        else :
            return None


    #____________________________________________________
    def getDivUnit( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return int( self.__runlist[index]['unit'] )
        else :
            return None


    #____________________________________________________
    def getNEvents( self, index ) :

        index = self.__getIndex( index )
        if not index is None :
            return int( self.__runlist[index]['nevents'] )
        else :
            return None


    #____________________________________________________
    def __getIndex( self, index ) :

        if index in self.__keys :
            return self.__keys.index( index )
        elif isinstance( index, int ) \
            and len( self.__runlist ) > index :
            return index
        else :
            return None


    #____________________________________________________
    def getWorkDir( self, path ) :

        if self.__flReady is False :
            return None

        data = dict()
        with open( path, 'r' ) as f :
            data = yaml.load( f.read() )

        tmp = os.path.expanduser( data['WORKDIR'] )
        workdir = tmp if os.path.exists( tmp ) and os.path.isdir( tmp ) \
                  else utility.ExitFailure( 'Cannot find directory: ' + tmp )

        return workdir


    #____________________________________________________
    def makeRunlist( self, path ) :

        if self.__flReady is False :
            return

        cdir = os.getcwd()
        raw_runlist = self.decodeRunlist( path )
        os.chdir( self.__workdir )

        runlist = list()
        for item in raw_runlist :

            run = dict()

            run['key'] = item[0]

            pbin = None
            if os.path.exists( item[1]['bin'] ) \
                and os.path.isfile( item[1]['bin'] ) :
                pbin = item[1]['bin']
            else :
                utility.ExitFailure( 'Cannot find file: ' + item[1]['bin'] )
            run['bin'] = pbin

            runno = None
            if os.path.exists( item[1]['data'] ) and os.path.isfile( item[1]['data'] ) :
                tmp = os.path.splitext( os.path.basename( item[1]['data'] ) )[0]
                runno = int( tmp[3:8] ) if tmp[3:8].isdigit() else None
            else :
                runno = item[0] if isinstance( item[0], int ) else None

            pdata = self.makeDataPath( item[1]['data'], runno )
            run['data']    = pdata
            run['nevents'] = self.getNEvents( os.path.dirname( os.path.abspath( pdata ) ), runno )

            pconf = None
            if  os.path.exists( item[1]['conf'] ) :
                if os.path.isfile( item[1]['conf'] ) :
                    pconf = item[1]['conf']
                elif os.path.isdir( item[1]['conf'] ) and not runno is None :
                    pconf = item[1]['conf'] + '/analyzer_%05d.conf' % runno
            if pconf is None :
                utility.ExitFailure( 'Cannot decide conf file path' )
            run['conf'] = pconf

            base = item[0] + os.path.basename( pbin ) if runno is None \
                    else 'run' + '%05d_' % runno + os.path.basename( pbin )
            run['root'] = self.makeRootPath( item[1]['root'], base )

            run['unit']  = item[1]['unit']  if isinstance( item[1]['unit'], int )  else 0
            run['queue'] = item[1]['queue'] if isinstance( item[1]['queue'], str ) else 's'
            run['nproc'] = item[1]['nproc'] if isinstance( item[1]['nproc'], int ) else 1

            if item[1]['buff'] is None :
                pbuff = None
            elif os.path.exists( item[1]['buff'] ) and os.path.isdir( item[1]['buff'] ) :
                pbuff = item[1]['buff']
            else :
                utility.ExitFailure( 'Cannot decide buffer file path' )
            run['buff'] = pbuff

            self.__runlist.append( run )
            self.__keys.append( run['key'] )

        os.chdir( self.__workdir )

        return runlist


    #____________________________________________________
    def decodeRunlist( self, path ) :

        if self.__flReady is False :
            return None

        if self.__workdir is None :
            self.__workdir = self.getWorkdir( path )

        data = dict()
        with open( path, 'r' ) as f :
            data = yaml.load( f.read() )

        defset = data['DEFAULT']

        runlist = list()
        for key, parsets in data['RUN'].items() :

            if parsets is None :
                runlist.append( [ key, defset ] )
            else:
                tmp = copy.deepcopy( defset )
                tmp.update( parsets )
                runlist.append( [ key, tmp ] )

        return runlist


    #____________________________________________________
    def makeDataPath( self, path, runno = None ) :

        if self.__flReady is False :
            return None

        data_path = None

        if not os.path.exists( path ) :
            utility.ExitFailure( 'Cannot find file: ' + path )
        else :
            if os.path.isfile( path ) :
                data_path = os.path.realpath( path )
            elif ( os.path.isdir( path )
                   and not runno is None
                   and isinstance( runno, int ) ) :
                base = path + '/run{0:05d}'.format( runno )
                tmp = base + '.dat.gz'
                if not os.path.isfile( tmp ) :
                  tmp = base + '.dat'
                data_path = ( os.path.realpath( tmp )
                              if os.path.isfile( tmp )
                              else utility.ExitFailure( 'Cannot find file: ' + tmp ) )
            else :
                utility.ExitFailure( 'Cannot decide deta file path' )

        return data_path


    #____________________________________________________
    def makeRootPath( self, path, base = None ) :

        if self.__flReady is False :
            return None

        root_path = None

        if not os.path.exists( path ) :
            dir_path = os.path.dirname( path )
            if os.path.exists( dir_path )\
                    and os.path.isdir( dir_path ):
                root_path = os.path.realpath( path )
            else :
                utility.ExitFailure( 'Cannot decide root file path' )
        elif os.path.isfile( path ) :
            root_path = os.path.realpath( path )
        elif os.path.isdir( path ) and not base is None :
            root_path = os.path.realpath( path + '/' + base + '.root' )
        else :
            utility.ExitFailure( 'Cannot decide root file path' )

        return root_path


    #__________________________________________________
    def getNEvents( self, path, runno = None ) :

        if self.__flReady is False :
            return None

        nevents = None

        if os.path.exists( path ) and os.path.isdir( path ) :
            reclog_path = path + '/recorder.log'
            if os.path.exists( reclog_path ) and os.path.isfile( reclog_path ) :
                cand = list()
                freclog = open( reclog_path, 'r' )
                for line in freclog :
                    words = line.split()
                    if len( words ) > 2 and runno == int( words[1] ) :
                        cand.append( words[15] ) if len( words ) > 15 else -1
                freclog.close()

                nevents = int( cand[0] ) if len( cand ) == 1 else None

        return nevents
