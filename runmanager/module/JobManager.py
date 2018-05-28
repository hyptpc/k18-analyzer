#!/usr/bin/env python3.5

#____________________________________________________

__author__  = 'Y.Nakada <nakada@km.phys.sci.osaka-u.ac.jp>'
__version__ = '2.1'
__date__    = '12 May 2018'

#____________________________________________________

import os
import sys
import shutil
import time
import copy
import shlex
import subprocess

import xml.etree.ElementTree

#____________________________________________________

SCRIPT_DIR = os.path.dirname( os.path.abspath( sys.argv[0] ) )
XML_NAMESPACE = 'http://www.w3.org/2001/XMLSchema-instance'
SCHEMA_LOC_ATTRIB = 'noNamespaceSchemaLocation'

#____________________________________________________

class JobManager :

    #__________________________________________________
    def __init__( self,
                  tag,
                  key,
                  fexec_path, fconf_path, fdata_path, fout_path,
                  div_unit = -1, nevents = None ) :

        self.__flDone  = False
        self.__finStat = None

        self.__tag = tag

        self.__key = key

        self.__stime = time.time()
        self.__ctime = time.time()
        self.__diffTime = 0

        self.__baseName = self.__tag + '_' \
                          + os.path.splitext( os.path.basename( fout_path ) )[0]

        self.__fExecPath     = fexec_path
        self.__fConfPath     = fconf_path
        self.__fDataPath     = fdata_path
        self.__fOutPath      = fout_path
        self.__fPreFetchPath = None
        self.__fUnpackPath   = None
        self.__fSchemaPath   = None
        self.__fLogPath      = None

        self.__nEvents = nevents
        self.__divUnit = div_unit

        self.__option = '-q s'    # for short que

        self.__elemList    = list()
        self.__jobIdList   = list()
        self.__fConfList   = list()
        self.__fUnpackList = list()
        self.__fOutList    = list()
        self.__fLogList    = list()

        self.__jobStat = dict()

        self.__procMerge = None

        self.__makeFLog()
        self.__makeFPreFetch()

        self.__dumpInitInfo()

        self.__makeElem()


    #__________________________________________________
    def getFinalStatus( self ) :

        if self.__flDone is False :
            self.checkMerged()

        return self.__finStat


    #__________________________________________________
    def getKey( self ) :

        return self.__key


    #__________________________________________________
    def getInfo( self ) :

        tag_exec = os.path.basename( self.__fExecPath )
        tag_conf = os.path.basename( self.__fConfPath )
        tag_data = os.path.basename( self.__fDataPath )
        tag_out  = os.path.basename( self.__fOutPath )

        return self.__key, self.__nEvents, self.__divUnit,\
               tag_exec, tag_conf, tag_data, tag_out


    #__________________________________________________
    def killBjob( self, jid = None ) :

        if jid is None :
            self.__updateJobStat()
            fl = False
            for key, stat in self.__jobStat.items() :
                if stat in [ 0, 1 ] :
                    sys.stdout.write( 'Killing job #%d\n' % key )
                    tmp = subprocess.run( [ 'bkill', str( key ) ],
                                          stdout = subprocess.DEVNULL,
                                          stderr = subprocess.DEVNULL )
                    fl = True
                    buff = 'bsub process was killed [JID: %d]' % key
                    self.__dumpLog( 'killBjob', buff )
            if fl :
                self.__dumpLog( None,  64 * '_' )

        else :
            if jid in self.__jobIdList :
                self.__updateJobStat( jid )
                if self.__jobStat[jid] in [ 0, 1 ] :
                    sys.stdout.write( 'Killing job #%d\n' % jid )
                    tmp = subprocess.run( [ 'bkill', str( jid ) ],
                                          stdout = subprocess.DEVNULL,
                                          stderr = subprocess.DEVNULL )
            buff = 'bsub process was killed [JID: %d]' % jid
            self.__dumpLog( 'killBjob', buff )
            self.__dumpLog( None,  64 * '_' )

    #__________________________________________________
    def killMerge( self ) :

        self.checkMerged()
        if not self.__procMerge is None \
           and self.__procMerge.poll() is None :

            sys.stdout.write( 'Killing merging process\n' )
            self.__procMerge.kill()
            buff = 'merging process was killed'
            self.__dumpLog( 'killMerge', buff )
            self.__dumpLog( None,  64 * '_' )


    #__________________________________________________
    def killAll( self ) :

        self.killBjob()
        self.killMerge()


    #__________________________________________________
    def getStartTime( self ) :

        return time.ctime( self.__stime )


    #__________________________________________________
    def getDiffTime( self ) :

        self.__ctime = time.time()
        self.__diffTime = self.__ctime - self.__stime

        return self.__diffTime


    #__________________________________________________
    def execute( self ) :

        self.__makeFConflist()
        self.__makeFOutlist()
        self.__makeFLoglist()

        for item in zip( self.__elemList,
                         self.__fConfList,
                         self.__fOutList,
                         self.__fLogList ) :

            if os.path.exists( item[3] ):
                os.remove( item[3] )

            bsubcomm = 'bsub' + ' ' \
                       + self.__option + ' ' + \
                       '-o' + ' ' + item[3] + ' ' \
                       '-a \"prefetch (' + self.__fPreFetchPath + ')\"'
            command =  shlex.split( bsubcomm )

            target = [ self.__fExecPath,
                       item[1],
                       self.__fDataPath,
                       item[2] ]
            command.extend( target )

            try :
                proc = subprocess.run( command,
                                       stdout = subprocess.PIPE,
                                       stderr = subprocess.PIPE,
                                       check = True )
            except subprocess.CalledProcessError as e :
                sys.stderr.write( 'ERROR: command \'%s\' returned error code (%d)\n' \
                                  % ( ' '.join( e.cmd ), e.returncode ) )
                sys.stderr.write( proc.stderr )

            jid  = readJobId( proc.stdout )
            self.__jobIdList.append( jid )
            self.__jobStat[jid] = 0

        self.__dumpLog( 'job', self.__jobIdList )
        self.__dumpLog( None,  64 * '_' )


    #__________________________________________________
    def __dumpLog( self, key = None, buff = None ) :

        with open( self.__fLogPath, 'a' ) as flog :

            if key is None :
                if buff is None :
                    flog.write( '\n' )
                else :
                    flog.write( str( buff ) + '\n' )
            else :
                if buff is None :
                    flog.write( str( key ) + ':\n' )
                elif isinstance( buff, list ) or isinstance( buff, tuple ) :
                    for i, item in enumerate( buff ) :
                        flog.write( ('%s(%d):' % ( key, i )).ljust(16) )
                        flog.write( str( item ) + '\n' )
                else :
                    flog.write( ('%s:' % key).ljust( 16 ) )
                    flog.write( str( buff ) + '\n' )


    #__________________________________________________
    def __dumpInitInfo( self ) :

        self.__dumpLog( 'start',    time.ctime( self.__stime ) )
        self.__dumpLog( None,       64 * '_' )
        self.__dumpLog( 'key',      self.__key )
        self.__dumpLog( None,       64 * '_' )
        self.__dumpLog( 'bin',      self.__fExecPath )
        self.__dumpLog( 'conf',     self.__fConfPath )
        self.__dumpLog( 'data',     self.__fDataPath )
        self.__dumpLog( 'unit',     self.__divUnit )
        self.__dumpLog( 'nevent',   self.__nEvents )
        self.__dumpLog( 'prefetch', self.__fPreFetchPath )
        self.__dumpLog( 'out',      self.__fOutPath )
        self.__dumpLog( None,       64 * '_' )


    #__________________________________________________
    def setOption( self, option ) :

        self.__option = option


    #__________________________________________________
    def __updateJobStat( self, jid = None ) :

        if jid is None :

            for key in self.__jobStat.keys() :

                command = 'bjobs %d' % key
                try :
                    proc = subprocess.run( shlex.split( command ),
                                           stdout = subprocess.PIPE,
                                           stderr = subprocess.PIPE,
                                           check = True )
                except subprocess.CalledProcessError as e :
                    sys.stderr.write( 'ERROR: command \'%s\' returned error code (%d)\n' \
                                      % ( ' '.join( e.cmd ), e.returncode ) )
                    sys.stderr.write( proc.stderr )

                buff = proc.stdout.splitlines()[1]
                self.__jobStat[key] = readBsubStatus( key, buff.decode() )

        else :

            if jid in self.__jobStat.keys() :

                command = 'bjobs %d' % jid
                try :
                    proc = subprocess.run( shlex.split( command ),
                                           stdout = subprocess.PIPE,
                                           stderr = subprocess.PIPE,
                                           check = True )
                except subprocess.CalledProcessError as e :
                    sys.stderr.write( 'ERROR: command \'%s\' returned error code (%d)'\
                                      % ( ' '.join( e.cmd ), e.returncode ) )
                    sys.stderr.write( proc.stderr )

                buff = proc.stdout.splitlines()[1]
                self.__jobStat[jid] = readBsubStatus( jid, buff.decode() )


    #__________________________________________________
    def getJobResult( self, jid = None ) :

        fl = None

        if jid is None :
            self.__updateJobStat()
            if 0 < len( self.__jobStat ) :
                if 3 in self.__jobStat.values() : # in the case of 'EXIT'
                    for key, value in self.__jobStat.items() :
                        if value == 3 :
                            sys.stderr.write( 'ERROR: Process has unexpectedly exited [JID: %d]\n' \
                                              % key )
                    fl = False
                elif len( self.__jobStat ) == list( self.__jobStat.values() ).count( 2 ) : # in the case of 'DONE'
                    fl = True
        else :
            self.__updateJobStat( jid )
            if jid in self.__jobStat.keys() :
                if 3 == self.__jobStat[jid] : # in the case of 'EXIT'
                     fl = False
                elif 2 == self.__jobStat[jid] == 2 : # in the case of 'DONE'
                    fl = True

        return fl


    #__________________________________________________
    def getNSegs( self ) :
        return len( self.__elemList )


    #__________________________________________________
    def getProgress( self ) :

        numer = 0

        for jid in self.__jobStat.keys() :
            if self.__jobStat[jid] in [ 2, 3 ] :
                numer += 1

        return numer


    #__________________________________________________
    def __makeElem( self ) :

        if self.__divUnit <= 0 :
            nsegs = 1
        else :
            nsegs = 10 if self.__nEvents is None else\
                    self.__nEvents // self.__divUnit + 1

        for i in range( nsegs ) :
            self.__elemList.append( self.__baseName + '_%d' % i )

        self.__dumpLog( 'elem', self.__elemList )
        self.__dumpLog( None,   64 * '_' )

        return nsegs


    #__________________________________________________
    def __makeFLog( self ) :

        dlog = SCRIPT_DIR + '/log'
        if not os.path.exists( dlog ) :
            os.mkdir( dlog )

        self.__fLogPath = dlog + '/' + self.__baseName + '.log'

        if os.path.exists( self.__fLogPath ) :
            os.remove( self.__fLogPath )


    #__________________________________________________
    def __makeFPreFetch( self ) :

        dPreFetch = SCRIPT_DIR + '/prefetch'
        if not os.path.exists( dPreFetch ) :
            os.mkdir( dPreFetch )

        self.__fPreFetchPath = dPreFetch + '/' + self.__baseName + '.pf'

        with open( self.__fPreFetchPath, 'w' ) as f :
            f.write( self.__fDataPath )


    #__________________________________________________
    def __makeFConflist( self ) :

        dummy_dconf = SCRIPT_DIR + '/conf'
        if not os.path.exists( dummy_dconf ) :
            os.mkdir( dummy_dconf )

        with open( self.__fConfPath, 'r' ) as f :
            lines = f.readlines()

        buff = list()
        buff_unpack = list()

        for line in lines :
            if line[0] == '#' :
                continue
            elif 0 == line.find( 'UNPACK:' ) :
                words = line.split()
                buff_unpack.append( words[1] )
            else :
                buff.append( line )

        if 1 == len( buff_unpack ) :
            self.__fUnpackPath = buff_unpack[0]
        else :
            sys.stderr.write( 'ERROR: Invalid unpacker declaration was found in %s'\
                              % self.__fConfPath )
            sys.exit( 1 )

        for item in self.__elemList :

            path_conf = dummy_dconf + '/' + item + '.conf'
            f = open( path_conf, 'w' )
            path_unpack = '%s/unpack/%s.xml' % ( SCRIPT_DIR, item )
            f.write( 'UNPACK: ' + path_unpack + '\n' )
            f.writelines( buff )
            f.close()

            self.__fConfList.append( path_conf )
            self.__fUnpackList.append( path_unpack )

        self.__dumpLog( 'conf',   self.__fConfList )
        self.__dumpLog( 'unpack', self.__fUnpackList )
        self.__dumpLog( None,     64 * '_' )

        self.__generateFUnpack()



    #__________________________________________________
    def __generateFUnpack( self ) :

        dummy_dunpack = SCRIPT_DIR + '/unpack'
        if not os.path.exists( dummy_dunpack ) :
            os.mkdir( dummy_dunpack )

        tree = xml.etree.ElementTree.parse( self.__fUnpackPath )

        self.__generateSchema( os.path.dirname( self.__fUnpackPath ),
                               tree )

        for i, item in enumerate( self.__fUnpackList ) :
            tmp = copy.deepcopy( tree )

            root = tmp.getroot()
            key = '{%s}%s' % ( XML_NAMESPACE, SCHEMA_LOC_ATTRIB )
            root.set( key, self.__fSchemaPath )

            for node in tmp.findall( 'control/skip' ) :
                node.text = str( i * self.__divUnit )
            for node in tmp.findall( 'control/max_loop' ) :
                node.text = str( -1 if i == len( self.__fUnpackList ) - 1 \
                                 else self.__divUnit )

            tmp.write( item )


    #__________________________________________________
    def __generateSchema( self, path, tree ) :

        root = tree.getroot()
        tmp = root.get( '{%s}%s' % ( XML_NAMESPACE, SCHEMA_LOC_ATTRIB ) )

        orig = path + '/' + tmp
        self.__fSchemaPath = SCRIPT_DIR + '/unpack/'\
                             + self.__baseName + '.xsd'
        shutil.copyfile( orig, self.__fSchemaPath )


    #__________________________________________________
    def __makeFOutlist( self ) :

        for item in self.__elemList :
            path = os.path.dirname( self.__fOutPath ) + '/' + item + '.root'
            self.__fOutList.append( path )

        self.__dumpLog( 'out', self.__fOutList )
        self.__dumpLog( None,  64 * '_' )

    #__________________________________________________
    def __makeFLoglist( self ) :

        dlog = SCRIPT_DIR + '/log'
        if not os.path.exists( dlog ) :
            os.mkdir( dlog )

        for item in self.__elemList :
            path = dlog + '/' + item + '.log'
            self.__fLogList.append( path )

        self.__dumpLog( 'log', self.__fLogList )
        self.__dumpLog( None,  64 * '_' )


    #__________________________________________________
    def mergeFOut( self ) :

        command = shlex.split( 'hadd -ff' )
        command.append( self.__fOutPath )
        command.extend( self.__fOutList )

        self.__procMerge = subprocess.Popen( command,
                                             stdout = subprocess.PIPE,
                                             stderr = subprocess.PIPE )


    #__________________________________________________
    def checkMerged( self ) :

        if not self.__procMerge is None\
           and not self.__procMerge.poll() is None :

            if self.__procMerge.poll() != 0 :
                sys.stderr.write( 'ERROR: Merging root file failed\n' )
                self.__finStat = False
            else :
                self.__finStat = True

            self.__flDone = True

            outs, errs = self.__procMerge.communicate()

            self.__dumpLog( 'hadd_out', None )
            self.__dumpLog( None, outs.decode() )
            self.__dumpLog( 'hadd_err', None )
            self.__dumpLog( None, errs.decode() )
            self.__dumpLog( None, 64 * '_' )

        else :
            self.__flDone = False


    #__________________________________________________
    def clear( self ) :

        self.clearFOut()
        self.clearFConf()
        self.clearFUnpack()


    #__________________________________________________
    def clearAll( self ) :

        self.clearFOut()
        self.clearFConf()
        self.clearFUnpack()

        self.clearAllLog()


    #__________________________________________________
    def clearFOut( self ) :

        for item in self.__fOutList :
            if os.path.exists( item ) :
                os.remove( item )

    #__________________________________________________
    def clearFPreFetch( self ) :

        if not self.__fLogPath is None :
            if os.path.exists( self.__fPreFetchPath ) :
                os.remove( item )


    #__________________________________________________
    def clearAllLog( self ) :

        self.clearFLog()
        self.clearMainLog()


    #__________________________________________________
    def clearFLog( self ) :

        for item in self.__fLogList :
            if os.path.exists( item ) :
                os.remove( item )


    #__________________________________________________
    def clearMainLog( self ) :

        if not self.__fLogPath is None :
            if os.path.exists( self.__fLogPath ) :
                os.remove( self.__fLogPath )


    #__________________________________________________
    def clearFConf( self ) :

        for item in self.__fConfList :
            if os.path.exists( item ) :
                os.remove( item )


    #__________________________________________________
    def clearFUnpack( self ) :

        if not self.__fSchemaPath is None :
            if os.path.exists( self.__fSchemaPath ) :
                os.remove( self.__fSchemaPath )

        for item in self.__fUnpackList :
            if os.path.exists( item ) :
                os.remove( item )


#__________________________________________________

def readJobId( buff ) :

    words = buff.split()
    jid = int( words[1][1:-1] )

    return jid


#__________________________________________________

def readBsubStatus( jid, buff ) :

    status = -1

    words = buff.split()
    if int( words[0] ) == jid :
        if words[2] == 'PEND' :
            status = 0
        elif words[2] == 'RUN' :
            status = 1
        elif words[2] == 'DONE' :
            status = 2
        elif words[2] == 'EXIT' :
            status = 3

    return status
