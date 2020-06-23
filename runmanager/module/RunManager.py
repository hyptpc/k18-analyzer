#!/usr/bin/env python3.5

#____________________________________________________

__author__  = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.0'
__date__    = '2 April 2019'

#____________________________________________________

import os
import sys
import fcntl
import resource
import psutil
import shutil
import subprocess

import time
import datetime
import copy

import shlex

import tempfile

import json
import configparser
import xml.etree.ElementTree

import utility

#____________________________________________________

import Singleton

#____________________________________________________

rsrc = resource.RLIMIT_NPROC
soft, hard = resource.getrlimit( rsrc )
resource.setrlimit( rsrc, ( hard, hard ) )
MAX_NPROC = hard

rsrc = resource.RLIMIT_NOFILE
soft, hard = resource.getrlimit( rsrc )
resource.setrlimit( rsrc, ( hard, hard ) )
MAX_NOFILE = hard

#____________________________________________________

SCRIPT_DIR = os.path.dirname( os.path.abspath( sys.argv[0] ) )
XML_NAMESPACE = 'http://www.w3.org/2001/XMLSchema-instance'
SCHEMA_LOC_ATTRIB = 'noNamespaceSchemaLocation'
BSUB_RESPONCE = 'Job <JobID> is submitted to queue <queue>'

#____________________________________________________

HSTAGE_PATH = '/ghi/fs02/hstage/requests'
HSTAGE_NMIN = 20
HSTAGE_NMAX = 10000

#____________________________________________________

RUN_PERIOD = 5 # unit: second

#____________________________________________________

#__________________________________________________
#
# BJob Class
#__________________________________________________

class BJob( object ) :

    #__________________________________________________
    def __init__( self, jid ) :

        self.__jid = jid
        # 0: STAGE, 1: PEND, 2: RUN, 3: DONE, 4: EXIT, 4: killed, -1: unknown
        self.__status = None
        self.getStatus()


    #__________________________________________________
    def getJobId( self ) :

        return self.__jid

    #__________________________________________________
    def getStatus( self ) :

        if self.__status is 2 \
           or self.__status is 3 \
           or self.__status is 4 :
            return self.__status

        cmd = 'bjobs {}'.format( self.__jid )
        try :
            proc = subprocess.run( shlex.split( cmd ),
                                   stdout = subprocess.PIPE,
                                   stderr = subprocess.PIPE,
                                   check = True )
        except subprocess.CalledProcessError as e :
            sys.stderr.write( 'ERROR: command \'{}\' returned error code ({})'\
                              .format( ' '.join( e.cmd ), e.returncode ) )
            sys.stderr.write( proc.stderr )

        buff = proc.stdout.splitlines()[1].decode().split()

        status = None
        if int( buff[0] ) == self.__jid :
            if buff[2] ==   'PEND' :
                status = 0
            elif buff[2] == 'RUN' :
                status = 1
            elif buff[2] == 'DONE' :
                status = 2
            elif buff[2] == 'EXIT' :
                status = 3
            else :
                status = -1

            self.__status = status

        return status

    #__________________________________________________
    def kill( self ) :

        self.getStatus()
        if self.__status is 0 \
           or self.__status is 1 :

            cmd = 'bkill {}'.format( self.__jid )
            proc = subprocess.run( shlex.split( cmd ),
                                   stdout = subprocess.DEVNULL,
                                   stderr = subprocess.DEVNULL )
            self.__status = 4


    #__________________________________________________
    @staticmethod
    def readJobId( buff ) :

        jid = None
        fl = True

        words = buff.split()

        for i, item in enumerate( BSUB_RESPONCE.split() ) :
            if i == 1 or i == 6 : continue
            if item != words[i] :
                fl = False

        if fl is True :
            jid = int( words[1][1:-1] )

        return jid



#__________________________________________________
#
# Analysis Job Class
#__________________________________________________

class AnalysisJob( object ) :

    #__________________________________________________
    def __init__( self,
                  manager,
                  tag, conf, out, log ) :

        self.__mainProc = psutil.Process()
        self.__manager = manager

        self.__tag = tag

        self.__conf = conf
        self.__out  = out
        self.__log  = log

        self.__pid  = None
        self.__proc = None

        self.__jid  = None
        self.__bjob = None

        self.__stime = None
        self.__rtime = None

        # True: success, False: failure,
        # 0: process thrown, 1: process return and bsub running, 2: killed, -1: unknown
        self.__status   = None

        # True: success, False: failure, 0: processing, 1: killed, -1: unknown
        self.__statProc = None # UNIX process.
        self.__statBjob = None # bsub process.


    #__________________________________________________
    def getTag( self ) :

        return self.__tag


    #__________________________________________________
    def getProcessId( self ) :

        return self.__pid


    #__________________________________________________
    def getJobId( self ) :

        self.__updateProcessStatus()

        return self.__jid


    #__________________________________________________
    def getRunTime( self ) :

        self.__updateJobStatus()

        if self.__rtime is None :
            return 0
        else :
            return self.__rtime


    #__________________________________________________
    def getStatus( self ) :

        if self.__status is True \
            or self.__status is False \
            or self.__status is 2 :
            return self.__status

        self.__updateStatus()

        return self.__status

    #__________________________________________________
    def __updateStatus( self ) :

        if self.__status is True \
            or self.__status is False \
            or self.__status is 2 :
            return

        self.__updateProcessStatus()
        self.__updateJobStatus()

        if self.__statProc is None :
            pass
        elif self.__statProc is True :
            self.__status = 0
        elif self.__statProc is False :
            self.__status = False
        elif self.__statProc is 1 :
            self.__status = 2
        else :
            self.__status = -1

        if self.__statBjob is None :
            pass
        elif self.__statBjob is True :
            self.__status = True
        elif self.__statBjob is False :
            self.__status = False
        elif self.__statBjob is 0 :
            self.__status = 1
        elif self.__statBjob is 1 :
            self.__status = 2
        else :
            self.__status = -1


    #__________________________________________________
    def __updateProcessStatus( self ) :

        if not self.__statProc is 0 :
            return

        if self.__proc.poll() is None :
            return
        elif self.__proc.poll() == 0 :
            self.__statProc = True
            if self.__bjob is None :
                self.__registerJob()
                self.__stime = time.time()
                self.__rtime = self.__stime
        elif self.__proc.poll() == 1 :
            sys.stderr.write( 'ERROR: bsub command failed at {}\n'\
                              .format( self.__tag ) )
            outs, errs = proc.communicate()
            self.__statProc = False
        else :
            self.__statProc = -1
            # self.__status   = -1


    #__________________________________________________
    def __updateJobStatus( self ) :

        if not self.__statBjob is 0 :
            return

        status = self.__bjob.getStatus()
        if status is None :     # initial
            pass
        elif status is 0 :      # PEND
            pass
        elif status is 1 :      # RUN
            self.__rtime = time.time() - self.__stime
            pass
        elif status is 2 :      # DONE
            self.__statBjob = True
            # self.__status   = True
        elif status is 3 :      # EXIT
            self.__statBjob = False
            # self.__status   = False
        elif status is 4 :      # killed
            pass
        else :                  # unknow case
            self.__statBjob = -1
            # self.__status   = -1


    #__________________________________________________
    def __registerJob( self ) :

        self.__updateStatus()

        if self.__jid is None and self.__statProc is True :
            buff  = self.__proc.communicate()[0]
            self.__jid = BJob.readJobId( buff.decode() )
            self.__bjob = BJob( self.__jid )
            self.__statBjob = 0
            # self.__status   = 1


    #__________________________________________________
    def execute( self ) :

        self.__updateStatus()

        if not self.__statProc is None :
            return

        while True :
            nfds = self.__mainProc.num_fds()
            nproc = len( self.__mainProc.children( recursive=True ) )
            if nfds < MAX_NOFILE and nproc < MAX_NPROC :
                break

        # pf = self.__manager.getPreFetchPath()

        cmd = shlex.split( 'bsub' + ' ' \
                           + '-q' + ' ' + self.__manager.getQueue() + ' ' \
                           + self.__manager.getOption() + ' ' \
                           + '-o' + ' ' + self.__log + ' ' \
                           # + '-a \"prefetch (' + pf + ')\"' \
                           )

        cmd.extend( [ self.__manager.getExecutingPath(),
                      self.__conf,
                      self.__manager.getDataPath(),
                      self.__out ] )

        self.__proc = subprocess.Popen( cmd,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE )

        self.__pid = self.__proc.pid

        self.__statProc = 0
        # self.__status   = 0


    #__________________________________________________
    def kill( self ) :

        self.__updateProcessStatus()
        if self.__statProc is 0 \
           or self.__statProc is -1 :
            sys.stdout.write( 'Killing process [pid: {}]\n'.format( self.__pid ) )
            self.__proc.kill()
            self.__statProc = 1

        self.__updateJobStatus()
        if self.__statBjob is 0 \
           or self.__statBjob is -1 :
            sys.stdout.write( 'Killing job [jid: {}]\n'.format( self.__pid ) )
            self.__bjob.kill()
            self.__statBjob = 1



#__________________________________________________
#
# Single Run Manager Class
#__________________________________________________

class SingleRunManager( object ) :

    #__________________________________________________
    def __init__( self,
                  tag,
                  runinfo ) :

        # None: initial
        # True: success,      False: failure,
        #   10: unstaged,        11: staged,
        #   20: bjob running,    21: bjob complete
        #   30: merging,
        #   99: killed,          -1: unknown
        self.__status = None

        # true: staged, false: unstaged
        self.__statStage = None
        # true: success, false: failure, 0: running, 1: killed
        self.__statBjob  = None
        self.__statMerge = None

        self.__tag = tag

        self.__key           = runinfo['key']
        self.__fExecPath     = runinfo['bin']
        self.__fConfPath     = runinfo['conf']
        self.__fDataPath     = runinfo['data']
        self.__fOutPath      = runinfo['root']

        # self.__fPreFetchPath = None
        self.__fUnpackPath   = None
        self.__fSchemaPath   = None
        self.__fLogPath      = None
        self.__fMergeLogPath = None

        self.__nProc    = runinfo['nproc']
        self.__buffPath = runinfo['buff']

        self.__nEvents = runinfo['nevents']
        self.__divUnit = runinfo['unit']

        self.__queue  = runinfo['queue']
        self.__option = ''

        self.__stime = time.time()
        self.__diffTime = 0

        self.__baseName = self.__tag + '_' \
                          + os.path.splitext( os.path.basename( self.__fOutPath ) )[0]

        dtmp_path = '{}/tmp'.format( SCRIPT_DIR )
        if not os.path.exists( dtmp_path ) :
            os.mkdir( dtmp_path )
        self.__dDummy = tempfile.TemporaryDirectory( dir=dtmp_path )

        self.__elemList    = list()
        self.__fConfList   = list()
        self.__fUnpackList = list()
        self.__fOutList    = list()
        self.__fLogList    = list()
        self.__bjobList    = list()
        self.__jidList     = list()
        self.__jstatList   = list()

        self.__procStage = None
        self.__procMerge = None
        self.__jobMerge  = None

        self.__makeFLog()
        # self.__makeFPreFetch()

        self.__dumpInitInfo()

        self.__makeElem()

        self.__nComplete = 0


    #__________________________________________________
    def setOption( self, option ) :

        self.__option = option


    #__________________________________________________
    def getOption( self ) :

        return self.__option


    #__________________________________________________
    def getKey( self ) :

        return self.__key


    #__________________________________________________
    # def getPreFetchPath( self ) :

    #     return self.__fPreFetchPath


    #__________________________________________________
    def getExecutingPath( self ) :

        return self.__fExecPath


    #__________________________________________________
    def getDataPath( self ) :

        return self.__fDataPath


    #__________________________________________________
    def getQueue( self ) :

        return self.__queue


    #__________________________________________________
    def getNSegs( self ) :

        return len( self.__elemList )


    #__________________________________________________
    def getProgress( self ) :

        return self.__nComplete


    #__________________________________________________
    def getStartTime( self ) :

        return time.ctime( self.__stime )


    #__________________________________________________
    def getDiffTime( self ) :

        self.__updateStatus()

        if not self.__status is True \
           and not self.__status is False \
           and not self.__status is 99 :
            self.__diffTime = time.time() - self.__stime

        return self.__diffTime


    #__________________________________________________
    def getInfo( self ) :

        self.__updateStatus()

        info = dict()
        info['queue'] = self.__queue
        info['nproc'] = self.__nProc
        info['unit']  = self.__divUnit
        info['nev']   = self.__nEvents
        info['bin']   = self.__fExecPath
        info['conf']  = self.__fConfPath
        info['data']  = self.__fDataPath
        info['root']  = self.__fOutPath
        info['time']  = self.getDiffTime()
        info['stat']  = self.getStatus()
        info['nseg']  = self.getNSegs()
        info['prog']  = self.getProgress()

        return { format( self.__key ): info }


    #__________________________________________________
    def getStatus( self ) :

        if self.__status is 99 \
           or self.__status is True \
           or self.__status is False :
            return self.__status

        self.__updateStatus()

        return self.__status


    #__________________________________________________
    def __updateStatus( self ) :

        self.__updateStagingStatus()
        self.__updateJobStatus()
        self.__updateMergingStatus()

        if self.__statStage is None : # initial
            pass
        elif self.__statStage is False : # unstaged
            self.__status = 10
        elif self.__statStage is True : # staged
            self.__status = 11
        elif self.__statStage is -1 : # killed
            pass
        else :  # unknown
            self.__status = -1

        if self.__statBjob is None : # initial
            pass
        elif self.__statBjob is 0 : # running
            self.__status = 20
        elif self.__statBjob is True : # complete
            self.__status = 21
        elif self.__statBjob is False : # failure
            self.__status = False
            return
        elif self.__statBjob is 1 : # killed
            self.__status = 99
            return
        else : # unknown
            self.__status = -1

        if self.__statMerge is None : # initial
            pass
        elif self.__statMerge is 0 : # merging
            self.__status = 30
        elif self.__statMerge is True : # complete
            self.__status = True
            return
        elif self.__statMerge is False : # failure
            self.__status = False
            return
        elif self.__statMerge is 1 : # killed
            self.__status = 99
            return
        else : # unknown
            self.__status = -1


    #__________________________________________________
    def __updateStagingStatus( self ) :

        if self.__statStage is not True :
            if self.isStaged() :
                self.__statStage = True
            else :
                self.__statStage = False


    #__________________________________________________
    def isStaged( self ) :

        if self.__statStage is True :
            return True

        cmd = 'ghils {}'.format( self.__fDataPath )
        try :
            proc = subprocess.run( shlex.split( cmd ),
                                   stdout = subprocess.PIPE,
                                   stderr = subprocess.PIPE,
                                   check = True )
        except subprocess.CalledProcessError as e :
            sys.stderr.write( 'ERROR: command \'{}\' returned error code ({})'\
                              .format( ' '.join( e.cmd ), e.returncode ) )
            sys.stderr.write( proc.stderr )

        buff = proc.stdout.decode().split()[0]

        if buff == 'G' or buff == 'B' :
            return True
        else :
            return False


    #__________________________________________________
    def accessDataStream( self ) :

        if self.__statStage is True :
            return

        if self.__procStage is not None :
            return

        cmd = 'head {}'.format( self.__fDataPath )
        print(cmd)
        try :
            self.__procStage = subprocess.Popen( shlex.split( cmd ),
                                                 stdout = subprocess.DEVNULL,
                                                 stderr = subprocess.DEVNULL )
                                                 # check = True )
            self.__dumpLog( 'pid[stage]', self.__procStage.pid )
            self.__dumpLog( None,  64 * '_' )
        except subprocess.CalledProcessError as e :
            sys.stderr.write( 'ERROR: command \'{}\' returned error code ({})'\
                              .format( ' '.join( e.cmd ), e.returncode ) )
            sys.stderr.write( proc.stderr )


    #__________________________________________________
    def __updateJobStatus( self ) :

        if not self.__statBjob is 0 :
            return

        n_complete = 0
        for i, job in enumerate( self.__bjobList ) :

            stat = job.getStatus()
            if stat is True :
                n_complete += 1
                if self.__jstatList[i] is 1 :
                    self.__dumpLog( 'time[bjob({})]'.format( i ), \
                        self.decodeTime( job.getRunTime() ) )
                    self.__dumpLog( None,  64 * '_' )
            elif stat is False :
                self.__statBjob = False
                self.__dumpLog( 'error', 'error at {}'.format( job.getTag() ) )
                self.__dumpLog( None,  64 * '_' )
            elif stat is 1 :
                if self.__jidList[i] is None :
                    jid = job.getJobId()
                    self.__jidList[i] = jid
                    self.__dumpLog( 'jid[bjob({})]'.format( i ), jid )
                    self.__dumpLog( None,  64 * '_' )

            self.__jstatList[i] = stat

        if len( self.__elemList ) == n_complete :
            self.__statBjob = True

        self.__nComplete = n_complete


    #__________________________________________________
    def mergeFOut( self ) :

        if not self.__statMerge is None :
            return

        if 1 == len( self.__elemList ) :
            os.rename( self.__fOutList[0], self.__fOutPath )
            self.__statMerge = True
            return

        size = 0
        for item in self.__fOutList :
            size += os.path.getsize( item )

        # qOpt = '-q sx' if size > 2 * 1024**3 + 5 * 1024**2 else '-q s'
        qOpt = '-q s'
        pOpt = '-j {}'.format( self.__nProc ) if self.__nProc > 1 else ''
        bOpt = '-d {}'.format( self.__buffPath ) if not self.__buffPath is None else ''

        cmd = shlex.split( 'bsub {} -o {} hadd -ff {} {}'\
                           .format( qOpt, self.__fMergeLogPath, pOpt, bOpt ) )
        cmd.append( self.__fOutPath )
        cmd.extend( self.__fOutList )

        self.__procMerge = subprocess.run( cmd,
                                           stdout = subprocess.PIPE,
                                           stderr = subprocess.PIPE,
                                           check = True )
        buff = self.__procMerge.stdout
        self.__jobMerge = BJob( BJob.readJobId( buff.decode() ) )

        self.__dumpLog( 'jid[merge]', self.__jobMerge.getJobId() )
        self.__dumpLog( None,  64 * '_' )

        self.__statMerge = 0


    #__________________________________________________
    def __updateMergingStatus( self ) :

        if not self.__statMerge is 0 :
            return

        stat = self.__jobMerge.getStatus()
        if stat == 0 or stat == 1 : # PEND or RUN
            return
        elif stat == 2 :    # DONE
            self.__statMerge = True
        elif stat == 3 :    # EXIT
            self.__statMerge = False
            self.__dumpLog( 'error', 'merging error' )
            self.__dumpLog( None,  64 * '_' )
        else :              # unknown
            self.__statMerge = -1


    #__________________________________________________
    def execute( self ) :

        if not self.__statBjob is None :
            return

        self.__makeFConflist()
        self.__makeFOutlist()
        self.__makeFLoglist()
        if self.__statBjob is False :
            return

        for i, elem in enumerate( self.__elemList ) :

            if os.path.exists( self.__fOutList[i] ) :
                os.remove( self.__fOutList[i] )

            job = AnalysisJob( self, elem,
                               self.__fConfList[i],
                               self.__fOutList[i],
                               self.__fLogList[i] )
            job.execute()

            self.__bjobList.append( job )
            self.__jidList.append( None )
            self.__jstatList.append( job.getStatus() )

        for i, item in enumerate( self.__bjobList ) :
            self.__dumpLog( 'pid[bjob({})]'.format( i ), item.getProcessId() )
        self.__dumpLog( None,  64 * '_' )

        self.__statBjob = 0
        # self.__status   = 0

    #__________________________________________________
    def killStaging( self ) :

        if self.__procStage is None :
            return

        if self.__procStage.poll() is None :
            pid = self.__procStage.pid
            self.__procStage.kill()
            buff = 'Process was killed [pid: {}]'\
                   .format( pid )
            self.__dumpLog( 'killStage', buff )
            self.__dumpLog( None,  64 * '_' )
            self.__statStage = -1


    #__________________________________________________
    def killBjob( self ) :

        if self.__statBjob is True :
            return

        for job in self.__bjobList :

            stat = job.getStatus()
            if stat is 0 or stat is 1 :

                job.kill()

                self.__statBjob = 1

                jpid = job.getProcessId() if stat is 0 \
                       else job.getJobId()
                buff = 'Process was killed [jid/pid: {}]'\
                       .format( jpid )
                self.__dumpLog( 'killBjob', buff )

        if self.__statBjob is 1 :
            self.__dumpLog( None,  64 * '_' )


    #__________________________________________________
    def killMerge( self ) :

        self.__updateMergingStatus()

        if self.__statMerge is True \
           or self.__statMerge is False :
            return

        if self.__statMerge is 0 :

            sys.stdout.write( 'Killing merging job\n' )
            self.__jobMerge.kill()
            self.__statMerge = 1
            buff = 'merging job was killed [jid: {}]'\
                   .format( self.__jobMerge.getJobId() )
            self.__dumpLog( 'killMerge', buff )
            self.__dumpLog( None,  64 * '_' )


    #__________________________________________________
    def killAll( self ) :

        self.__updateStatus()

        if self.__status is True :
            return

        self.killStaging()
        self.killBjob()
        self.killMerge()
        self.__updateStatus()


    #__________________________________________________
    def __makeElem( self ) :

        if self.__divUnit <= 0 :
            nsegs = 1
        else :
            nsegs = 1 if self.__nEvents is None else \
                    self.__nEvents // self.__divUnit + 1

        for i in range( nsegs ) :
            self.__elemList.append( '{}_{}'.format( self.__baseName, i ) )

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
    # def __makeFPreFetch( self ) :

    #     self.__fPreFetchPath = self.__dDummy.name + '/' + self.__baseName + '.pf'

    #     with open( self.__fPreFetchPath, 'w' ) as f :
    #         f.write( self.__fDataPath )


    #__________________________________________________
    def __makeFConflist( self ) :

        with open( self.__fConfPath, 'r' ) as f :
            buff = '[dummy]\n' + f.read()

        config = configparser.ConfigParser( delimiters=':', comment_prefixes='#' )
        config.optionxform = lambda option: option
        try :
            config.read_string( buff )
            tmp_unpack = config.get( 'dummy', 'UNPACK' )
        except configparser.Error as e :
            sys.stderr.write( 'ERROR: Invalid format in {}\n'\
                              .format( self.__fConfPath ) )
            self.__statBjob = False
            return

        if os.path.exists( tmp_unpack ) :
            self.__fUnpackPath = os.path.abspath( tmp_unpack )
        else :
            sys.stderr.write( 'ERROR: Cannot find file > {}'\
                              .format( tmp_unpack ) )
            self.__statBjob = False
            return

        for item in self.__elemList :

            path_unpack = self.__dDummy.name + '/' + item + '.xml'
            config['dummy']['UNPACK'] = path_unpack

            path_conf = self.__dDummy.name + '/' + item + '.conf'
            with open( path_conf, 'w' ) as f :
                for option in config.options( 'dummy' ) :
                    f.write( option + ':\t' + config.get( 'dummy', option ) + '\n' )

            self.__fConfList.append( path_conf )
            self.__fUnpackList.append( path_unpack )

        self.__dumpLog( 'conf',   self.__fConfList )
        self.__dumpLog( 'unpack', self.__fUnpackList )
        self.__dumpLog( None,     64 * '_' )

        self.__generateFUnpack()


    #__________________________________________________
    def __generateFUnpack( self ) :

        tree = xml.etree.ElementTree.parse( self.__fUnpackPath )

        self.__generateSchema( os.path.dirname( self.__fUnpackPath ),
                               tree )

        for i, item in enumerate( self.__fUnpackList ) :
            tmp = copy.deepcopy( tree )

            root = tmp.getroot()
            key = '{{{}}}{}'.format( XML_NAMESPACE, SCHEMA_LOC_ATTRIB )
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
        fsrc = root.get( '{{{}}}{}'.format( XML_NAMESPACE, SCHEMA_LOC_ATTRIB ) )
        src = path + '/' + fsrc

        self.__fSchemaPath = self.__dDummy.name + '/' + self.__baseName + '.xsd'
        shutil.copyfile( src, self.__fSchemaPath )


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
            if os.path.exists( path ) :
                os.remove( path )
            self.__fLogList.append( path )

        if 1 < len( self.__elemList ) :
            self.__fMergeLogPath = dlog + '/' + self.__baseName + '_merge.log'
            if os.path.exists( self.__fMergeLogPath ) :
                os.remove( self.__fMergeLogPath )

        self.__dumpLog( 'log', self.__fLogList )
        self.__dumpLog( 'mergelog', self.__fMergeLogPath )
        self.__dumpLog( None,  64 * '_' )


    #__________________________________________________
    def finalize( self, flag = False ) :

        self.__updateStatus()
        self.__dumpLog( 'end', time.ctime( time.time() ) )

        if self.__status is True :
            self.clearAll()
        elif 99 == self.__status :
            if flag is True :
                self.clear()
            else :
                self.clearAll()
        else :
            self.clear()


    #__________________________________________________
    def clear( self ) :

        self.clearFOut()


    #__________________________________________________
    def clearAll( self ) :

        self.clearFOut()
        self.clearAllLog()


    #__________________________________________________
    def clearFOut( self ) :

        for item in self.__fOutList :
            if os.path.exists( item ) :
                os.remove( item )


    #__________________________________________________
    def clearAllLog( self ) :

        self.clearFLog()
        self.clearMainLog()
        self.clearMergeLog()


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
    def clearMergeLog( self ) :

        if not self.__fMergeLogPath is None :
            if os.path.exists( self.__fMergeLogPath ) :
                os.remove( self.__fMergeLogPath )


    #__________________________________________________
    def __dumpLog( self, key = None, msg = None ) :

        buff = str()

        if key is None :
            if not msg is None :
                buff = str( msg )
        else :
            if msg is None :
                buff = str( key ) + ':'
            elif isinstance( msg, list ) or isinstance( msg, tuple ) :
                for i, item in enumerate( msg ) :
                    buff += '{}({}):'.format( key, i ).ljust(16) \
                            + str( item )
                    if i < len( msg ) - 1 :
                        buff += '\n'
            else :
                buff = '{}:'.format( key ).ljust( 16 ) \
                        + str( msg )

        with open( self.__fLogPath, 'a' ) as flog :
            flog.write( buff + '\n' )



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
        self.__dumpLog( 'out',      self.__fOutPath )
        self.__dumpLog( 'nproc',    self.__nProc )
        self.__dumpLog( 'buffPath', self.__buffPath )
        self.__dumpLog( 'queue',    self.__queue )
        # self.__dumpLog( 'prefetch', self.__fPreFetchPath )
        self.__dumpLog( 'dirdummy', self.__dDummy.name )
        self.__dumpLog( None,       64 * '_' )


    #__________________________________________________
    @staticmethod
    def decodeStatus( data ) :

        buff = ''

        stat = data['stat']
        nseg = data['nseg']
        prog = data['prog']

        if stat is None :
            buff = 'init'
        elif stat is 10 :
            buff = 'unstaged'
        elif stat is 11 :
            buff = 'staged'
        elif stat is 20 :
            buff = 'running({}/{})'.format( prog, nseg )
        elif stat is 21 :
            buff = 'running({}/{})'.format( prog, nseg )
        elif stat is 30 :
            buff = 'merging({})'.format( nseg )
        elif stat is 99 :
            buff = 'terminated'.format( nseg )
        elif stat is True :
            buff = 'done({})'.format( nseg )
        elif stat is False :
            buff = 'error'.format( nseg )
        else :
            buff = 'unknown'.format( nseg )

        return buff


    #____________________________________________________
    @staticmethod
    def decodeTime( data ) :

        second = 0
        if isinstance( data, dict ) :
            second = int( data['time'] )
        elif isinstance( data, float ) :
            second = int( data )
        elif isinstance( data, int ) :
            second = int( data )

        hour    = second // 3600
        second -= hour    * 3600
        minute  = second // 60
        second -= minute  * 60

        return '{}:{:02d}:{:02d}'.format( hour, minute, second )


#__________________________________________________
#
# Run Manager Class
#__________________________________________________

class RunManager( metaclass = Singleton.Singleton ) :

    #__________________________________________________
    def __init__( self ) :

        self.stime = None

        self.__runlistMan = None
        self.__fStatPath  = None

        self.__tag   = None
        self.__nruns = None

        self.__fStageList = None

        self.__flReady    = False

        self.__flDone     = list()

        self.__runJobList = list()

        self.__stageList  = list()

    #__________________________________________________
    def setRunlistManager( self, runlistman ) :

        self.__runlistMan = runlistman
        self.__tag   = self.__runlistMan.getTag()
        self.__nruns = self.__runlistMan.getNRuns()


    #__________________________________________________
    def setStatusOutputPath( self, fstat_path ) :

        self.__fStatPath  = fstat_path


    #__________________________________________________
    def initialize( self ) :

        if self.__runlistMan is None \
            or self.__fStatPath is None :
            return

        self.stime = time.time()
        self.__fStageList = HSTAGE_PATH + '/' \
                            + self.__tag + '.lst.'\
                            + datetime.datetime.now().strftime( '%Y%m%d%H%M%S' )

        for i in range( self.__nruns ) :
            run = self.__runlistMan.getRunInfo( i )
            runJob = SingleRunManager( self.__tag, run )
            self.__runJobList.append( runJob )

        self.__stageList = list()
        for runJob in self.__runJobList :
            self.__flDone.append( 0 )
            if not runJob.isStaged() :
                self.__stageList.append( runJob.getDataPath() )

        self.__flReady = True


    #__________________________________________________
    def registerStaging( self ) :

        if self.__flReady is False :
            return

        self.__fStageList = HSTAGE_PATH + '/' \
                            + self.__tag + '.lst.'\
                            + datetime.datetime.now().strftime( '%Y%m%d%H%M%S' )

        if len( self.__stageList ) > HSTAGE_NMAX:
            self.__flReady = False
            utility.ExitFailure('too much runlist!!!')
            return
        if len( self.__stageList ) > HSTAGE_NMIN:
            print('put {}'.format(self.__fStageList))
            with open( self.__fStageList, 'w' ) as f :
                for item in self.__stageList :
                    f.write( item + '\n' )


    #__________________________________________________
    def __runSingleCycle( self ) :

        for i, runJob in enumerate( self.__runJobList ) :

            status = runJob.getStatus()
            if status is None :
                pass
            elif status is 10 :
                if len( self.__stageList ) <= HSTAGE_NMIN:
                    runJob.accessDataStream()
            elif status is 11 :
                runJob.execute()
            elif status is 20 :
                pass
            elif status is 21 :
                runJob.mergeFOut()
            elif status is 30 :
                pass
            elif status is True :
                if self.__flDone[i] == 0 :
                    runJob.finalize()
                    self.__flDone[i] = 1
            elif status is False :
                if self.__flDone[i] == 0 :
                    runJob.finalize()
                    self.__flDone[i] = 1
            else :
                sys.stderr.write( 'ERROR: unknown status was detected.\n' )


    #__________________________________________________
    def run( self ) :

        if self.__flReady is False :
            return

        ptime = time.time()
        while sum( self.__flDone ) != self.__nruns :
            self.__runSingleCycle()
            self.dumpStatus()
            dtime = RUN_PERIOD - ( time.time() - ptime )
            if dtime > 0 :
                time.sleep( dtime )
            ptime = time.time()


    #__________________________________________________
    def dumpStatus( self ) :

        if self.__flReady is False :
            return

        info = dict()
        for runJob in self.__runJobList :
            info.update( runJob.getInfo() )
        buff = json.dumps( info, indent=4 )
        # with open( self.__fStatPath, 'w+' ) as f :
        with open( self.__fStatPath, 'r+' ) as f :
            try :
                fcntl.flock( f.fileno(), fcntl.LOCK_EX )
            except IOError :
                # sys.stderr.write( 'ERROR: I/O error was detected.\n' )
                pass
            else :
                f.write( buff )
                f.truncate()
                f.flush()
            finally :
                fcntl.flock( f.fileno(), fcntl.LOCK_UN )


    #__________________________________________________
    def finalize( self, flag = False ) :

        if self.__flReady is False :
            return

        for runJob in self.__runJobList :
            runJob.finalize( flag )
        self.dumpStatus()


    #__________________________________________________
    def kill( self ) :

        if self.__flReady is False :
            return

        for runJob in self.__runJobList :
            status = runJob.getStatus()
            if not status is True \
                or not status is False :
                runJob.killAll()
        self.dumpStatus()
