#!/usr/bin/env python3.5

#____________________________________________________

__author__  = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.0'
__date__    = '2 April 2019'

#____________________________________________________

#__________________________________________________
#
# Singleton Class
#__________________________________________________

class Singleton( type ) :

    _instances = {}

    def __call__( cls, *args, **kwargs ) :

        if cls not in cls._instances :
            cls._instances[cls] = super( Singleton, cls ).__call__( *args, **kwargs )
        return cls._instances[cls]
