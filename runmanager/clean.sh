#!/bin/sh


cd `dirname $0`

rm -rf prefetch conf unpack log stat
rm -f *~ \#*\# runlist/*~ runlist/\#*\#
rm -rf module/__pycache__
rm -f module/*.pyc
