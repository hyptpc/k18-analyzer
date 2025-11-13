#!/bin/sh

bjobs | while read jid buff; do

    if [ $jid == 'JOBID' ]; then
	continue
    fi

    bkill -r $jid

done

killall hadd
killall run.py
killall dst_run.sh
