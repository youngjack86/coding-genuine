#!/bin/bash
set -e;
set -u;

#############################
# Note: this script is only a demo for large file download monitor, if you need a monitor like this, the download link and other parameters should be modified
#############################

delta=1000000; # this is the initial value for delta downloaded date due the monitor time interval (in bytes)
newSize0=0; # the initial value for data size in bytes
newSize1=0; # the initial value for data size in bytes
while [ 1 ];do
	status=$(tail -n 5 nt-download.err 2>/dev/null | awk '{if(NF==9 && $8!~/M$/){print $8}}' | wc -l); # the status of last several 500K blocks, to get the download speed
	pid=$(ps uxf | grep 'https://ftp.ncbi.nlm.nih.gov/blast/db/nt.gz' | grep -v grep | awk '{print $2}'); # the running downloading thread PID
	running=$(ps uxf | grep 'https://ftp.ncbi.nlm.nih.gov/blast/db/nt.gz' | grep -v grep | wc -l); # whether is running
	if [ $running -eq 1 -a $status -ge 4 -o $delta -lt 500000 ];then # only for running jobs that are in low speed
		echo -e "[$(date)] status is $status and delta download is $delta ($newSize1 - $newSize0), detected low speed downloading, try to restart ...\n";
		echo -e "kill -9 $pid\nnohup wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/nt.gz 1>nt-download.log 2>nt-download.err &\n" | sh;
		sleep 30;
		delta=1000000;
	else
		if [ $running -eq 0 ];then # whether finished (not running, not low speed downloading)
			echo -e "[$(date)] the downloading maybe finished! Exit now ...\n";
			break;
		fi
		echo -e "[$(date)] status is $status and delta download is $delta ($newSize1 - $newSize0), downloading peacefully ...\n";
		newSize0=$(du -sb nt.gz 2>/dev/null | awk '{print $1}');
		sleep 30;
		newSize1=$(du -sb nt.gz 2>/dev/null | awk '{print $1}');
		delta=$[$newSize1-$newSize0];
	fi
done
