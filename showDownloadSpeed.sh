#!/bin/bash

set -u;
set -e;

if [ $# -lt 1 -o $# -gt 2 ];then
	echo -e "\n\tTo run:  sh  $(basename $0)  file  [interval in seconds] \n\n";
	exit 1;
else
	if [ $# -eq 1 ];then
		file1=$1;
		interVar=3;
	else
		file1=$1;
		interVar=$2;
	fi
fi

#size0=$(ls -la $file1 | awk '{print $5}');
size0=$(du -sbL $file1 | awk '{print $1}');

echo -e  "\n----- Testing Started at [ $(date) ] ------\n";
while [ 1 ];do 
  interval=$interVar;
  size1=$(du -sbL $file1 | awk '{print $1}');
  if [ "-"$size0 == "-" ];then 
    size0=$size1;
  else 
    awk -v size1=$size1 -v size0=$size0 -v interV=$interval -v fileX=$file1 '
    BEGIN{
      len1=split("B,KB,MB,GB,TB",unit3,",");
      printf sprintf("\rDownloading \033[01;34m%s\033[0;00m: \033[01;35m%.2f\033[0;00m %s/s [ %'\''d ]                ",fileX,(size1-size0)/interV/(1024**(split(sprintf("%'\''d,",(size1-size0)/interV),tmpArr,",")-2)),unit3[split(sprintf("%'\''d,",(size1-size0)/interV),tmpArr,",")-1],size1-size0)
    }';size0=$size1;
    sleep $interval;
  fi
done
