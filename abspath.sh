#!/bin/sh
if [ $# -eq 0 ]
	then
	echo -e "\n\e[01;35m"$(pwd)"\e[0;00m\t\t[ Current Working DIR ]\n"
	else
	echo $1 | awk -F "/" 'BEGIN{
			"pwd" | getline current_path;
			len=split(current_path,cp,"/");
			} 
			{
				a=len;
				i=1;
				b="";
				p="";
			if($0~/\.\./) {
				while($i~/\.\./){
					a=a-1;
					i=i+1;
					}
				}else if($0~/\.\//){
				i=i+1;	
				}
			for(j=i;j<=NF;j++){
				if($j~/^\.$/){
					next;
					}
				b=b"/"$j;
				}
			for(k=2;k<=a;k++){
				p=p"/"cp[k];
				}
			gsub(/\/+/,"/",b);
			gsub(/\/+/,"/",p);
			if(($1 ne "") || (length($1)>0)){
				"file "p""b | getline retrn;
				len = split(retrn,arr,": ");
				if(arr[len]~/No such file or directory/){
					print "\n\033[05;45m"p""b"\033[0;00m\t\t[ Not Exist ]\n\n";
					}else if(arr[len]~/directory$/){
					print "\n\033[01;35m"p""b"\033[0;00m\t\t[ A Directory ]\n\n";
					}else{
					split(arr[len],disc,",");
					print "\n\033[01;35m"p""b"\033[0;00m\n\n[ A File: "disc[1]" ]\n";
					}
				}else{
				"file "b | getline retrn;
				len = split(retrn,arr,": ");
                                if(arr[len]~/No such file or directory/){
                                        print "\n\033[05;45m"b"\033[0;00m\t\t[ Not Exist ]\n\n";
                                        }else if(arr[len]~/directory$/){
                                        print "\n\033[01;35m"b"\033[0;00m\t\t[ A Directory ]\n\n";
                                        }else{
                                        split(arr[len],disc,",");
                                        print "\n\033[01;35m"b"\033[0;00m\n\n[ A File: "disc[1]" ]\n";
                                        }
				}
		}'
	fi
