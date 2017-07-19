#!/bin/sh
if [ $# -eq 0 ]
then
echo -e "\n\e[01;35m"$(pwd)"\e[0;00m\t\t[ Current Working DIR ]\n"
else
for i in $@;do
	echo $i | awk -F "/" 'BEGIN{
	"pwd" | getline current_path;
} 
{
	if(($1 ne "") || (length($1)>0)){
		pathall = current_path"/"$0;
	}else{
		pathall = $0;
	}
	gsub(/\/+/,"/",pathall);
	len=split(pathall,arr,"/");
	if(pathall~/\.\./){ 
		for(x=1;x<=len;x++){
			tmp = x;
			if(arr[x]~/\.\./){
				arr_cons[x] = 2;
				while(tmp>0){
					if(arr_cons[tmp-1]==2){
						tmp --;
					}else{
						arr_cons[tmp-1] = 2;
						tmp = 0;
					}
				}
			}
			if(arr[x] == "."){
				arr_cons[x] = 2;
			}
		}
	}
	for(y=1;y<=len;y++){
		if(arr_cons[y] == 2){
			arr[y] = "";
		}
		pathFinal = pathFinal"/"arr[y];
	}
	gsub(/\/+/,"/",pathFinal);
	"file "pathFinal | getline retrn;
	len = split(retrn,arr,": ");
	if(arr[len]~/No such file or directory/){
		print "\n\033[05;45m"pathFinal"\033[0;00m\t\t[ Not Exist ]\n\n";
	}else if(arr[len]~/directory$/){
		print "\n\033[01;35m"pathFinal"\033[0;00m\t\t[ A Directory ]\n\n";
	}else{
		split(arr[len],disc,",");
		print "\n\033[01;35m"pathFinal"\033[0;00m\n\n[ A File: "disc[1]" ]\n";
	}
}'
done
fi
