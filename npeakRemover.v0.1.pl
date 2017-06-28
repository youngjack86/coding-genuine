#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $usage = << "HELP";

#################################################################
This script works to trim or drop N-containing reads from fastq  
files according to user-defined threshold and strategy;
=================================================================
To RUN:								
   perl  < npeakRemover.pl >  [hfrstmopkz]  < fq > .. 
=================================================================
\033[00;34m
-h | --help    Print this help infomation;

-f | --fq1     The fastq file of Read1;

-r | --fq2     The fastq file of Read2;

-s | --offset  The number of bases offset to skip low-qual region;
	       Note: 10 bp by default
-t | --thresh  The threshold (by %) of per-base N% along reads;    
               e.g. -t 2, >2% N-peak treated!
               Default: -t 0.5, ≤ 0.5 % N-peak will pass.

-o | --output  Output file name prefix, as original by default;

-p | --path    Target path for output, current dir by default;

-k | --keepdr  Output all to-be-dropped/trimmed-reads to file; 
	       -k to keep, no args required,discard by default.
-m | --manpks  The user-specified peaks, seperated by "," 
	       Only specified peaks will be treated, -t acceptable
-z | --zipped  Force output gzipped, as input format if it is not
               specified.   no: for unzipped      yes: for zipped
-g | --args    The search string or pattern. || "N"
               regex supported, but input should in " "
\033[00;00m
=================================================================
\033[01;35m
Author: zhangyongjian\@novogene.com                   June 20, 2017

* For personal use only, not for pipeline or intensified tests. *
\033[00;00m
#################################################################

HELP

my %options;
sub initArgs{
### Get global options
	GetOptions (
		"h|help!" 	=> \$options{h},
		"f|fq1=s" 	=> \$options{f},
		"r|fq2=s" 	=> \$options{r},
		"m|manpks=s"	=> \$options{m},
		"s|offset=i"	=> \$options{s},
		"t|thresh=f"	=> \$options{t},
		"o|output=s"	=> \$options{o},
		"p|path=s"	=> \$options{p},
		"k|keepdr!"	=> \$options{k},
		"z|zipped=s"	=> \$options{z},
		"g|args=s"	=> \$options{g},
	);
	#die $usage unless($options{f} and $options{a} and ($options{n} or $options{t}) and !$options{h});
	die $usage unless($options{f} and !$options{h});
### Initialize key parameters
	if(!defined($options{z})){
		$options{z} = 'yes';
	}
	if(!defined($options{o})){
		my @main_tmp1 = split /\//,$options{f};
		$options{o} = pop @main_tmp1;
		$options{o} =~ s/\.fq(\.gz)$//g;
		$options{o} =~ s/_[12]$//g;
	}
	print "[initArgs]: The output prefix will be [ $options{o} ]!\n";
	my ($main_prefix,undef) = getPrefixAndSuffix($options{f});
	if(!defined($options{p})){
		chomp($options{p} = `pwd`);
		$options{p} .= "/".$main_prefix;
	}
	print "[initArgs]: The output path will be [ $options{p} ]!\n";
	(!system("mkdir -p $options{p}") || die "[initArgs]: make output path dir error, ".$!) if(! -e $options{p});
	if(!defined($options{r})){
		print "[initArgs]: Read2 fastq file not defined! Program will treat the input as single-end sequencing results!\n";
	}else{
		die "[initArgs]: [ $options{r} ] file not exist!\n" unless(-e $options{r});
	}
### The threshold [ -t ] and the cutoff [ -c ] and the offset [ -s ]
	if(!defined($options{t})){
		$options{t} = 0.5/100;
	}elsif($options{t}=~/^[1-9]?[0-9]+\.*[0-9]+$/ or $options{t}=~/^[1-9][0-9]?$/){
		$options{t} = $options{t}/100;
	}else{
		die "[initArgs]: [ $options{t} ] is not allowed!\n";
	}
	if(defined($options{m})){
		my @init_arr = split /,/,$options{m};
		if(@init_arr>0){
			$options{m} = join(",",@init_arr);
			$options{m} =~ s/^,//;
			$options{m} =~ s/[^,0-9]//g;
		}else{
			die "[initArgs]: [ $options{m} ] is not allowed!\n";
		}
	}
	if(!defined($options{s})){
		$options{s} = 10;
	}elsif($options{s}>250 or $options{s}<0){
		die "[initArgs]: [ $options{s} ] is not allowed! Reset it in 0-250 bp please!\n";
	}
	if(!defined($options{g})){
		$options{g} = 'N';
	}
}
sub fastqNparser{
## This function works only for one read of either read1/read2
# input: fR_seq_arr_ref (from fastqReader function's last output, for only one read)
# output: INFO_ARR: info_hash_ref (readlength => length, num_Ns => number_of_Ns, pos_arr_ref => N_pos_arr, pos_hash_ref => N_pos-content_hash )
	my $fNp_seq_arr_ref = shift;
	my $fNp_pattern = $options{g};
	if(@_){
		$fNp_pattern = $_[0];
	}
	my $fNp_totalN = 0;
	my @fNp_info_arr;
	foreach my $fNp_seq (@{$fNp_seq_arr_ref}){
		my %fNp_hash;
		#print "[][][] $fNp_seq\n";
		$fNp_hash{'readlength'} = length($fNp_seq);
		$fNp_hash{'num_Ns'} = 0;
		$fNp_hash{'pos_arr_ref'} = ();
		$fNp_hash{'pos_arr_ref'} = ();
		my $fNp_init_offset = 0;
		#print "+++++\n";if($fNp_seq=~/N/i);
		while($fNp_seq =~ /$fNp_pattern/gi){
			$fNp_hash{'num_Ns'} ++;
			$fNp_totalN ++;
			my $fNp_match = $&;
			my $fNp_pos = index($fNp_seq,$fNp_match,$fNp_init_offset);
			push @{$fNp_hash{'pos_arr_ref'}},$fNp_pos;
			$fNp_hash{'pos_hash_ref'}->{$fNp_pos} = $fNp_match;
			$fNp_init_offset = $fNp_pos;
		}
		#print "[*][*][*]\t".$fNp_hash{'num_Ns'}."\t".$fNp_hash{'readlength'}."\n";
		push @fNp_info_arr,\%fNp_hash;
	}
	return($fNp_totalN,\@fNp_info_arr);
}
sub main(){
	## initialize the default
	initArgs;
	my %stat;
	$stat{'total'} = 0;
	$stat{'keep'} = 0;
	$stat{'drop'} = 0;
#=block # shared part
	my ($FHfwd,$FHrvs,$FHfwdw,$FHrvsw,$FHfwdadptw,$FHrvsadptw,$FHfwdDropw,$FHrvsDropw);
	if(isZippedFile($options{f})){
		$FHfwd = FileHandleManager("zr",$options{f},"FHfwd");
		$FHrvs = FileHandleManager("zr",$options{r},"FHrvs") if($options{r});
	}else{
		$FHfwd = FileHandleManager("r",$options{f},"FHfwd");
		$FHrvs = FileHandleManager("r",$options{r},"FHrvs") if($options{r});
	}
	if($options{z} =~ /^no$/i){
		$FHfwdw = FileHandleManager("w",$options{p}."/".$options{o}."_1.fq","FHfwdw");
		$FHrvsw = FileHandleManager("w",$options{p}."/".$options{o}."_2.fq","FHrvsw") if($options{r});
		if($options{k}){
			$FHfwdDropw = FileHandleManager("w",$options{p}."/drop_".$options{o}."_1.fq","FHfwdDropw");
			$FHrvsDropw = FileHandleManager("w",$options{p}."/drop_".$options{o}."_2.fq","FHrvsDropw") if($options{r});
		}
	}elsif(isZippedFile($options{f}) or ($options{z} !~ /^no$/i)){
		$FHfwdw = FileHandleManager("zw",$options{p}."/".$options{o}."_1.fq.gz","FHfwdw");
		$FHrvsw = FileHandleManager("zw",$options{p}."/".$options{o}."_2.fq.gz","FHrvsw") if($options{r});
		if($options{k}){
			$FHfwdDropw = FileHandleManager("zw",$options{p}."/drop_".$options{o}."_1.fq.gz","FHfwdDropw");
			$FHrvsDropw = FileHandleManager("zw",$options{p}."/drop_".$options{o}."_2.fq.gz","FHrvsDropw") if($options{r});
		}

	}
#=cut # End of shared part
#### Set random seed
	srand(1927);
#=block # Neo part
	my %hash_count;
	if(defined($options{m})){
		my @pos_arrs = split /,/,$options{m};
		my %pos_hashs;
		foreach(@pos_arrs){
			$pos_hashs{$_} = 1;
		}
		if(defined($options{r})){
			my ($read1_max,$read2_max)=($options{s}+1,$options{s}+1);
			while(1){
				my @read1_arr = fastqReader($FHfwd);
				my @read2_arr = fastqReader($FHrvs);
				if($read1_arr[0]>0 and $read2_arr[0]>0){	# Read fastq file successfully
					if(${$read1_arr[1]}[0] eq ${$read2_arr[1]}[0]){	# Read1 and read2 are in pairs
						$stat{'total'} ++;
						my ($read1_num,$read1_info_arr_ref) = &fastqNparser($read1_arr[3]);
						my ($read2_num,$read2_info_arr_ref) = &fastqNparser($read2_arr[3]);
						my $read1_len = ${$$read1_info_arr_ref[0]}{'readlength'};
						my $read2_len = ${$$read2_info_arr_ref[0]}{'readlength'};
						if($read1_num==0 and $read2_num==0){
							fastqWriter($FHfwdw,$read1_arr[2]);
							fastqWriter($FHrvsw,$read2_arr[2]);
							$stat{'keep'} ++;
						}else{
							my ($deter1,$deter2) = (0,0);
							if($read1_num>0 and defined($$read1_info_arr_ref[0]{'pos_hash_ref'}->{$read1_max}) and defined($pos_hashs{$read1_max})){
								$deter1 = &determine($stat{'keep'},$options{t},$read1_max,\%{$hash_count{'read1'}});
							}
							if($read2_num>0 and defined($$read2_info_arr_ref[0]{'pos_hash_ref'}->{$read2_max}) and (defined($pos_hashs{$read2_max}) or defined($pos_hashs{$read2_max+$read1_len}))){
								$deter2 = &determine($stat{'keep'},$options{t},$read2_max,\%{$hash_count{'read2'}});
							}
							if($deter1>0 or $deter2>0){
								if($options{k}){	# Need to keep the N-containing reads
									fastqWriter($FHfwdDropw,$read1_arr[2]);
									fastqWriter($FHrvsDropw,$read2_arr[2]);
								}
								$stat{'drop'} ++;
							}else{
								fastqWriter($FHfwdw,$read1_arr[2]);
								fastqWriter($FHrvsw,$read2_arr[2]);
								$stat{'keep'} ++;
								if($read1_num>0){
									$read1_max = &getMaxPos($read1_max,\%{$hash_count{'read1'}},\%{$$read1_info_arr_ref[0]{'pos_hash_ref'}},$options{s});
								}
								if($read2_num>0){
									$read2_max = &getMaxPos($read2_max,\%{$hash_count{'read2'}},\%{$$read2_info_arr_ref[0]{'pos_hash_ref'}},$options{s});
								}
							}
						}
					}else{
						next;	# Skip the unpaired lines
					}
				}else{
					last;	# Read to file end
				}
			}
		}else{
			my $read1_max=$options{s}+1;
			while(1){
				my @read1_arr = fastqReader($FHfwd);
				if($read1_arr[0]>0){	# Read fastq file successfully
					$stat{'total'} ++;
					my ($read1_num,$read1_info_arr_ref) = &fastqNparser($read1_arr[3]);
					my $read1_len = $$read1_info_arr_ref[0]{'readlength'};
					if($read1_num==0){
						fastqWriter($FHfwdw,$read1_arr[2]);
						$stat{'keep'} ++;
					}else{
						my $deter1 = 0;
						if($read1_num>0 and defined($$read1_info_arr_ref[0]{'pos_hash_ref'}->{$read1_max}) and defined($pos_hashs{$read1_max})){
							$deter1 = &determine($stat{'keep'},$options{t},$read1_max,\%{$hash_count{'read1'}});
						}
						if($deter1>0){
							if($options{k}){	# Need to keep the N-containing reads
								fastqWriter($FHfwdDropw,$read1_arr[2]);
							}
							$stat{'drop'} ++;
						}else{
							fastqWriter($FHfwdw,$read1_arr[2]);
							$stat{'keep'} ++;
							if($read1_num>0){
								$read1_max = &getMaxPos($read1_max,\%{$hash_count{'read1'}},\%{$$read1_info_arr_ref[0]{'pos_hash_ref'}},$options{s});
							}
						}
					}
				}else{
					last;	# Read to file end
				}
			}
		}
	}else{
		if(defined($options{r})){
			my ($read1_max,$read2_max)=($options{s}+1,$options{s}+1);
			while(1){
				my @read1_arr = fastqReader($FHfwd);
				my @read2_arr = fastqReader($FHrvs);
				if($read1_arr[0]>0 and $read2_arr[0]>0){	# Read fastq file successfully
					if(${$read1_arr[1]}[0] eq ${$read2_arr[1]}[0]){	# Read1 and read2 are in pairs
						$stat{'total'} ++;
						my ($read1_num,$read1_info_arr_ref) = &fastqNparser($read1_arr[3]);
						my ($read2_num,$read2_info_arr_ref) = &fastqNparser($read2_arr[3]);
						my $read1_len = $$read1_info_arr_ref[0]{'readlength'};
						my $read2_len = $$read2_info_arr_ref[0]{'readlength'};
						if($read1_num==0 and $read2_num==0){
							fastqWriter($FHfwdw,$read1_arr[2]);
							fastqWriter($FHrvsw,$read2_arr[2]);
							$stat{'keep'} ++;
						}else{
							my ($deter1,$deter2) = (0,0);
							if($read1_num>0 and defined($$read1_info_arr_ref[0]{'pos_hash_ref'}->{$read1_max})){
								$deter1 = &determine($stat{'keep'},$options{t},$read1_max,\%{$hash_count{'read1'}});
							}
							if($read2_num>0 and defined($$read2_info_arr_ref[0]{'pos_hash_ref'}->{$read2_max})){
								$deter2 = &determine($stat{'keep'},$options{t},$read2_max,\%{$hash_count{'read2'}});
							}
							if($deter1>0 or $deter2>0){
								if($options{k}){	# Need to keep the N-containing reads
									fastqWriter($FHfwdDropw,$read1_arr[2]);
									fastqWriter($FHrvsDropw,$read2_arr[2]);
								}
								$stat{'drop'} ++;
							}else{
								fastqWriter($FHfwdw,$read1_arr[2]);
								fastqWriter($FHrvsw,$read2_arr[2]);
								$stat{'keep'} ++;
								if($read1_num>0){
									$read1_max = &getMaxPos($read1_max,\%{$hash_count{'read1'}},\%{$$read1_info_arr_ref[0]{'pos_hash_ref'}},$options{s});
								}
								if($read2_num>0){
									$read2_max = &getMaxPos($read2_max,\%{$hash_count{'read2'}},\%{$$read2_info_arr_ref[0]{'pos_hash_ref'}},$options{s});
								}
							}
						}
					}else{
						next;	# Skip the unpaired lines
					}
				}else{
					last;	# Read to file end
				}
			}
		}else{
			my $read1_max=$options{s}+1;
			while(1){
				my @read1_arr = fastqReader($FHfwd);
				if($read1_arr[0]>0){	# Read fastq file successfully
					$stat{'total'} ++;
					my ($read1_num,$read1_info_arr_ref) = &fastqNparser($read1_arr[3]);
					if($read1_num==0){
						fastqWriter($FHfwdw,$read1_arr[2]);
						$stat{'keep'} ++;
					}else{
						my $deter1 = 0;
						if($read1_num>0 and defined($$read1_info_arr_ref[0]{'pos_hash_ref'}->{$read1_max})){
							$deter1 = &determine($stat{'keep'},$options{t},$read1_max,\%{$hash_count{'read1'}});
						}
						if($deter1>0){
							if($options{k}){	# Need to keep the N-containing reads
								fastqWriter($FHfwdDropw,$read1_arr[2]);
							}
							$stat{'drop'} ++;
						}else{
							fastqWriter($FHfwdw,$read1_arr[2]);
							$stat{'keep'} ++;
							if($read1_num>0){
								$read1_max = &getMaxPos($read1_max,\%{$hash_count{'read1'}},\%{$$read1_info_arr_ref[0]{'pos_hash_ref'}},$options{s});
							}
						}
					}
				}else{
					last;	# Read to file end
				}
			}
		}
	}
#=cut # End of Neo part
	close $FHfwd;
	close $FHfwdw;
	if($options{k}){
		close $FHfwdDropw;
	}
	if($options{r}){
		close $FHrvs;
		close $FHrvsw;
		if($options{k}){
			close $FHrvsDropw;
		}
	}
	print  "\033[01;34m$stat{'total'}\033[0;00m reads processed, with \033[01;34m$stat{'keep'}\033[0;00m kept, \033[01;34m$stat{'drop'}\033[0;00m reads dropped!\n";
}
sub determine($$$$){
# input: num_total threshold pos_max count_hash
# output: boolean [1,0]
	my ($de_num,$de_thres,$de_pos_max,$de_hash_ref) = @_[0..3];
	if($de_num>0 and defined($$de_hash_ref{$de_pos_max})){
		if($$de_hash_ref{$de_pos_max}/$de_num <= $de_thres){
			return 0;
		}else{
			return 1;
		}
	}else{
		return 0;
	}
}
sub getMaxPos($$$$){
# input: ori_max_pos count_hash new_pos_hash offset
# output: max_pos
	my ($gMP_oriMaxPos,$gMP_hash_ref,$gMP_posHash_ref,$gMP_offset) = @_[0..3];
	if(!defined($gMP_offset)){
		$gMP_offset = 0;
	}
	if(!defined($$gMP_hash_ref{$gMP_oriMaxPos})){
		$$gMP_hash_ref{$gMP_oriMaxPos} = 0;
	}
	foreach my $gMP_key (keys %{$gMP_posHash_ref}){
		if(defined($$gMP_hash_ref{$gMP_key})){
			$$gMP_hash_ref{$gMP_key}+=1;
		}else{
			$$gMP_hash_ref{$gMP_key} = 1;
		}
		next if($gMP_key <= $gMP_offset);
		if($$gMP_hash_ref{$gMP_oriMaxPos}<$$gMP_hash_ref{$gMP_key}){
			$gMP_oriMaxPos=$gMP_key;
		}
	}
	return $gMP_oriMaxPos;
}
sub FileHandleManager($$){
# input: r/w/zr/zw, filename
# output: FHref
	my $FHM_mode = shift;
	my $FHM_file = shift;
	my $FHM_name;
	if($FHM_mode eq "r"){
		open($FHM_name,"<",$FHM_file) || die "[FileHandleManager]: Read Plain TXT Mode, ".$!;
	}elsif($FHM_mode eq "w"){
		open($FHM_name,">",$FHM_file) || die "[FileHandleManager]: Write Plain TXT Mode, ".$!;
	}elsif($FHM_mode eq "zr"){
		open($FHM_name,"zcat $FHM_file |") || die "[FileHandleManager]: Read Gzipped Plain TXT Mode, ".$!;
	}elsif($FHM_mode eq "zw"){
		open($FHM_name," | gzip > $FHM_file") || die "[FileHandleManager]: Read Gzipped Plain TXT Mode, ".$!;
	}else{
		die "[FileHandleManager]: unrecognized mode!\n";
	}
	return $FHM_name;
}

sub fastqReader{
# input: reference_of_opened_file_handle, number_of_reads
# output: number_of_seq_cached, ref_to_header_arr(machine_flowcell_tile_x_y, no read info or index info),reference_of_cached_arr(corresponding to the 4 lines),ref_to_sequence_arr
	my $fR_FH = shift;
	my $fR_num = shift;
	defined($fR_num) ? ($fR_num = $fR_num*4) : ($fR_num = 4);
	my (@fR_id,@fR_arr,@fR_seq);
	my $fR_cached_num = 0;
#=block
	while($fR_num > 0){
		my $fR_line = <$fR_FH>;
		last if(!defined($fR_line));
		chomp($fR_line);
#=cut
		if($fR_line =~ /^@/ and ($fR_num%4 == 0)){
			push @fR_arr,$fR_line;
			my @fR_tmp_arr = split(/\s/,$fR_line);
			push @fR_id,$fR_tmp_arr[0];
			$fR_num --;
		}elsif($fR_line =~ /^[AUGCTN]+$/i and ($fR_num%4 == 3)){
			push @fR_arr,$fR_line;
			push @fR_seq,$fR_line;
			$fR_num --;
		}elsif($fR_line =~ /^[^\@AUGCTN]+$/i and ($fR_num%4 == 2)){
			push @fR_arr,$fR_line;
			$fR_num --;
		#}elsif($fR_line =~ /^[^@].+$/i and (length($fR_line) == length($fR_seq[$fR_cached_num])) and ($fR_num%4 == 1)){
		}elsif($fR_line =~ /^.+$/i and (length($fR_line) == length($fR_seq[$fR_cached_num])) and ($fR_num%4 == 1)){
			push @fR_arr,$fR_line;
			$fR_num --;
			$fR_cached_num ++;
		}else{
			warn "[fastqReader]: [$fR_line] unrecognized record encountered!\n";
		}
	}
	$fR_num = undef;
	return($fR_cached_num,\@fR_id,\@fR_arr,\@fR_seq);
}
sub fastqWriter{
# input: reference_of_opened_file_handle, ref_fastq_arr
# output: success or not? number of lines written
	my $fW_FH = shift;
	my $fW_arr_ref = shift;
	my $fW_stat = 0;
	my $fW_num = @{$fW_arr_ref};
	for(1..($fW_num%4)){
		pop @{$fW_arr_ref};
	}
#=block
	foreach(@{$fW_arr_ref}){
		chomp;
		(print $fW_FH $_."\n" and ($fW_stat=1)) || (warn "[fastqWriter]: ".$! and ($fW_stat=0));
	}

#=cut
	$fW_num = @{$fW_arr_ref};
	return($fW_stat,$fW_num);
}

sub isZippedFile($){
# input: file_to_be_test
# output: return 1 if zipped, or return 0
	my $file = shift;
	my $stat = 0;
	if($file =~ /\.gz$/){
		if(-e $file){
			(open(Test,"zcat $file |")  and ($stat = 1)) || die "[isZippedFile]: ".$!;
			close Test if($stat);
		}else{
			$stat = -1;
		}
	}else{
		if(-e $file){
			(open(Test,"zcat $file |")  and ($stat = 2)) || die "[isZippedFile]: ".$!;
			close Test if($stat);
		}else{
			$stat = -2;
		}
	}
	if($stat == 1 or $stat == -1){
		return 1;
		}else{
			return 0;
		}
}

sub getPrefixAndSuffix($){
# input: a string
# output: prefix suffix
	my $PAS_input = shift;
	$PAS_input = (split(/[\/\\]/,$PAS_input))[-1];
	my @PAS_tmp_arr = split /[\._]/,$PAS_input;
	return($PAS_tmp_arr[0],$PAS_tmp_arr[-1]);
}

main();
