#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;

my $info = << "INFO";

  ###########################################################
  This script works to compare the per base nucleotide compo-
  sition and pickup new barcode sets;
  ###########################################################
  -h | --help   Print this help infomation;
  -f | --file   The input file for initial comparison;
  -l | --len    The length of barcode, 8 bp by default;
  -c | --gc     The G/C base numer limit within a barcode;
                default: 3-5
  -r | --repeat The maximum limit of continuous bases;
                default: 3, i.e. ≤ 3 bp (3 is acceptable)
  -i | --ihomo  The max num. of per base nucleotide composi-
                tion allowed with the original barcodes;
                default: 4, i.e. ≤ 4 bp
  -n | --nhomo  The max num. of same bases allowed among the
                new designed barcode set;
                default: 4, i.e. ≤ 4 bp
  -t | --numSet Num. of barcode sets to generate, default 10;
  -o | --output Output name, "barcode.{stat,xls}" by default;
  ###########################################################
   zhangyongjian\@novogene.com  Wed, Mar 15, 2017  12:27:30 PM
  ###########################################################

INFO

my %options=('l'=>'8','c'=>'3-5','r'=>'3','i'=>'4','n'=>'4','t'=>'10','o'=>'barcode');

GetOptions(
	'help|h!'   =>  \$options{'h'},
	'file|f=s'  =>  \$options{'f'},
	'len|l=i'   =>  \$options{'l'},
	'gc|c=s'    =>  \$options{'c'},
	'repeat|r=i'=>  \$options{'r'},
	'ihomo|i=i' =>  \$options{'i'},
	'nhomo|n=i' =>  \$options{'n'},
	'numSet|t=i'=>  \$options{'t'},
	'output|o=s'  =>  \$options{'o'},
	) or die $info;

die "\nAn input file is required!\n".$info if(!$options{'f'});
die $info if($options{'h'});

sub readFile{
	#input: $file
	#output: $arr_ref,$hash_ref,$linesHash_ref
	my $filein = shift;
	my (@name_arr,%seqhash,%lineshash);
	open(FH,"<",$filein) or die "[readFile] ".$!;
	while(my $line=<FH>){
		chomp($line);
		$line =~ s/\W+$//;
		my @tmp_arr = split /\t/,$line;
		next if(@tmp_arr != 3);
		push @name_arr,$tmp_arr[0];
		$seqhash{$tmp_arr[0]} = $tmp_arr[2];
		$lineshash{$tmp_arr[0]} = $line;
	}
	close FH;
	return(\@name_arr,\%seqhash,\%lineshash);
}

sub comparePerBaseComposition{
	## This part could be re-writtern to calculate the distance;
	#input: $scalar_id,$scalar_seq,\%hash
	#output: $arr_ref
	die "[comparePerBaseComposition] need a scalar and a hash!\n" if(@_ != 3);
	my $query_id = shift;
	my $query_seq = shift;
	my $hash_ref = shift;
	my @query_seq_arr;
	my @stat;
	if(defined($$hash_ref{$query_id})){
		$$hash_ref{$query_id}=~s/\W//g;
		@query_seq_arr = split //, lc($$hash_ref{$query_id});
		@stat = [(0) x (length($$hash_ref{$query_id})+1)];
	}else{
		$query_seq=~s/\W//gi;
		@query_seq_arr = split //, lc($query_seq);
		@stat = [(0) x (length($query_seq) + 1)];
	}
	foreach my $id (keys %{$hash_ref}){
		chomp($id);
		if($id eq $query_id){	# check whether the query is in this hash
			next;		# skip the query sequence itself
		}else{
			my $count = 0;
			my @tmp_seq_arr;
			if(defined($$hash_ref{$id})){
				@tmp_seq_arr = split //,lc($$hash_ref{$id});
				}else{
				@tmp_seq_arr = split //,lc($id);
				}
			for(0..@query_seq_arr-1){	# count from the left side (may be revised if necessary)
				if($query_seq_arr[$_] eq $tmp_seq_arr[$_]){
					$count ++;
				}
			}
			${$stat[0]}[$count] ++;
		}
	}
	return(\@{$stat[0]});
}

sub constructRawBarcodes{
	#input: \@contentArr,$length,$gc_limit,$repeat_limit
	#output: \%{$rawBarcodeHash{$length}}
	die "[constructRawBarcodes] need \@contentArray, \$length, \$gc_limit and \$repeat_limit!\n" if(@_ != 4);
	my ($contextArr_ref,$length,$gc_limit_str,$repeat_limit) = @_[0..3];
	die "[constructRawBarcodes] \$gc_limit_str:$gc_limit_str is not right!\n" if($gc_limit_str!~/[1-9][0-9]*-[1-9][0-9]*/);
	my @gc_limit = sort {$a<=>$b} (split /-/,$gc_limit_str);
	my $contents = uc(join("",@{$contextArr_ref}));
	my %rawBarcodeHash;
	my ($total,$cached,$skipped) = (0,0,0);
	for my $tmp_len (1..$length){
		foreach my $tmp_base (@{$contextArr_ref}){
			if($tmp_len>1){
				foreach my $tmp_seq (keys %{$rawBarcodeHash{$tmp_len-1}}){
					my $tmp_count_gc = (uc($tmp_seq)=~tr/GC/GC/);
					my $tmp_count_at = (uc($tmp_seq)=~tr/AT/AT/);
					if($tmp_seq=~/([\Q$contents\E])\1{\Q$repeat_limit\E,}/gi or ($tmp_count_gc>$gc_limit[1] or $tmp_count_at>($length-$gc_limit[0]))){
						print "\t\t[constructRawBarcodes] warn: [$tmp_seq] is not allowed\n";
						delete $rawBarcodeHash{$tmp_len-1}->{$tmp_seq};
					}else{
						my $tmp_neo_seq = $tmp_seq.$tmp_base;
						$tmp_count_gc = (uc($tmp_neo_seq)=~tr/GC/GC/);
						$tmp_count_at = (uc($tmp_neo_seq)=~tr/AT/AT/);
						$total ++;
						if($tmp_neo_seq=~/([\Q$contents\E])\1{\Q$repeat_limit\E,}/gi or ($tmp_count_gc>$gc_limit[1] or $tmp_count_at>($length-$gc_limit[0]))){
							if($tmp_len == $length){
								$skipped ++;
							}
							next;
						}else{
							$rawBarcodeHash{$tmp_len}->{$tmp_neo_seq} = undef;
							if($tmp_len == $length){
								$cached ++;
							}
						}
					}
				}
			}else{
				$rawBarcodeHash{$tmp_len}->{$tmp_base} = undef;
			}
		}
	%{$rawBarcodeHash{$tmp_len-1}} = ();	# clean up to save memory
	}
	print "[constructRawBarcodes] $total tried, with $cached cached and $skipped skipped!\n";
	return(\%{$rawBarcodeHash{$length}});
}

sub generateCandidateBarcodes{
	## This part will call the defined functions: comparePerBaseComposition;
	## I use a random generation of seeds
	#input: $oriBarcodeHash_ref,$rawBarcodeHash_ref,$interSetsHomology,$intraSetsHomology
	#output: $num,$neoBarcodeSetsHash_ref
	die "[generateCandidateBarcodes] need \$oriBarcodeHash_ref, \$rawBarcodeHash_ref, \$interSetsHomology, and \$intraSetsHomology\n" if(@_ != 4);
	my ($oriBarcodeHash_ref,$rawBarcodeHash_ref,$interSetsHomo,$intraSetsHomo) = @_[0..3];
	$interSetsHomo += 1;	# take = to consideration
	$intraSetsHomo += 1;	# take = to consideration
	my %neoBarcodeSetsHash;
	srand(`date +%s`);	# generate random seed by time
	my @id_arr = keys %{$rawBarcodeHash_ref};
	my $num = 0;
	while(@id_arr){
		my $random = int(rand(@id_arr));
		my $tmp_determine = 0;
		my $raw_barcode = $id_arr[$random];
		splice(@id_arr,$random,1);
		my $tmp_score_inter = comparePerBaseComposition($raw_barcode,$raw_barcode,$oriBarcodeHash_ref);
		for my $tmp_index ($interSetsHomo..@{$tmp_score_inter}-1){
			if($$tmp_score_inter[$tmp_index]>0){	# if the per-base composition reached the threshold,take into consideration
				$tmp_determine ++;
			}
		}
		next if($tmp_determine > 0);	# timely stop query when the inter-sets-limit exceeds
		if(%neoBarcodeSetsHash){
			my $tmp_score_intra = comparePerBaseComposition($raw_barcode,$raw_barcode,\%neoBarcodeSetsHash);
			for my $tmp_index ($intraSetsHomo..@{$tmp_score_intra}-1){
				if($$tmp_score_intra[$tmp_index]>0){
					$tmp_determine ++;
				}
			}
		}
		next if($tmp_determine > 0);	# skip saving this barcode when conflict with the existing barcodes
		$neoBarcodeSetsHash{$raw_barcode} = $raw_barcode;
		$num ++;
	}
	return($num,\%neoBarcodeSetsHash);
}

sub main{
	## This is the main function of this script, which works to call the sub-functions and write to files
	#input: global variables
	#output: output files
	my ($barcodeID_arr_ref,$oriBarcodeHash_ref,$linesHash_ref) = readFile($options{'f'});
	open(OUT,">",$options{'o'}.".stat") or die "[main] error in output the barcode statistics\n\t ".$!;
	print OUT (".\t" x 2).".";
	for(0..length($$oriBarcodeHash_ref{$$barcodeID_arr_ref[0]})){
		print OUT "\tN".(length($$oriBarcodeHash_ref{$$barcodeID_arr_ref[0]}) - $_);
	}
	print OUT "\n";
	foreach my $barcodeid (@{$barcodeID_arr_ref}){
		my $compare_arr_ref = comparePerBaseComposition($barcodeid,$$oriBarcodeHash_ref{$barcodeid},$oriBarcodeHash_ref);
		print OUT $$linesHash_ref{$barcodeid};
		for my $index (0..@{$compare_arr_ref}-1){
			print OUT "\t$$compare_arr_ref[ @{$compare_arr_ref} - 1 - $index ]";
		}
		print OUT "\n";
	}
	close OUT;

	print "\n[main] now try to construct raw barcode hash\n\n";
	my @bases = ("A","T","G","C");
	my $rawBarcodeHash_ref = constructRawBarcodes(\@bases,$options{'l'},$options{'c'},$options{'r'});
	open(FINAL,">",$options{'o'}.".xls") or die "[main] error in write new barcode sets\n\t".$!;
	print FINAL "Program started at ".`date`."\n";
	for(1..$options{'t'}){
		my ($barcodeNum,$neoBarcodeSetHash_ref) = generateCandidateBarcodes($oriBarcodeHash_ref,$rawBarcodeHash_ref,$options{'i'},$options{'n'});
		print FINAL "RUN $_:\t$barcodeNum barcodes in this set besides the original barcodes!\n";
		print "\n[main] Round $_: $barcodeNum barcodes in this set besides the original barcodes!\n";
		foreach my $barcodeF (keys %{$neoBarcodeSetHash_ref}){
			print FINAL "$barcodeF\n";
		}
		print FINAL ("-" x 20)."\n";
	}
	print FINAL "Finished at ".`date`."\n";
	close FINAL;
	print "\n[main] all finished! The output files are [$options{o}.stat] and [$options{o}.xls]!\n";
}

main();
