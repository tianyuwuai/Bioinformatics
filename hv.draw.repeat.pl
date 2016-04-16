############################################
# ok, this is a perl script to draw some color bar 
# to show Hv's telomeric repeat and centromeric repeat
# add ltr and cds as pos reference
#########
# first demand can be found in the QQ's email at Apr 2, 2016
#########
# There are 4 tpyes of file, total 16 input files
# Will draw 7x4 color bar

#########
# second demand: 
# need a dotplot and some color bars abut the sides of dotplot
# this script need 2 all-gene-gff files, a BLAST M8 file for dotplot, 2 BLAST M8 files for 2 speci-telmeric,
# 2 BLAST M8 files for 2 speci-centromeric, 2 LTR files for LTR heat map.
#########
# files format:
##### gff:
# chr\tgeneid\tstart_pos\tend_pos
# hv1	Hv1g00001	4070	4229
# 4G	Hv1g00002	4311	4498
# 07g	Hv1g00003	4596	4787
##### M8:
# the same as BLAST M8 result
##### LTR:
# chr\tstart_pos\tend_pos
# hv1	246050	246207
# Hv05g	246703	247050
# 1G	247051	247225
##########
# in this script, user can select special chromosome, filter Blast M8 file only for dotplot 
# can set length of each speci-grid-line, can set length of block, can set 2 boundary values for each color bar
# here, join the chromosome number with "_", and also join the boundary values with "_".
# I also put some parameters ahead the script, for example colors array, input file names. It can be changed conveniently.

# It is a example input command below:
# perl hv.draw.repeat.pl -e 1e-5 -s 200 -n 10 -gs 10000000 -gq 10000000 -cs 1_3_6 -cq 2_5 -b 1000000 -vc 0.4_0.8 -vt 0.05_0.1 -vl 0.0005_0.003 -vg 0.005_0.01
####### parameters:
# -e	e-value of BLAST	default:5e-2
# -s	score of BLAST	default:0
# -n	hitnumber of BLAST	default:50
# -gs	grid line length of subject species	default:10M
# -gq	grid line lenght of query species	default:10M
# -b	length of block	default:1M
# -cs	selected subject chrmosome number	default:all
# -cq	selected query chromosome number	default:all
# -vc	boundary value of centromeric	default:0.04_0.08
# -vt	boundary value of telomeric	default:0.05_0.1
# -vl	boundary value of LTR	default:0.005_0.05
# -vg	boundary value of cds	default:0.005_0.05
############################################

use strict;
use Getopt::Long;
use GD;
use GD::Text::Align;
# set parameter of figure
my $frame_width=2000;
my $frame_height=2000;
my $left_curb=200;
my $top_curb=200;
my $img=GD::Image -> new($frame_width+1.5*$left_curb,$frame_height+1.5*$top_curb);

########################### total parameters below
##set variate
#files
my $tel_sfile; #telomeric file
my $cen_sfile; #centromeric file
my $ltr_sfile; #ltr file
my $gff_sfile="os.chr.gff"; #gff file
my $tel_qfile; #telomeric file
my $cen_qfile; #centromeric file
my $ltr_qfile; #ltr file
my $gff_qfile="os.chr.gff"; #gff file
my $plot_file; #blast file for dotplot
my $output_fig="draw.hv.tel.cen.jpg"; #output figure
#set color
my $white= $img->colorAllocate(255,255,255);
my $black= $img->colorAllocate(0,0,0);
my $mintcream= $img->colorAllocate(245,255,250);
my $dodgerblue= $img->colorAllocate(30,144,255);
my $red= $img->colorAllocate(255,0,0);
my $darkviolet= $img->colorAllocate(148,0,211);
my $orange= $img->colorAllocate(255,165,0);
#screen input
my $poro="p";
my $evalue=5e-2; #blast e-value for dotplot
my $score=0; #blast score for dotplot filter
my $hitnum=50; #blast hitnumber for dotplot
my $grid_s=10000000; #sbjct-speci-grid length
my $grid_q=10000000; #query-speci-grid length
my $chrn_s="all"; #selected chromosome of sbjct
my $chrn_q="all"; #selected chromosome of query
my $myblock=1000000; #length of block
my $val_cen="0.04_0.08"; #value of centromeric
my $val_tel="0.05_0.1"; #boundary value of telomeric
my $val_ltr="0.005_0.05"; #value of LTR
my $val_cds="0.005_0.05"; #value of cds
#get variate
Getopt::Long::GetOptions(
	'e=f' => \$evalue,
	's=i' => \$score,
	'n=i' => \$hitnum,
	'gs=i' => \$grid_s,
	'gq=i' => \$grid_q,
	'b=i' => \$myblock,
	'cs=s' => \$chrn_s,
	'cq=s' => \$chrn_q,
	'vc=s' => \$val_cen,
	'vt=s' => \$val_tel,
	'vl=s' => \$val_ltr,
	'vg=s' => \$val_cds);

########################################

###################### start dotplot part
####### start read gff file
## read query gff and get gene position, then sort hash and get gene order
open(QGFF,$gff_qfile) or die "can't open query gff file.\n";
my %querygene2pos; #gene id -> pos
my %querygene2order; #gene id -> order
my %querygene2chr; #gene id -> chr number
my %querychr2plen; #chr -> chr total pos
my %querychr2olen; #chr -> chr total order
my %tempsort; # temp has for sorting to get order
#read query gff and get gene pos, chr
while (<QGFF>) {
	my @a=split("\t",$_);
	if($a[2]>$a[3]){my $t=$a[3];$a[3]=$a[2];$a[2]=$t;} #make $a[3] as max
	$querygene2pos{$a[1]} = $a[2]; # gene->pos
	$a[0] =~ s/[^0-9]//g;
	my $chr=int($a[0]);
	$querygene2chr{$a[1]} = $chr; #gene->chr
	$tempsort{$a[1]}=$chr."g".$a[2]; #temp for sorting
	#chr->total
	if($querychr2plen{$chr} !~ /^\d/){$querychr2plen{$chr} = $a[3];}
	elsif($querychr2plen{$chr} < $a[3]){$querychr2plen{$chr} = $a[3];}
}
close($gff_qfile);
#sort hash and get order
my $geneNumOnAchr;
my $lastchr = "";
foreach my $key (sort {$tempsort{$a} <=> $tempsort{$b}} keys %tempsort)
{
	my @a=split('g',$tempsort{$key});
	if($a[0] eq $lastchr){$geneNumOnAchr ++;}
    else{$geneNumOnAchr = 1;}
    $querygene2order{$key}=$geneNumOnAchr; #gene->order
    #chr->order
    if($querychr2olen{$a[0]} !~ /^\d/){$querychr2olen{$a[0]} = $geneNumOnAchr;}
    elsif($querychr2olen{$a[0]} < $geneNumOnAchr){$querychr2olen{$a[0]} = $geneNumOnAchr;}
    $lastchr=$a[0];
}
undef(%tempsort);

## read sbjct gff and get gene position, then sort hash and get gene order
open(SGFF,$gff_sfile) or die "can't open sbjct gff file.\n";
my %sbjctgene2pos; #gene id -> pos
my %sbjctgene2order; #gene id -> order
my %sbjctgene2chr; #gene id -> chr number
my %sbjctchr2plen; #chr -> chr total pos
my %sbjctchr2olen; #chr -> chr total order
my %tempsort; # temp has for sorting to get order
#read query gff and get gene pos, chr
while (<SGFF>) {
	my @a=split("\t",$_);
	if($a[2]>$a[3]){my $t=$a[3];$a[3]=$a[2];$a[2]=$t;} #make $a[3] as max
	$sbjctgene2pos{$a[1]} = $a[2]; # gene->pos
	$a[0] =~ s/[^0-9]//g;
	my $chr=int($a[0]);
	$sbjctgene2chr{$a[1]} = $chr; #gene->chr
	$tempsort{$a[1]}=$chr."g".$a[2]; #temp for sorting
	#chr->total
	if($sbjctchr2plen{$chr} !~ /^\d/){$sbjctchr2plen{$chr} = $a[3];}
	elsif($sbjctchr2plen{$chr} < $a[3]){$sbjctchr2plen{$chr} = $a[3];}
}
close($gff_sfile);
#sort hash and get order
my $geneNumOnAchr;
my $lastchr = "";
foreach my $key (sort {$tempsort{$a} <=> $tempsort{$b}} keys %tempsort)
{
	my @a=split('g',$tempsort{$key});
	if($a[0] eq $lastchr){$geneNumOnAchr ++;}
    else{$geneNumOnAchr = 1;}
    $sbjctgene2order{$key}=$geneNumOnAchr; #gene->order
    #chr->order
    if($sbjctchr2olen{$a[0]} !~ /^\d/){$sbjctchr2olen{$a[0]} = $geneNumOnAchr;}
    elsif($sbjctchr2olen{$a[0]} < $geneNumOnAchr){$sbjctchr2olen{$a[0]} = $geneNumOnAchr;}
    $lastchr=$a[0];
}
undef(%tempsort);
####### end read gff file

############### get selected chromosome total length and order
my @querychrlen = ();
my %querychr2order;
#if selected chr is all chr, change chrn_ to all chr number
if ($chrn_q eq "all") {
	my @a=();
	foreach my $key (sort keys %querychr2plen){push(@a,$key);}
	$chrn_q=join("_",@a);
}
my @querychr=split('_',$chrn_q);

my @sbjctchrlen = ();
my %sbjctchr2order;
#if selected chr is all chr, change chrn_ to all chr number
if ($chrn_s eq "all") {
	my @a=();
	foreach my $key (sort keys %sbjctchr2plen){push(@a,$key);}
	$chrn_q=join("_",@a);
}
my @sbjctchr=split('_',$chrn_s);
if($poro eq "p")
{
	for(my $i=0; $i<=$#querychr; $i++)
	{
		$querychr[$i]=~s/[^0-9]//g;
		my $chr = int($querychr[$i]);
		$querychrlen[$#querychrlen+1] = $querychr2plen{$chr}; #sum length
		$querychr2order{$chr} = $#querychrlen; #total order
	}

	for(my $i=0; $i<=$#sbjctchr; $i++)
	{
		$sbjctchr[$i]=~s/[^0-9]//g;
		my $chr = int($sbjctchr[$i]);
		$sbjctchrlen[$#sbjctchrlen+1] = $sbjctchr2plen{$chr}; #sum length
		$sbjctchr2order{$chr} = $#sbjctchrlen; #total order
	}
}
else
{
	for(my $i=0; $i<=$#querychr; $i++)
	{
		$querychr[$i]=~s/[^0-9]//g;
		my $chr = int($querychr[$i]);
		$querychrlen[$#querychrlen+1] = $querychr2olen{$chr}; #sum length
		$querychr2order{$chr} = $#querychrlen; #total order
	}

	for(my $i=0; $i<=$#sbjctchr; $i++)
	{
		$sbjctchr[$i]=~s/[^0-9]//g;
		my $chr = int($sbjctchr[$i]);
		$sbjctchrlen[$#sbjctchrlen+1] = $sbjctchr2olen{$chr}; #sum length
		$sbjctchr2order{$chr} = $#sbjctchrlen; #total order
	}
}
## get sum length
my ($genome1_length, $genome2_length) = (sum_chr_len(@querychrlen), sum_chr_len(@sbjctchrlen));
## get total selected chromosome quantity
my $querychrno = $#querychrlen + 1;
my $sbjctchrno = $#sbjctchrlen + 1;

###### calculate scale ratio
my $scale_ratio1 = $genome1_length/($frame_width); ## horizontal
my $scale_ratio2 = $genome2_length/($frame_height);## vertical

###############################OK, finished all parameters setting

############################### Let's draw dotplot
# draw the frame and the saprating lines corresponding to chromosome borders
$img -> interlaced('true');
$img -> rectangle($left_curb*1.25,$top_curb*1.25,$left_curb*1.25+$frame_width,$top_curb*1.25+$frame_height,$black);









open(FIG,">".$output_fig) or die "can't open figure!\n";
binmode FIG;
print FIG $img -> jpeg;
close(FIG);






sub sum_chr_len()
{
   my @chrlen = @_;
#print "@chrlen[0..$#chrlen]\n";

   my $sumlen = $chrlen[0];
   for(my $i=1; $i<=$#chrlen; $i++)
   {
#    print "input chro len ".$chrlen[$i]."\n";
    $sumlen = $sumlen + $chrlen[$i];
   }
   return $sumlen;
}

# my $bar_num=28;
# my $chr_quan=7;
# my $linespacing=8;


# #print $tel_file;

# # read lens file to get length
# my %hash_len;
# my $maxlens=0;
# open(LEN,$len_file) or die "can't open len file.\n";
# while(<LEN>)
# {
# 	my ($chr,$len)=split("\t",$_);
# 	$hash_len{$chr}=$len;
# 	if($len>$maxlens){$maxlens=$len;}
# }
# close($len_file);
# #
# my $single=($frame_width-2*$left_curb)/$maxlens;#length of each pix
# print 1/$single." length\n";
# my $line_height=($frame_height-2*$top_curb-($bar_num+$chr_quan-2)*$linespacing)/$bar_num;#hight of each color bar



# my $align = GD::Text::Align->new($img, valign => 'center', halign => 'center');
# $align->set_font('Arial.ttf',34);
# $img -> interlaced('true');

# # draw chromosome
# foreach my $key (sort keys %hash_len)
# {
# 	for (my $i = 0; $i < $bar_num/$chr_quan; $i++) {
# 		$img -> filledRectangle($left_curb,$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing)+$i*($line_height+$linespacing),$hash_len{$key}*$single+$left_curb,$line_height+$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing)+$i*($line_height+$linespacing),$dodgerblue);
# 	}
# }
# # draw tick
# for(my $i=0;$i<int($maxlens/10000000)+2;$i++)
# {
# 	$img -> line($left_curb+$i*10000000*$single,$top_curb-20,$left_curb+$i*10000000*$single,$top_curb-8,$black);
# }

# my %hash_staarr;
# my %hash_endarr;
# my %hash_merge;
# my %hash_cover;

# # read tel file and put the location into hash_staarr and hash_endarr
# open(TF,$tel_file) or die "can't open telomeric file";
# while(<TF>)
# {
# 	my @linearr=split(/\s+/,$_);
# 	my $chr=$linearr[0];
# 	my $start=$linearr[6];
# 	my $end=$linearr[7];
# 	if($start>$end){$start=$linearr[7];$end=$linearr[6];}
# 	push(@{$hash_staarr{$chr}},$start);
# 	push(@{$hash_endarr{$chr}},$end);
# }
# close($tel_file);
# merge_hash();
# undef %hash_staarr; undef %hash_endarr;
# #print("@{$hash_merge{2}}\n");
# drawblock(1,1);

# #read cen file
# open(CF,$cen_file) or die "can't open centromeric file.";
# while (<CF>) 
# {
# 	my @linearr=split(/\s+/,$_);
# 	my $chr=$linearr[1];
# 	my $start=$linearr[8];
# 	my $end=$linearr[9];
# 	if($start>$end){$start=$linearr[9];$end=$linearr[8];}
# 	push(@{$hash_staarr{$chr}},$start);
# 	push(@{$hash_endarr{$chr}},$end);
# }
# close($cen_file);
# merge_hash();
# undef %hash_staarr; undef %hash_endarr;
# drawblock(2,2);
# #read ltr file
# open(LTR,$ltr_file) or die "can't open ltr file.";
# while (<LTR>) 
# {
# 	my @linearr=split(/\s+/,$_);
# 	my $chr=$linearr[0];
# 	my $start=$linearr[1];
# 	my $end=$linearr[2];
# 	if($start>$end){$start=$linearr[2];$end=$linearr[1];}
# 	push(@{$hash_staarr{$chr}},$start);
# 	push(@{$hash_endarr{$chr}},$end);
# }
# close($ltr_file);
# merge_hash();
# undef %hash_staarr; undef %hash_endarr;
# drawblock(3,3);
# #read cds file
# open(CDS,$cds_file) or die "can't open cds file.";
# while (<CDS>) 
# {
# 	my @linearr=split(/\s+/,$_);
# 	my $chr=$linearr[0];
# 	my $start=$linearr[1];
# 	my $end=$linearr[2];
# 	if($start>$end){$start=$linearr[2];$end=$linearr[1];}
# 	push(@{$hash_staarr{$chr}},$start);
# 	push(@{$hash_endarr{$chr}},$end);
# }
# close($cds_file);
# merge_hash();
# undef %hash_staarr; undef %hash_endarr;
# drawblock(4,4);

# binmode FIG;
# print FIG $img -> jpeg;
# close(FIG);
# #set color by values
# sub setcolor
# {
# 	my $type=$_[0]; my $value=$_[1];
# 	#print $type." ".$value;
# 	if ($type eq "1") {
# 		if($value<0.000028){return $orange;}
# 		elsif($value<1){return $firebrick;}
# 	}
# 	elsif($type eq "2"){
# 		if($value<0.00009){return $orange;}
# 		elsif($value<0.005){return $firebrick;}
# 	}
# 	elsif($type eq "3"){
# 		if($value<0.005){return $dodgerblue;}
# 		elsif($value<0.01){return $orange;}
# 		elsif($value<1){return $firebrick;}
# 	}
# 	elsif($type eq "4"){
# 		if($value<0.005){return $dodgerblue;}
# 		elsif($value<0.01){return $orange;}
# 		elsif($value<1){return $firebrick;}
# 	}
# }
# sub merge_hash()
# {
# 	foreach my $key (sort keys %hash_staarr)
# 	{
# 		# sort each location array
# 		@{$hash_staarr{$key}}=sort {$a<=>$b} @{$hash_staarr{$key}};
# 		@{$hash_endarr{$key}}=sort {$a<=>$b} @{$hash_endarr{$key}};
# 		# push all braskpoints
# 		my $min=int(@{$hash_staarr{$key}}[0]/$myblock)+1;
# 		my $max=int(@{$hash_endarr{$key}}[-1]/$myblock)+1;
# 		#print $min." ".$max."\n";
# 		for($min;$min<$max;$min++)
# 		{
# 			push(@{$hash_staarr{$key}},$min*$myblock+1.5);
# 		    push(@{$hash_endarr{$key}},$min*$myblock);
# 		}
# 		# re-sort array
# 		@{$hash_staarr{$key}}=sort {$a<=>$b} @{$hash_staarr{$key}};
# 		@{$hash_endarr{$key}}=sort {$a<=>$b} @{$hash_endarr{$key}};
# 		# merge array
# 		my $value_s=@{$hash_staarr{$key}}[0];
# 		for(my $i=0;$i<@{$hash_staarr{$key}};$i++)
# 		{
# 			if($value_s<=@{$hash_endarr{$key}}[$i])
# 			{
# 				if(@{$hash_endarr{$key}}[$i]<@{$hash_staarr{$key}}[$i+1]-1||$i+1==@{$hash_staarr{$key}}||@{$hash_endarr{$key}}[$i]%$myblock==0)
# 				{
# 					push(@{$hash_merge{$key}},$value_s);
# 					push(@{$hash_merge{$key}},@{$hash_endarr{$key}}[$i]);
# 					if(@{$hash_endarr{$key}}[$i]%$myblock==0)
# 					{
# 					    $value_s=@{$hash_endarr{$key}}[$i]+1;next;
# 				    }
# 					$value_s=@{$hash_staarr{$key}}[$i+1];
# 				}
# 			}
# 			else
# 			{
# 				$value_s=@{$hash_staarr{$key}}[$i+1];
# 			}
# 		}
# 	}
# 	#print("@{$hash_merge{2}}\n");#test print
# }
# sub drawblock
# {
# 	my $index=$_[0]; my $type=$_[1];
# 	foreach my $key (sort keys %hash_merge)
# 	{
# 		for (my $i = 0; $i+1 < @{$hash_merge{$key}};$i+=2) {
# 			 $hash_cover{$key}{int(@{$hash_merge{$key}}[$i]/$myblock)+1}+=@{$hash_merge{$key}}[$i+1]-@{$hash_merge{$key}}[$i]+1;		 
# 		}
# 	}
# 	undef %hash_merge;
# 	foreach my $key (sort keys %hash_cover)
# 	{
# 		foreach my $sec_key (sort keys %{$hash_cover{$key}})
# 		{
# 			my $cover=%{$hash_cover{$key}}{$sec_key}/$myblock;
# 			#print $cover." ";
# 			my $color=$dodgerblue;
# 			$color=setcolor($type,$cover);
# 			$img -> filledRectangle($left_curb+($sec_key-1)*$myblock*$single,$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing)+($index-1)*($line_height+$linespacing),$left_curb+$sec_key*$myblock*$single,$line_height+$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing)+($index-1)*($line_height+$linespacing),$color);
# 		}
# 		#$img -> filledRectangle($left_curb,$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing),$hash_len{$key}*$single+$left_curb,$line_height+$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing),"blue");
# 	}
# 	undef %hash_cover;
# }


