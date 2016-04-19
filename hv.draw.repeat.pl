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
# hv1 Hv1g00001 4070  4229
# 4G  Hv1g00002 4311  4498
# 07g Hv1g00003 4596  4787
##### M8:
# the same as BLAST M8 result
##### LTR:
# chr\tstart_pos\tend_pos
# hv1 246050  246207
# Hv05g 246703  247050
# 1G  247051  247225
##########
# in this script, user can select special chromosome, filter Blast M8 file only for dotplot 
# can set length of each speci-grid-line, can set length of block, can set 2 boundary values for each color bar
# here, join the chromosome number with "_", and also join the boundary values with "_".
# I also put some parameters ahead the script, for example colors array, input file names. It can be changed conveniently.

# It is a example input command below:
# perl hv.draw.repeat.pl -e 1e-5 -s 200 -n 10 -gs 10000000 -gq 10000000 -cs 1_3_6 -cq 2_5 -b 1000000 -vc 0.4_0.8 -vt 0.05_0.1 -vl 0.0005_0.003 -vg 0.005_0.01
####### parameters:
# -e  e-value of BLAST  default:5e-2
# -s  score of BLAST  default:0
# -n  hitnumber of BLAST  default:50
# -gs grid line length of subject,if none,input "0" default:10M
# -gq grid line lenght of query,if none,input "0" default:10M
# -b  length of block default:1M
# -cs selected subject chrmosome number default:all
# -cq selected query chromosome number  default:all
# -vc boundary value of centromeric default:0.04_0.08
# -vt boundary value of telomeric default:0.05_0.1
# -vl boundary value of LTR default:0.005_0.05
# -vg boundary value of cds default:0.005_0.05
############################################

use strict;
use Getopt::Long;
use GD;
use GD::Text::Align;
# set parameter of figure
my $frame_width=2000;
my $frame_height=2000;
my $left_curb=300;
my $top_curb=300;
my $img=GD::Image -> new($frame_width+1.5*$left_curb,$frame_height+1.5*$top_curb);


########################### total parameters below
##set variate
#files
my $tel_sfile="telomeric.hv.1e-3.blastn.m8"; #telomeric file
my $cen_sfile="hv.genome.hv.centromeric.blastn.m8"; #centromeric file
my $ltr_sfile="hv.ltr.txt"; #ltr file
my $gff_sfile="os.chr.gff"; #gff file
my $tel_qfile="os.telomeric.telomere.associated.1e-3.blastn.m8"; #telomeric file
my $cen_qfile="os.genome.os.centro.1e-3.blastn.m8"; #centromeric file
my $ltr_qfile="os.ltr.txt"; #ltr file
my $gff_qfile="os.chr.gff"; #gff file
my $plot_file="os.hv.cds.blastn.m8.new"; #blast file for dotplot
my $output_fig="hv.draw.repeat.jpg"; #output figure
#set color
my $white= $img->colorAllocate(255,255,255);
my $black= $img->colorAllocate(0,0,0);
my $gray= $img->colorAllocate(100,100,100);
my $red= $img->colorAllocate(255,0,0);
my $green= $img->colorAllocate(0,255,0);
my $blue= $img->colorAllocate(0,0,255);
my $mintcream= $img->colorAllocate(245,255,250);
my $dodgerblue= $img->colorAllocate(30,144,255);
my $darkviolet= $img->colorAllocate(148,0,211);
my $orange= $img->colorAllocate(255,165,0);
#font
my $align = GD::Text::Align->new($img, valign => 'center', halign => 'center', color => $black);
$align->set_font('Arial.ttf',34);
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
#get values array of hot map color standar
my @vcarr=split('_',$val_cen);
my @vtarr=split('_',$val_tel);
my @vlarr=split('_',$val_ltr);
my @vgarr=split('_',$val_cds);
#set color parameters
my $bar_qua=4;
my $linespacing=8;
my $lineheight=($top_curb-($bar_qua-1)*$linespacing)/$bar_qua;
my $columewidth=($left_curb-($bar_qua-1)*$linespacing)/$bar_qua;
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
while (<QGFF>){
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
foreach my $key (sort {$tempsort{$a} <=> $tempsort{$b}} keys %tempsort){
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
while (<SGFF>){
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
foreach my $key (sort {$tempsort{$a} <=> $tempsort{$b}} keys %tempsort){
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
if ($chrn_q eq "all"){
  my @a=();
  foreach my $key (sort{$a<=>$b} keys %querychr2plen){push(@a,$key);}
  $chrn_q=join("_",@a);
}
my @querychr=split('_',$chrn_q);

my @sbjctchrlen = ();
my %sbjctchr2order;
#if selected chr is all chr, change chrn_ to all chr number
if ($chrn_s eq "all"){
  my @a=();
  foreach my $key (sort{$a<=>$b} keys %sbjctchr2plen){push(@a,$key);}
  $chrn_s=join("_",@a);
}
my @sbjctchr=split('_',$chrn_s);
if($poro eq "p"){
  for(my $i=0; $i<=$#querychr; $i++){
    $querychr[$i]=~s/[^0-9]//g;
    my $chr = int($querychr[$i]);
    $querychrlen[$#querychrlen+1] = $querychr2plen{$chr}; #sum length
    $querychr2order{$chr} = $#querychrlen; #total order
  }

  for(my $i=0; $i<=$#sbjctchr; $i++){
    $sbjctchr[$i]=~s/[^0-9]//g;
    my $chr = int($sbjctchr[$i]);
    $sbjctchrlen[$#sbjctchrlen+1] = $sbjctchr2plen{$chr}; #sum length
    $sbjctchr2order{$chr} = $#sbjctchrlen; #total order
  }
}
else{
  for(my $i=0; $i<=$#querychr; $i++){
    $querychr[$i]=~s/[^0-9]//g;
    my $chr = int($querychr[$i]);
    $querychrlen[$#querychrlen+1] = $querychr2olen{$chr}; #sum length
    $querychr2order{$chr} = $#querychrlen; #total order
  }

  for(my $i=0; $i<=$#sbjctchr; $i++){
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
# draw the frame of dotplot part and the saprating lines corresponding to chromosome borders
$img -> interlaced('true');
#frame
$img -> rectangle($left_curb*1.25,$top_curb*1.25,$left_curb*1.25+$frame_width,$top_curb*1.25+$frame_height,$black);
$img -> rectangle($left_curb*1.25-1,$top_curb*1.25-1,$left_curb*1.25+$frame_width+1,$top_curb*1.25+$frame_height+1,$black);
#lines
my @query_chro_pos = ();
my @sbjct_chro_pos = ();

my $accumulated_length = 0;
for(my $i=0; $i<=$#querychrlen; $i++){
   my $posy1 = $top_curb*1.25;   
   my $posy2 = $top_curb*1.25 + $frame_height;
#draw grid line
   if($grid_q != 0){     
     for(my $sum_len1 = $grid_q;$sum_len1<$querychrlen[$i];$sum_len1+=$grid_q){
      my $posx0 = $left_curb*1.25 + int(($accumulated_length+$sum_len1)/$scale_ratio1);
      $img -> line($posx0, $posy1, $posx0, $posy2, $green);
     }
   }
   #draw lines between chromosome
   $accumulated_length += $querychrlen[$i];
   my $length = int($querychrlen[$i]/$scale_ratio1);
   my $posx1 = $left_curb*1.25 + int($accumulated_length/$scale_ratio1);
   $query_chro_pos[$i] = $posx1;
   my $posx2 = $posx1;
   $img -> line($posx1, $posy1, $posx2, $posy2, $black);
   $img -> line($posx1+1, $posy1, $posx2+1, $posy2, $black);
   # plot chromosome number
   my $chr = $querychr[$i];
   $align->set_text($chr);
   $align->draw($posx1-int($length/2),50,0);   
}
#sbjct lines and grid lines
$accumulated_length = 0;
for(my $i=0; $i<=$#sbjctchrlen; $i++){
   my $posx1 = $left_curb*1.25;
   my $posx2 = $left_curb*1.25+$frame_width;
   #grid lines
   if($grid_s != 0){     
     for(my $sum_len2 = $grid_s;$sum_len2<$sbjctchrlen[$i];$sum_len2+=$grid_s){
      my $posy0 = $top_curb*1.25 + int(($accumulated_length+$sum_len2)/$scale_ratio2);
      $img -> line($posx1, $posy0, $posx2, $posy0, $green);
     }
   }
   #lines between chromosomes
   $accumulated_length += $sbjctchrlen[$i];
   my $length = int($sbjctchrlen[$i]/$scale_ratio2);
   my $posy1 = $top_curb*1.25 + int($accumulated_length/$scale_ratio2);
   $sbjct_chro_pos[$i] = $posy1;
   my $posy2 = $posy1;
   $img -> line($posx1, $posy1, $posx2, $posy2, $black);
   $img -> line($posx1, $posy1+1, $posx2, $posy2+1, $black);
   #chromosome numbers
   my $chr = $sbjctchr[$i];
   $align->set_text($chr);
   $align->draw(50, $posy1-int($length/2), 1.57);
}
###########
###########start draw dot
open(DOT,$plot_file) or die "can't open blast file for dotplot.\n";
my $lastquery = "";
my $hitnum = 0;
while(<DOT>){
   $_ =~ s/[\n\r]//g;
   my ($querygene, $sbjctgene, $identity, $matchlen, $mismatchnum, $gaplen, $querystart, $queryend, $sbjctstart, $sbjctend, $evalue, $score) = split(/\t/, $_);
   my $querychr = $querygene2chr{$querygene};
   my $sbjctchr = $sbjctgene2chr{$sbjctgene};
# only selected chromosome;
   my $is2skip1 = 1;
   for(my $i=0; $i<=$#querychr; $i++){if($querychr eq $querychr[$i]){$is2skip1 = 0; last;}}
   my $is2skip2 = 1;
   for(my $i=0; $i<=$#sbjctchr; $i++){if($sbjctchr eq $sbjctchr[$i]){$is2skip2 = 0; last;}}
   if($is2skip1 eq 1 || $is2skip2 eq 1){next;}
# max hit can not more than HITNUM;
   if($lastquery ne $querygene){$hitnum = 1;$lastquery = $querygene;}
   else{$hitnum ++;}
   if($hitnum > $HITNUM){next;}
# blast score should be more than SCORE
   if($score < $SCORE){next;}
# e-value should be less than EVALUE
   if($evalue > $EVALUE){next;}
# get x, y value of dot
   my ($posx1, $posy1, $posx2, $posy2, $selfhit1x, $selfhit1y, $selfhit2x, $selfhit2y);
   if($ARGV[5] eq "order"){
      if($querychr2order{$querychr} eq 0){$posx1 = $left_curb + $querygene2order{$querygene}/$scale_ratio1;}
      else{$posx1 = $query_chro_pos[$querychr2order{$querychr}-1] + $querygene2order{$querygene}/$scale_ratio1;}
      if($sbjctchr2order{$sbjctchr} eq 0){$posy1 = $top_curb + $sbjctgene2order{$sbjctgene}/$scale_ratio2;}
      else{$posy1 = $sbjct_chro_pos[$sbjctchr2order{$sbjctchr}-1] +  $sbjctgene2order{$sbjctgene}/$scale_ratio2;}
   }
   elsif($ARGV[5] eq "pos"){
      if($querychr2order{$querychr} eq 0){$posx1 = $left_curb + $querygene2pos{$querygene}/$scale_ratio1;}
      else{$posx1 = $query_chro_pos[$querychr2order{$querychr}-1] + $querygene2pos{$querygene}/$scale_ratio1;}
      if($sbjctchr2order{$sbjctchr} eq 0){$posy1 = $top_curb + $sbjctgene2pos{$sbjctgene}/$scale_ratio2;}
      else{$posy1 = $sbjct_chro_pos[$sbjctchr2order{$sbjctchr}-1] +  $sbjctgene2pos{$sbjctgene}/$scale_ratio2;}
    }  
# set dot's color
   my $color = $gray;
   if($hitnum == 1){$color = $red;}
   elsif($hitnum == 2){$color = $blue;}
# draw dot
   $img -> filledRectangle($posx1, $posy1, $posx1+1, $posy1+1, $color);
}
close($plot_file);
############################################################################
#########Bravo, now, finished all dotplot part drawing!
############################################################################
#########From here, start to add some others into figure
####add color bar abut two sides of dotplot

#draw total selected chromosome background
for(my $i=0;$i<$bar_qua;$i++){
  #draw top
  $img -> filledRectangle($left_curb*1.25,$top_curb*0.25+$i*($lineheight+$linespacing),$left_curb*1.25+$frame_width,$top_curb*0.25+$i*($lineheight+$linespacing)+$lineheight,$dodgerblue);
  #draw left
  $img -> filledRectangle($left_curb*0.25+$i*($columewidth+$linespacing),$top_curb*1.25,$left_curb*0.25+$i*($columewidth+$linespacing)+$frame_height,$top_curb*1.25+$columewidth,$dodgerblue);
}

my %hash_staarr;
my %hash_endarr;
my %hash_merge;
my %hash_cover;
# # read query tel file and draw hot map
open(QTF,$tel_qfile) or die "can't open query telomeric file";
while(<QTF>){
  my @linearr=split(/\s+/,$_);
  $linearr[0]=~s/[^0-9]//g;
  my $chr=int($linearr[0]);
  my $start=$linearr[6];
  my $end=$linearr[7];
  if($start>$end){$start=$linearr[7];$end=$linearr[6];}
  push(@{$hash_staarr{$chr}},$start);
  push(@{$hash_endarr{$chr}},$end);
}
close($tel_qfile);
#union sets
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#draw map
drawblock("vt",1,"q");

# # read query cen file and draw hot map
open(QCF,$cen_qfile) or die "can't open query centromeric file.\n";
while(<QCF>){
  my @linearr=split(/\s+/,$_);
  $linearr[1]=~s/[^0-9]//g;
  my $chr=int($linearr[1]);
  my $start=$linearr[8];
  my $end=$linearr[9];
  if($start>$end){$start=$linearr[9];$end=$linearr[8];}
  push(@{$hash_staarr{$chr}},$start);
  push(@{$hash_endarr{$chr}},$end);
}
close($cen_qfile);
#union sets
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#draw map
drawblock("vc",2,"q");

# # read query ltr file and draw hot map
open(QLF,$ltr_qfile) or die "can't open query ltr file.\n";
while(<QLF>){
  my @linearr=split(/\s+/,$_);
  $linearr[0]=~s/[^0-9]//g;
  my $chr=int($linearr[0]);
  my $start=$linearr[1];
  my $end=$linearr[2];
  if($start>$end){$start=$linearr[2];$end=$linearr[1];}
  push(@{$hash_staarr{$chr}},$start);
  push(@{$hash_endarr{$chr}},$end);
}
close($ltr_qfile);
#union sets
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#draw map
drawblock("vl",3,"q");

# # read query cds file and draw hot map
open(QGF,$cds_qfile) or die "can't open query cds file.\n";
while(<QGF>){
  my @linearr=split(/\s+/,$_);
  $linearr[0]=~s/[^0-9]//g;
  my $chr=int($linearr[0]);
  my $start=$linearr[2];
  my $end=$linearr[3];
  if($start>$end){$start=$linearr[3];$end=$linearr[2];}
  push(@{$hash_staarr{$chr}},$start);
  push(@{$hash_endarr{$chr}},$end);
}
close($cds_qfile);
#union sets
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#draw map
drawblock("vg",4,"q");


# # read sbjct tel file and draw hot map
open(STF,$tel_sfile) or die "can't open sbjct telomeric file";
while(<STF>){
  my @linearr=split(/\s+/,$_);
  $linearr[0]=~s/[^0-9]//g;
  my $chr=int($linearr[0]);
  my $start=$linearr[6];
  my $end=$linearr[7];
  if($start>$end){$start=$linearr[7];$end=$linearr[6];}
  push(@{$hash_staarr{$chr}},$start);
  push(@{$hash_endarr{$chr}},$end);
}
close($tel_sfile);
#union sets
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#draw map
drawblock("vt",1,"s");

# # read sbjct cen file and draw hot map
open(SCF,$cen_sfile) or die "can't open sbjct centromeric file.\n";
while(<SCF>){
  my @linearr=split(/\s+/,$_);
  $linearr[1]=~s/[^0-9]//g;
  my $chr=int($linearr[1]);
  my $start=$linearr[8];
  my $end=$linearr[9];
  if($start>$end){$start=$linearr[9];$end=$linearr[8];}
  push(@{$hash_staarr{$chr}},$start);
  push(@{$hash_endarr{$chr}},$end);
}
close($cen_sfile);
#union sets
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#draw map
drawblock("vc",2,"s");

# # read sbjct ltr file and draw hot map
open(SLF,$ltr_sfile) or die "can't open sbjct ltr file.\n";
while(<SLF>){
  my @linearr=split(/\s+/,$_);
  $linearr[0]=~s/[^0-9]//g;
  my $chr=int($linearr[0]);
  my $start=$linearr[1];
  my $end=$linearr[2];
  if($start>$end){$start=$linearr[2];$end=$linearr[1];}
  push(@{$hash_staarr{$chr}},$start);
  push(@{$hash_endarr{$chr}},$end);
}
close($ltr_sfile);
#union sets
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#draw map
drawblock("vl",3,"s");

# # read sbjct cds file and draw hot map
open(SGF,$cds_sfile) or die "can't open sbjct cds file.\n";
while(<SGF>){
  my @linearr=split(/\s+/,$_);
  $linearr[0]=~s/[^0-9]//g;
  my $chr=int($linearr[0]);
  my $start=$linearr[2];
  my $end=$linearr[3];
  if($start>$end){$start=$linearr[3];$end=$linearr[2];}
  push(@{$hash_staarr{$chr}},$start);
  push(@{$hash_endarr{$chr}},$end);
}
close($cds_sfile);
#union sets
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#draw map
drawblock("vg",4,"s");


open(FIG,">".$output_fig) or die "can't open figure!\n";
binmode FIG;
print FIG $img -> jpeg;
close(FIG);


sub sum_chr_len(){
   my @chrlen = @_;
   my $sumlen = $chrlen[0];
   for(my $i=1; $i<=$#chrlen; $i++){
     $sumlen = $sumlen + $chrlen[$i];
   }
   return $sumlen;
}
sub merge_hash(){
  foreach my $key (sort keys %hash_staarr)  {
    # sort each location array
    @{$hash_staarr{$key}}=sort {$a<=>$b} @{$hash_staarr{$key}};
    @{$hash_endarr{$key}}=sort {$a<=>$b} @{$hash_endarr{$key}};
    # push all braskpoints
    my $min=int(@{$hash_staarr{$key}}[0]/$myblock)+1;
    my $max=int(@{$hash_endarr{$key}}[-1]/$myblock)+1;
    #print $min." ".$max."\n";
    for($min;$min<$max;$min++)    {
      push(@{$hash_staarr{$key}},$min*$myblock+1.5);
      push(@{$hash_endarr{$key}},$min*$myblock);
    }
    # re-sort array
    @{$hash_staarr{$key}}=sort {$a<=>$b} @{$hash_staarr{$key}};
    @{$hash_endarr{$key}}=sort {$a<=>$b} @{$hash_endarr{$key}};
    # merge array
    my $value_s=@{$hash_staarr{$key}}[0];
    for(my $i=0;$i<@{$hash_staarr{$key}};$i++){
      if($value_s<=@{$hash_endarr{$key}}[$i]){
        if(@{$hash_endarr{$key}}[$i]<@{$hash_staarr{$key}}[$i+1]-1||$i+1==@{$hash_staarr{$key}}||@{$hash_endarr{$key}}[$i]%$myblock==0){
          push(@{$hash_merge{$key}},$value_s);
          push(@{$hash_merge{$key}},@{$hash_endarr{$key}}[$i]);
          if(@{$hash_endarr{$key}}[$i]%$myblock==0){
            $value_s=@{$hash_endarr{$key}}[$i]+1;next;
          }
          $value_s=@{$hash_staarr{$key}}[$i+1];
        }
      }
      else{
        $value_s=@{$hash_staarr{$key}}[$i+1];
      }
    }
  }
}
sub drawblock{
  my $q_or_s=$_[2];#query or sbjct
  my $index=$_[1];#which line
  my $type=$_[0];#what is it
  #get hot values
  foreach my $key (sort keys %hash_merge){
    for (my $i = 0; $i+1 < @{$hash_merge{$key}};$i+=2){
      $hash_cover{$key}{int(@{$hash_merge{$key}}[$i]/$myblock)+1}+=@{$hash_merge{$key}}[$i+1]-@{$hash_merge{$key}}[$i]+1;     
    }
  }
  undef %hash_merge;
  foreach my $key (sort keys %hash_cover){
    foreach my $sec_key (sort keys %{$hash_cover{$key}}){
      my $cover=%{$hash_cover{$key}}{$sec_key}/$myblock;
      my $color=$dodgerblue;
      #get color by hot value
      $color=setcolor($type,$cover);
      #draw hot rectangle by block width
      #query
      if($q_or_s eq "q"){
        my $posx1;
        if($querychr2order{$key} eq 0){$posx1 = $left_curb*1.25 + ($sec_key-1)*$myblock/$scale_ratio1;}
        else{$posx1 = $query_chro_pos[$querychr2order{$key}-1] + ($sec_key-1)*$myblock/$scale_ratio1;}
        my $posx2=$posx1+$myblock/$scale_ratio1;
        $img -> filledRectangle($posx1,$top_curb*0.25+($index-1)*($line_height+$linespacing),$posx2,$top_curb*0.25+($index-1)*($line_height+$linespacing)+$lineheight,$color);
      }
      #sbjct
      elsif($q_or_s eq "s"){
        my $posy1;
        if($sbjctchr2order{$key} eq 0){$posy1 = $top_curb*1.25 + ($sec_key-1)*$myblock/$scale_ratio2;}
        else{$posy1 = $sbjct_chro_pos[$sbjctchr2order{$key}-1] +  ($sec_key-1)*$myblock/$scale_ratio2;}
        my $posy2=$posy1+$myblock/$scale_ratio2;
        $img -> filledRectangle($left_curb*0.25+($index-1)*($columewidth+$linespacing),$posy1,$left_curb*0.25+($index-1)*($columewidth+$linespacing)+$columewidth,$posy2,$color);
      }
    }    
  }
  undef %hash_cover;
}
sub setcolor{
  my $type=$_[0]; my $value=$_[1];
  if ($type eq "vt") {
    if($value<=$vtarr[0]){return $dodgerblue;}
    elsif($value<=$vtarr[1]){return $orange;}
    elsif($value<=$vtarr[2]){return $red;}
  }
  elsif($type eq "vc"){
    if($value<=$vcarr[0]){return $dodgerblue;}
    elsif($value<=$vcarr[1]){return $orange;}
    elsif($value<=$vcarr[2]){return $red;}
  }
  elsif($type eq "vl"){
    if($value<=$vlarr[0]){return $dodgerblue;}
    elsif($value<=$vlarr[1]){return $orange;}
    elsif($value<=$vlarr[2]){return $red;}
  }
  elsif($type eq "vg"){
    if($value<=$vgarr[0]){return $dodgerblue;}
    elsif($value<=$vgarr[1]){return $orange;}
    elsif($value<=$vgarr[2]){return $red;}
  }
}