############################################
# ok, this is a perl script to draw some color bar 
# to show Hv's telomeric repeat and centromeric repeat
# add ltr and cds as pos reference
#########
# first demand can be found in the QQ's email at Apr 2, 2016
#########
# There are 4 tpyes of file, total 16 input files
# Will draw 7x4 color bar
############################################


use strict;
use Getopt::Long;
use GD;
use GD::Text::Align;

### total parameter below
#set variate
my $tel_file; #telomeric file
my $cen_file; #centromeric file
my $ltr_file; #ltr file
my $cds_file; #cds file
my $len_file; #lens file
my $myblock=1000000;
#get variate
Getopt::Long::GetOptions(
	'tf=s' => \$tel_file,
	'cf=s' => \$cen_file,
	'l=s' => \$ltr_file,
	'c=s' => \$cds_file,
	'len=s' => \$len_file,
	'b=i' => \$myblock);
open(FIG,">draw.hv.tel.cen.jpg") or die "can't open figure!\n";
my $bar_num=28;
my $chr_quan=7;
my $linespacing=8;

#########################
#print $tel_file;

# read lens file to get length
my %hash_len;
my $maxlens=0;
open(LEN,$len_file) or die "can't open len file.\n";
while(<LEN>)
{
	my ($chr,$len)=split("\t",$_);
	$hash_len{$chr}=$len;
	if($len>$maxlens){$maxlens=$len;}
}
close($len_file);
# set parameter of figure
my $frame_width=2400;
my $frame_height=2400;
my $left_curb=200;
my $top_curb=200;
my $single=($frame_width-2*$left_curb)/$maxlens;#length of each pix
print 1/$single." length\n";
my $line_height=($frame_height-2*$top_curb-($bar_num+$chr_quan-2)*$linespacing)/$bar_num;#hight of each color bar

my $img=GD::Image -> new($frame_width,$frame_height);

my $white= $img->colorAllocate(255,255,255);
my $black= $img->colorAllocate(0,0,0);
my $mintcream= $img->colorAllocate(245,255,250);
my $dodgerblue= $img->colorAllocate(30,144,255);
my $firebrick= $img->colorAllocate(255,0,0);
my $darkviolet= $img->colorAllocate(148,0,211);
my $orange= $img->colorAllocate(255,165,0);
my $white= $img->colorAllocate(255,255,255);
my $white= $img->colorAllocate(255,255,255);

my $align = GD::Text::Align->new($img, valign => 'center', halign => 'center');
$align->set_font('Arial.ttf',34);
$img -> interlaced('true');

# draw chromosome
foreach my $key (sort keys %hash_len)
{
	for (my $i = 0; $i < $bar_num/$chr_quan; $i++) {
		$img -> filledRectangle($left_curb,$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing)+$i*($line_height+$linespacing),$hash_len{$key}*$single+$left_curb,$line_height+$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing)+$i*($line_height+$linespacing),$dodgerblue);
	}
}
# draw tick
for(my $i=0;$i<int($maxlens/10000000)+2;$i++)
{
	$img -> line($left_curb+$i*10000000*$single,$top_curb-20,$left_curb+$i*10000000*$single,$top_curb-8,$black);
}

my %hash_staarr;
my %hash_endarr;
my %hash_merge;
my %hash_cover;

# read tel file and put the location into hash_staarr and hash_endarr
open(TF,$tel_file) or die "can't open telomeric file";
while(<TF>)
{
	my @linearr=split(/\s+/,$_);
	my $chr=$linearr[0];
	my $start=$linearr[6];
	my $end=$linearr[7];
	if($start>$end){$start=$linearr[7];$end=$linearr[6];}
	push(@{$hash_staarr{$chr}},$start);
	push(@{$hash_endarr{$chr}},$end);
}
close($tel_file);
merge_hash();
undef %hash_staarr; undef %hash_endarr;
#print("@{$hash_merge{2}}\n");
drawblock(1,1);

#read cen file
open(CF,$cen_file) or die "can't open centromeric file.";
while (<CF>) 
{
	my @linearr=split(/\s+/,$_);
	my $chr=$linearr[1];
	my $start=$linearr[8];
	my $end=$linearr[9];
	if($start>$end){$start=$linearr[9];$end=$linearr[8];}
	push(@{$hash_staarr{$chr}},$start);
	push(@{$hash_endarr{$chr}},$end);
}
close($cen_file);
merge_hash();
undef %hash_staarr; undef %hash_endarr;
drawblock(2,2);
#read ltr file
open(LTR,$ltr_file) or die "can't open ltr file.";
while (<LTR>) 
{
	my @linearr=split(/\s+/,$_);
	my $chr=$linearr[0];
	my $start=$linearr[1];
	my $end=$linearr[2];
	if($start>$end){$start=$linearr[2];$end=$linearr[1];}
	push(@{$hash_staarr{$chr}},$start);
	push(@{$hash_endarr{$chr}},$end);
}
close($ltr_file);
merge_hash();
undef %hash_staarr; undef %hash_endarr;
drawblock(3,3);
#read cds file
open(CDS,$cds_file) or die "can't open cds file.";
while (<CDS>) 
{
	my @linearr=split(/\s+/,$_);
	my $chr=$linearr[0];
	my $start=$linearr[1];
	my $end=$linearr[2];
	if($start>$end){$start=$linearr[2];$end=$linearr[1];}
	push(@{$hash_staarr{$chr}},$start);
	push(@{$hash_endarr{$chr}},$end);
}
close($cds_file);
merge_hash();
undef %hash_staarr; undef %hash_endarr;
drawblock(4,4);

binmode FIG;
print FIG $img -> jpeg;
close(FIG);
#set color by values
sub setcolor
{
	my $type=$_[0]; my $value=$_[1];
	#print $type." ".$value;
	if ($type eq "1") {
		if($value<0.000028){return $orange;}
		elsif($value<1){return $firebrick;}
	}
	elsif($type eq "2"){
		if($value<0.00009){return $orange;}
		elsif($value<0.005){return $firebrick;}
	}
	elsif($type eq "3"){
		if($value<0.005){return $dodgerblue;}
		elsif($value<0.01){return $orange;}
		elsif($value<1){return $firebrick;}
	}
	elsif($type eq "4"){
		if($value<0.005){return $dodgerblue;}
		elsif($value<0.01){return $orange;}
		elsif($value<1){return $firebrick;}
	}
}
sub merge_hash()
{
	foreach my $key (sort keys %hash_staarr)
	{
		# sort each location array
		@{$hash_staarr{$key}}=sort {$a<=>$b} @{$hash_staarr{$key}};
		@{$hash_endarr{$key}}=sort {$a<=>$b} @{$hash_endarr{$key}};
		# push all braskpoints
		my $min=int(@{$hash_staarr{$key}}[0]/$myblock)+1;
		my $max=int(@{$hash_endarr{$key}}[-1]/$myblock)+1;
		#print $min." ".$max."\n";
		for($min;$min<$max;$min++)
		{
			push(@{$hash_staarr{$key}},$min*$myblock+1.5);
		    push(@{$hash_endarr{$key}},$min*$myblock);
		}
		# re-sort array
		@{$hash_staarr{$key}}=sort {$a<=>$b} @{$hash_staarr{$key}};
		@{$hash_endarr{$key}}=sort {$a<=>$b} @{$hash_endarr{$key}};
		# merge array
		my $value_s=@{$hash_staarr{$key}}[0];
		for(my $i=0;$i<@{$hash_staarr{$key}};$i++)
		{
			if($value_s<=@{$hash_endarr{$key}}[$i])
			{
				if(@{$hash_endarr{$key}}[$i]<@{$hash_staarr{$key}}[$i+1]-1||$i+1==@{$hash_staarr{$key}}||@{$hash_endarr{$key}}[$i]%$myblock==0)
				{
					push(@{$hash_merge{$key}},$value_s);
					push(@{$hash_merge{$key}},@{$hash_endarr{$key}}[$i]);
					if(@{$hash_endarr{$key}}[$i]%$myblock==0)
					{
					    $value_s=@{$hash_endarr{$key}}[$i]+1;next;
				    }
					$value_s=@{$hash_staarr{$key}}[$i+1];
				}
			}
			else
			{
				$value_s=@{$hash_staarr{$key}}[$i+1];
			}
		}
	}
	#print("@{$hash_merge{2}}\n");#test print
}
sub drawblock
{
	my $index=$_[0]; my $type=$_[1];
	foreach my $key (sort keys %hash_merge)
	{
		for (my $i = 0; $i+1 < @{$hash_merge{$key}};$i+=2) {
			 $hash_cover{$key}{int(@{$hash_merge{$key}}[$i]/$myblock)+1}+=@{$hash_merge{$key}}[$i+1]-@{$hash_merge{$key}}[$i]+1;		 
		}
	}
	undef %hash_merge;
	foreach my $key (sort keys %hash_cover)
	{
		foreach my $sec_key (sort keys %{$hash_cover{$key}})
		{
			my $cover=%{$hash_cover{$key}}{$sec_key}/$myblock;
			#print $cover." ";
			my $color=$dodgerblue;
			$color=setcolor($type,$cover);
			$img -> filledRectangle($left_curb+($sec_key-1)*$myblock*$single,$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing)+($index-1)*($line_height+$linespacing),$left_curb+$sec_key*$myblock*$single,$line_height+$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing)+($index-1)*($line_height+$linespacing),$color);
		}
		#$img -> filledRectangle($left_curb,$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing),$hash_len{$key}*$single+$left_curb,$line_height+$top_curb+($key-1)*(4*($line_height+$linespacing)+$linespacing),"blue");
	}
	undef %hash_cover;
}


