use strict;
open (IN,$ARGV[0]) or die "cannot open file1 due to $!\n";
open (OUT,">$ARGV[1]")or die "cannot open file2 due to $!\n";

my %hash;
while (<IN>)
{
	chomp;
	if($_=~/^>/)
	{
		my @array=split(/\s+/,$_);
		my @arr=split(/\=/,$array[1]);
		$arr[1]=~s/\]//g;
		#if($arr[1]=''){print OUT ">000\n";next;}
		if($hash{$arr[1]}){print OUT ">000\n";next;}
		$hash{$arr[1]}=1;
		print OUT ">$arr[1]\n";

	}
	else{print OUT "$_\n";}

}
