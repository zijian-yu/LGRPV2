use strict;
open (IN,$ARGV[0]) or die "cannot open file1 due to $!\n";
open (OUT,">$ARGV[1].1.gff")or die "cannot open file2 due to $!\n";

my $gene;
my %hash;
my $num=0;
my $chr;
while (<IN>)
{
	chomp;
	if($_=~/^#/){next;}
	if($_=~/mRNA-1:/){next;}

	my @arr=split(/\s+/,$_);
	if($arr[2]=~/mRNA/){
		my @a=split(/;/,$arr[8]);
		$a[0]=~s/ID=//g;
		if($hash{$a[0]}){next;}
		$hash{$a[0]}=1;

		#$arr[0]=~s/\.2//g;		 		#删除原始gff第一列的  .1
		$arr[0]=~s/HiC_scaffold_//g;   			#删除原始gff第一列多余的字母

		#if ($arr[0] ne $chr){$num++;$chr=$arr[0];}


		if($arr[0]=~/^Contig/){next;} 		#删除多余的gff数据

		print OUT "$arr[0]\t$ARGV[1]$arr[0]\t$arr[3]\t$arr[4]\t$arr[6]\t$a[0]\n";
				#目标物种拉丁文缩写
	}
}
