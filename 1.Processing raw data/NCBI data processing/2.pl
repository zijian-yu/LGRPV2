use strict;
open (IN,"$ARGV[0].1.gff") or die "cannot open file1 due to $!\n";
open (OUT,">$ARGV[0].2.gff")or die "cannot open file2 due to $!\n";
my $old;#针对scafford的gff进行过滤和基因的重命名；最后生成文件在Excel里按scafford名，5’，3’端进行升序排序。
my $num;
while (<IN>)
{chomp;
	my @arr=split(/\s+/,$_);
	if ($old eq $arr[0]){$num++}
	else{$num=1}	
print OUT "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[0]"."g";
	printf OUT ("%0*d\n",4,$num."\n");
$old= $arr[0];
}
close IN;
close OUT;
`perl 3.pl $ARGV[0]`;
