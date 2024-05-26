use strict;
open (IN,"$ARGV[0].1.gff") or die "cannot open file1 due to $!\n";
open (OUT,">$ARGV[0].2.gff")or die "cannot open file2 due to $!\n";
my $old;#Õë¶ÔscaffordµÄgff½øÐÐ¹ýÂËºÍ»ùÒòµÄÖØÃüÃû£»×îºóÉú³ÉÎÄ¼þÔÚExcelÀï°´scaffordÃû£¬5¡¯£¬3¡¯¶Ë½øÐÐÉýÐòÅÅÐò¡£
my $num;
while (<IN>)
{chomp;
	my @arr=split(/\s+/,$_);
	if ($old eq $arr[0]){$num++}
	else{$num=1}	
print OUT "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[0]"."g";
	printf OUT ("%0*d\n",5,$num."\n");
$old= $arr[0];
}
close IN;
close OUT;
`perl 3.pl $ARGV[0]`;
