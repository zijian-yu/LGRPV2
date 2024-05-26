use strict;
my $in="$ARGV[0].2.gff";#排过序的gff文件（excel打开按物种名染色体列和5`末端升序）
open(IN,$in)or die "can not open my file $in due to $!.\n";
my $out="$ARGV[0].new.gff";#newgff文件
open (OUT,">",$out) or  die "can not open my file $out due to $!.\n";
my $out2="$ARGV[0].lens";#lens文件
open (OUT2,">",$out2) or  die "can not open my file $out2 due to $!.\n";
my $num;
my $old;
my $len;
my @a;
my @b;
my $last;
while (<IN>)
{
$_=~s/[\n\r]//g;
my @array=split(/\t/,$_);


$last=$array[0];
if ($old eq $array[0])
{print OUT $_."\t".$num."\n";}

else {
	push @a, $num-1;
	my $len=$old=~ s/\D+//g;
	push @b,$old;
	$num=1;print OUT $_."\t".$num."\n";
}
$num++;
$old=$array[0];
}

#print "$len\t$old\n";
#$old=~s/\D+//g;
push @b, $old;
push @a, $num-1;
for (my $i=1;$i<=$#a;$i++)
{if (1)
	{$b[$i]=~s/\D+//;
	print OUT2 "$b[$i]\t$a[$i]\n";}
}	
