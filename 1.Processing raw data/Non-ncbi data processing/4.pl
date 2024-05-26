use strict;
my $gff3file=$ARGV[0];###�ĺõ�λ���ļ�gff������׺����

open(IN,$gff3file) or die "cannot open file due to file $gff3file $!\n";
my $profile=$ARGV[1];###pep�ļ���������׺����
open (PRO,$profile) or die "cannot open file due to file $profile $!\n";
my $outfile=$ARGV[2];###������ļ���������׺����
open (OUT,">",$outfile) or die "cannot open file due to file $outfile $!\n";
my %hash;
while(<IN>)
{	$_=~s/[\n\r]//g;
	my @array=split(/\t/,$_);


	$hash{$array[4]}=$array[5];
}

while(<PRO>)
{
	s/[\n\r]//g;
	if($_ !~ /^>/){print OUT $_."\n"; next;}
	$_=~s/>//g;
	my @arra=split(/\s+/,$_);
	$arra[0]=~s/\.m01\.CDS//g;
	#my $seqid=$hash{$arra[4]};
	if(exists $hash{$arra[0]})

	{print OUT ">".$hash{$arra[0]}."\n";}
	if(not exists $hash{$arra[0]}){print OUT ">"."$arra[0]"."\n";}
}
