use strict;
#use warnings;

my @file = glob "*.fasta";
for (my $i=0;$i<=$#file;$i++)
{
	open(IN, "$file[$i]") or die "cannot open the infile $!";
	my @name=split(/\./,$file[$i]);
    open (OUT, ">$name[0].$name[1]") or die;

    print OUT ">";
	$/ = ">";
	while(<IN>)
	{
	    #chomp;
	    if($_!~/$ARGV[0]/){next;}
	    else{print OUT "$_";}
	}
}