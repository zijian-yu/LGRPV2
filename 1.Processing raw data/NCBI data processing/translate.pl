use strict;
#use Bio::AlignIO;
use Bio::Seq;
use Bio::SeqIO;
#use Bio::Seq::EncodedSeq;
my $is_prot_aln = Bio::SeqIO -> new(-file=>"$ARGV[0].cds", -format=>"fasta");
my $outfile="$ARGV[0].pep";
#my $is_prot=$is_prot_aln->next_seq;

my $out=Bio::SeqIO->new(-file =>">".$outfile,-format =>'fasta');
   while (my $prot = $is_prot_aln -> next_seq)
   {
	$out->write_seq($prot->translate());
   }
