
#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::SeqIO;
use Getopt::Long;
use Bio::Seq;



my $file;	


GetOptions('f=s'=>\$file);



my $outfile = $file ;
$outfile =~ s/\.\w+$// ;
$outfile .= 'processed.fa' ;

create_fasta_file($file, $outfile);


sub create_fasta_file {

	my ($file, $outfile) = @_;
	my $count = 0;

	my $in  = Bio::SeqIO->new(-file => $file, -format => "fasta");
	my $protein = {};
	while(my $seqs = $in->next_seq) {
	
		my $id  = $seqs->display_id;	
		my $seq = $seqs->seq;
		my $desc = $seqs->desc;

		my ($start) = $desc =~ /\((\d+)-/;
		my ($stop) = $desc =~ /-(\d+)\)$/;
        my $id_2 = (split '\|',$id)[1];

        $protein->{$seq}->{id} = $id;
        $protein->{$seq}->{id_2} = $id_2."_".$start."-".$stop;
        $protein->{$seq}->{desc} = $desc;
    }

    # write predicted
    my $count_pred = 0;
	my $out_pred  = Bio::SeqIO->new(-file => ">".$outfile."_predicted.fasta", -format => "fasta");
    foreach my $seq (keys %$protein) {
        my $id = $protein->{$seq}->{id_2};
        my $desc = $protein->{$seq}->{desc};
        my $seq_obj = Bio::Seq->new(-seq => $seq,-display_id => $id, -desc => $desc);
        $out_pred->write_seq($seq_obj);
        $count_pred++;
    }
	print "Total unique peptides $count_pred\n";

}


