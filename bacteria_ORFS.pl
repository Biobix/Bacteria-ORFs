#!/usr/bin/perl -w
	
use strict;
use warnings;
use diagnostics;

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;

# perl /data/elvis/RESPIN/bacteria_ORFS.pl -f /data/elvis/REPARATION/1.Genome/fasta/SALT_SL1344_GCA_000210855.2.fa -o SALT_Sl1344 -m 30 -g /data/elvis/REPARATION/1.Genome/gtf/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_14028s.ASM2216v1.34.gtf


#----------------------------------------------------
#			VARIBLES
#----------------------------------------------------

my $MINORF 	= 90;	# Minimum ORF length
my $genome;			# Genome file
my $outprefix;		# Output file prefix
my $codons = "ATG,GTG,TTG,CTG";
my $gtf = "";


	##	GET command line arguments
GetOptions(
	'f=s'=>\$genome,
	'm=i'=>\$MINORF,
	'o=s'=>\$outprefix,
    's=s'=>\$codons,
    'g=s'=>\$gtf
);

my %params = (
	f => $genome,
);


	#check if mandatory parameters are properly initialized
my @invalid = grep uninitialized_param($params{$_}), keys %params;	
die "Not properly initialized: @invalid\n" if @invalid;


# perl ORF_features.pl -g genome.fa -m 90 -o outprefix


#----------------------------------------------------
#			EXECUTION
#----------------------------------------------------

my $startRun = time();

my %translationHash = 
	(GCA => "A", GCG => "A", GCT => "A", GCC => "A", GCN => "A",
     TGC => "C", TGT => "C",
     GAT => "D", GAC => "D",
     GAA => "E", GAG => "E",
     TTT => "F", TTC => "F",
     GGA => "G", GGG => "G", GGC => "G", GGT => "G", GGN => "G",
     CAT => "H", CAC => "H",
     ATA => "I", ATT => "I", ATC => "I",
     AAA => "K", AAG => "K",
     CTA => "L", CTG => "L", CTT => "L", CTC => "L", CTN => "L", TTA => "L", TTG => "L",
     ATG => "M",
     AAT => "N", AAC => "N",
     CCA => "P", CCT => "P", CCG => "P", CCC => "P", CCN => "P",
     CAA => "Q", CAG => "Q",
     CGA => "R", CGG => "R", CGC => "R", CGT => "R", CGN => "R",
     AGA => "R", AGG => "R",
     TCA => "S", TCG => "S", TCC => "S", TCT => "S", TCN => "S",
     AGC => "S", AGT => "S",
     ACA => "T", ACG => "T", ACC => "T", ACT => "T", ACN => "T",
     GTA => "V", GTG => "V", GTC => "V", GTT => "V", GTN => "V",
     TGG => "W",
     TAT => "Y", TAC => "Y",
     TAG => "*", TAA => "*", TGA => "*");
	 

my $start_codons = {};
my @scodons = split /,/, $codons;
foreach my $codon(@scodons) {
    if (length($codon) != 3 or !(exists $translationHash{$codon})) {
        print "Codon '$codon' not a valid start codon\n";
	    print "Start codons must be 3 nucleotides long and contain either of A,C,G or T [example: -cdn ATG,GTG,TTG]\n";
        exit(1);
    }

	$start_codons->{$codon} = 1;
}




my $annotated = annotation($gtf);


my $genomes = read_fasta($genome);

my $ORFs = find_all_ORFs($genomes);

ORF_fasta_file_all($ORFs, $outprefix);
ORF_fasta_file_longest($ORFs, $outprefix);

print "Done...\n";
timer($startRun);	# Get Run time


##################
##	SUBS
##################



sub find_all_ORFs {
	
	my $genome 		= $_[0];
    my $ORFs 		= {};
	my $count_orf	= 0;
	my $max_length 	= 21;

	foreach my $region (keys %$genome) {
		my $directStrand = $genome->{$region};
		my $reverseComplement = revdnacomp($directStrand);

		for (my $i = 0; $i < 2; $i = $i + 1) {
			my $sequenceEntry = "";
			my $strand = "";
			if ($i == 0) {		# forward strand.
				$sequenceEntry = $directStrand;
				$strand = "+";
			} else {			# reverse strand
				$sequenceEntry = $reverseComplement;
				$strand = "-";
			}

			my $starts = {};
			foreach my $start_cdn (keys %$start_codons) {
				my $offset = 0;
				my $start = index($sequenceEntry, $start_cdn, $offset);
				$starts->{$start} = 1;
				while ($start != -1) {
					$offset = $start + 1;
					$start = index($sequenceEntry, $start_cdn, $offset);
					$starts->{$start} = 1;
				}
			}

			# check if there is a corresponding stop to all start positions
			foreach my $start (sort {$a<=>$b} keys %$starts) {

				my $stop;
				my $aa_seq = "";
				my $tr_seq = "";
				my $strt_cdn = uc(substr($sequenceEntry, $start, 3));
				for (my $j = $start; $j <= (length($sequenceEntry) - 3); $j = $j + 3) {
					my $codon = uc(substr($sequenceEntry, $j, 3));
					last if (length($codon) < 3);
					my $aa = $translationHash{$codon};
					if ($aa eq "*") {
						$stop = $j;
						last;
					}
					$aa_seq = $aa_seq.$aa;
					$tr_seq = $tr_seq.$codon;
				}

				if (defined $stop and ($stop - $start) >= $MINORF ) {	
					my $length = $stop - $start;
					$aa_seq =~ s/^./M/;
					if ($strand eq '-') {
						my $start_tmp 	= length($directStrand) - $stop + 1;
						$stop 			= length($directStrand) - $start;
						$start = $start_tmp;
					} else {
						$start += 1;
					}

					my $id = $region.":".$start."-".$stop;
					$ORFs->{$id}->{region} 		= $region;
					$ORFs->{$id}->{start} 		= $start;
					$ORFs->{$id}->{stop} 		= $stop;
					$ORFs->{$id}->{len} 		= $length;
					$ORFs->{$id}->{strand} 		= $strand;
					$ORFs->{$id}->{tr_seq} 		= $tr_seq;
					$ORFs->{$id}->{aa_seq} 		= $aa_seq;
					$ORFs->{$id}->{start_codon} = $strt_cdn;
					$ORFs->{$id}->{gene} 		= ($strand eq "+") ? $region.":FAM+".$stop: $region.":FAM-".$start;;

					if ($max_length < $length) {$max_length = $length;}
					$count_orf++;
				}
			}
		}
	}
	print "Total number of ORFs $count_orf maximum ORF length $max_length\n";
	return $ORFs;
}




sub ORF_fasta_file_all {

	my $ORFs = $_[0];
	my $file = $_[1];

	open(F, ">", $file."_6FT_all.fasta") or die ("Can't creat file: orf_length_dist.txt: $!\n");
	foreach my $id (sort keys %$ORFs) {
		my $seq = $ORFs->{$id}->{aa_seq};
        my $gene =  $ORFs->{$id}->{gene};
        my $anno = (exists $annotated->{$gene}->{$id}->{name}) ? $annotated : "";
		print F ">generic|$id|strand:$ORFs->{$id}->{strand} start_codon:$ORFs->{$id}->{start_codon} length:$ORFs->{$id}->{len} $anno\n$seq\n";
	}
	close F;
}


sub ORF_fasta_file_longest {

	my $ORFs = $_[0];
	my $file = $_[1];

	my $ORFs_longest = {};
	foreach my $ORF (keys %$ORFs) {
		my $gene = $ORFs->{$ORF}->{gene};
		if (exists $ORFs_longest->{$gene}) {
			if ($ORFs->{$ORF}->{len} > $ORFs->{$ORFs_longest->{$gene}}->{len}) {
				$ORFs_longest->{$gene} = $ORF;
			}
		} else {
			$ORFs_longest->{$gene} = $ORF;
		}
	}
	
	open(F, ">", $file."_6FT_longest.fasta") or die ("Can't creat file: orf_length_dist.txt: $!\n");
	foreach my $gene (sort keys %$ORFs_longest) {
		my $id = $ORFs_longest->{$gene};
        my $anno = (exists $annotated->{$gene}->{$id}->{name}) ? $annotated : "";
		print F ">generic|$id|strand:$ORFs->{$id}->{strand} start_codon:$ORFs->{$id}->{start_codon} length:$ORFs->{$id}->{len} $anno\n$ORFs->{$id}->{aa_seq}\n";
	}
	close F;
}



sub revdnacomp {
    my $dna = shift;
    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}


sub read_fasta {
	my $file = $_[0];
	my $cdna = {};

	my $in  = Bio::SeqIO->new(-file => $file, -format => "fasta");
	while(my $seqs = $in->next_seq) {
		my $id  = $seqs->display_id;	
		my $seq = $seqs->seq;
		my $desc = $seqs->desc;
		$cdna->{$id} = $seq;
	}
	return $cdna;
}



sub annotation {

	my $file = $_[0];
	my $annotated = {};

	open (F, $file) or die  "Error reading file: $file";
	while (<F>) {
		next if (/^#/);
		chomp $_;

		my @line = (split '\t', $_);
		my $region = $line[0];
		my $feature = $line[2];
		my $start = $line[3];
		my $stop = $line[4];
		my $strand = $line[6];

		if ($feature eq 'CDS') {

			my ($biotype) = $_ =~ /gene_biotype."?([^";]+)"?/;
			my ($gene_name) = $_ =~ /gene_name."?([^";]+)"?/;
			my ($gene_id) = $_ =~ /gene_id."?([^";]+)"?/;
			my ($tr_id) = $_ =~ /transcript_id."?([^";]+)"?/;

			my $name = ($gene_name) ? $gene_name : $gene_id;
		 	
		    $gene_id = ($strand eq "+") ? $region.":FAM+".$stop: $region.":FAM-".$start;

		 	$annotated->{$gene_id}->{$tr_id}->{$tr_id}->{start} = $start;
		 	$annotated->{$gene_id}->{$tr_id}->{stop} = $stop;
		 	$annotated->{$gene_id}->{$tr_id}->{biotype} = $biotype;
		 	$annotated->{$gene_id}->{$tr_id}->{region} = $region;
		 	$annotated->{$gene_id}->{$tr_id}->{strand} = $strand;
		 	$annotated->{$gene_id}->{$tr_id}->{name} = $name;
		}
	}
	close F;


	return $annotated;
}

sub uninitialized_param {
	my ($v) = @_;
	not ( defined($v) and length $v );
}


sub timer {
	my $startRun = shift;
	my $endRun 	= time();
	my $runTime = $endRun - $startRun;
	printf("\nTotal running time: %02d:%02d:%02d\n\n", int($runTime / 3600), int(($runTime  % 3600) / 60), int($runTime % 60));
}




##########################################


sub find_all_ORFs_old {
	
	my $genome 		= $_[0];
    my $ORFs 		= {};
	my $count_orf	= 0;
	my $max_length 	= 21;

	foreach my $region (keys %$genome) {
		my $directStrand = $genome->{$region};
		my $reverseComplement = revdnacomp($directStrand);

		for (my $i = 0; $i < 2; $i = $i + 1) {
			my $sequenceEntry = "";
			my $strand = "";
			if ($i == 0) {		# forward strand.
				$sequenceEntry = $directStrand;
				$strand = "+";
			} else {			# reverse strand
				$sequenceEntry = $reverseComplement;
				$strand = "-";
			}

			my $starts = {};
			foreach my $start_cdn (keys %$start_codons) {
				my $offset = 0;
				my $start = index($sequenceEntry, $start_cdn, $offset);
				$starts->{$start} = 1;
				while ($start != -1) {
					$offset = $start + 1;
					$start = index($sequenceEntry, $start_cdn, $offset);
					$starts->{$start} = 1;
				}
			}

			# check if there is a corresponding stop to all start positions
			foreach my $start (sort {$a<=>$b} keys %$starts) {

				my $stop;
				my $aa_seq = "";
				my $tr_seq = "";
				my $strt_cdn = uc(substr($sequenceEntry, $start, 3));
				for (my $j = $start; $j <= (length($sequenceEntry) - 3); $j = $j + 3) {
					my $codon = uc(substr($sequenceEntry, $j, 3));
					last if (length($codon) < 3);
					my $aa = $translationHash{$codon};
					if ($aa eq "*") {
						$stop = $j;
						last;
					}
					$aa_seq = $aa_seq.$aa;
					$tr_seq = $tr_seq.$codon;
				}

				if (defined $stop and ($stop - $start) >= $MINORF ) {	
					my $length = $stop - $start;
					$aa_seq =~ s/^./M/;
					if ($strand eq '-') {
						my $start_tmp 	= length($directStrand) - $stop + 1;
						$stop 			= length($directStrand) - $start;
						$start = $start_tmp;
					} else {
						$start += 1;
					}

					my $id = $region.":".$start."-".$stop;
					$ORFs->{$id}->{region} 		= $region;
					$ORFs->{$id}->{start} 		= $start;
					$ORFs->{$id}->{stop} 		= $stop;
					$ORFs->{$id}->{len} 		= $length;
					$ORFs->{$id}->{strand} 		= $strand;
					$ORFs->{$id}->{tr_seq} 		= $tr_seq;
					$ORFs->{$id}->{aa_seq} 		= $aa_seq;
					$ORFs->{$id}->{start_codon} = $strt_cdn;
					$ORFs->{$id}->{gene} 		= ($strand eq "+") ? $region.":FAM+".$stop: $region.":FAM-".$start;;

					if ($max_length < $length) {$max_length = $length;}
					$count_orf++;
				}
			}
		}
	}
	print "Total number of ORFs $count_orf maximum ORF length $max_length\n";
	return $ORFs;
}


