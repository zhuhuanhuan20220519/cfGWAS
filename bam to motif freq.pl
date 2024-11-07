#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw/min max/;

my $samtools="samtools";    #path to samtools
my $fa = "hg38.p13.fa";     # reference genome
my $motif_list = "4mer.motif.list";       #list of 256 types of 4mer motif

# Step 1: Filter raw BAM file to obtain clean SAM file
# Usage: perl script.pl <input.bam> <min_map_qual> <max_mismatch> <motif_length> <sample_prefix>

if (@ARGV != 5) {
    die "Usage: perl script.pl <input.bam> <min_map_qual> <max_mismatch> <motif_length> <sample_prefix>\n";
}

my ($bam_file, $min_map_qual, $max_mismatch, $motif_length, $sample_prefix) = @ARGV;
my %hash;
my %hash_pair;

# Open BAM file and filter reads based on quality, mismatches, and alignment flags
open my $bam_in, "-|", "$samtools view $bam_file" or die "Could not open BAM file: $!\n";

while (<$bam_in>) {
    chomp;
    my $line = $_;
    my @F = split /\s+/, $line;

    # Apply filtering conditions
    next if ($F[6] ne "=");                  # Only keep pairs mapped to the same chromosome
    next if ($F[4] < $min_map_qual);         # Mapping quality cutoff
    next if ($F[1] & 0x0400);                # Exclude PCR/optical duplicates
    next if ($F[1] & 0x4);                   # Exclude unmapped segments
    next if ($F[1] & 0x8);                   # Exclude pairs with unmapped segments
    next if ($F[1] & 0x10 && $F[1] & 0x20);  # Exclude if both reads are reverse complemented
    next unless ($F[1] & 0x10 || $F[1] & 0x20); # Keep if only one read is reverse complemented
    next if ($F[1] & 0x100);                 # Exclude secondary alignments
    next if ($F[8] > 600 || $F[8] < -600);   # Exclude fragments >600 bp in length
    next if (exists $hash{$F[0]});           # Exclude fragments with multiple mismatches

    # Store pairs for further filtering and check mismatches
    $hash_pair{$F[0]}++;
    if ($hash_pair{$F[0]} == 2) {
        delete $hash_pair{$F[0]};
    }

    my $distance = 10;  # Initialize mismatch distance
    if ($F[5] !~ /^(\d+)M$/) {
        $hash{$F[0]}++;
    }
    if ($line =~ /\s+NM:i:(\d+)\s+/) {
        $distance = $1;
        if ($distance > $max_mismatch) {
            $hash{$F[0]}++;
        }
    }
}
close $bam_in;

open $bam_in, "-|", "$samtools view $bam_file" or die "Could not open BAM file: $!\n";
open my $sam_out, ">", "$sample_prefix.filtered.sam" or die "Could not open output SAM file: $!\n";
while (<$bam_in>) {
    chomp;
    my $line = $_;
    my @F = split /\s+/, $line;

    # Apply filtering conditions
    next if ($F[6] ne "=");                  # Only keep pairs mapped to the same chromosome
    next if ($F[4] < $min_map_qual);         # Mapping quality cutoff
    next if ($F[1] & 0x0400);                # Exclude PCR/optical duplicates
    next if ($F[1] & 0x4);                   # Exclude unmapped segments
    next if ($F[1] & 0x8);                   # Exclude pairs with unmapped segments
    next if ($F[1] & 0x10 && $F[1] & 0x20);  # Exclude if both reads are reverse complemented
    next unless ($F[1] & 0x10 || $F[1] & 0x20); # Keep if only one read is reverse complemented
    next if ($F[1] & 0x100);                 # Exclude secondary alignments
    next if ($F[8] > 600 || $F[8] < -600);   # Exclude fragments >600 bp in length
    next if (exists $hash{$F[0]} || exists $hash_pair{$F[0]});          # Exclude fragments with multiple mismatches

    # Final check and output clean SAM lines
    print $sam_out "$line\n";
}

close $bam_in;
close $sam_out;

# Step 2: Calculate fragment lengths for filtered reads from SAM file and output in BED format
open my $bed_out, ">", "$sample_prefix.filtered.bed" or die "Could not open output bed file: $!\n";
open my $sam_in, "<", "$sample_prefix.filtered.sam" or die "Could not open filtered SAM file: $!\n";

while (<$sam_in>) {
    chomp;
    my @F = split /\s+/;
    next if ($F[6] ne "=");
    if ($F[8] < 0 && $F[3] >= $F[7]) {
        my $s = min($F[3], $F[7]) - 1;
        my $len = abs($F[3] - $F[7]) + length($F[9]);
        my $e = $s + $len;
        print $bed_out join("\t", $F[2], $s, $e, $len), "\n";
    }
}

close $sam_in;
close $bed_out;

# Step 3: Extract motif sequences from reference genome

my %hg;
open my $HG, '<', $fa or die "Could not open genome file '$fa': $!\n";
my $S = "";

while (<$HG>) {
    chomp;
    if (/^>/) {
        $_ =~ s/>//g;
        my @fields = split;
        $S = $fields[0];
        $hg{$S} = "";
    } else {
        $hg{$S} .= uc $_;
    }
}
close $HG;

open my $bed_in, '<', "$sample_prefix.filtered.bed" or die "Could not open filtered bed file: $!\n";
open my $motif_out, ">", "$sample_prefix.motif.bed" or die "Could not open motif output file: $!\n";

my $for_count=0;
my $rev_count=0;

while (<$bed_in>) {
    chomp;
    my $line=$_;
    my @F=split/\s+/,$line;
    my $start_forward = $F[1] - $motif_length;
    next if ($F[0] eq "chrM" && $start_forward < 0);
    my $E5 = substr($hg{$F[0]}, $F[1], $motif_length);
    my $x = $F[2] - $motif_length;
    my $E3 = substr($hg{$F[0]}, $x, $motif_length);
    my $com_seq = $E5 . $E3;
    if ($com_seq!~/N/){
      print $motif_out join("\t", $line, $E5, $E3)."\n";
      }
}
close $bed_in;
close $motif_out;

# Step 4: Calculate motif frequencies
my %hash_motif;
my @MER;
open my $motif_in, '<', $motif_list or die "Could not open motif list: $!\n";
while (<$motif_in>) {
    chomp;
    next if (/Type/);
    push @MER, $_;
}
close $motif_in;

my $sum_site = 0;
open my $motif_file, '<', "$sample_prefix.motif.bed" or die "Could not open motif bed file: $!\n";
while (<$motif_file>) {
    chomp;
    my @F = split;
    $sum_site += 2;
    my $motif_up = substr($F[4], 0, 4);
    my $motif_down=revcom($F[5]);
    $motif_down=substr($motif_down,0,4);
    $hash_motif{$motif_up}++;
    $hash_motif{$motif_down}++;

}

close $motif_file;

open my $motif_freq_out, '>', "$sample_prefix.4mer.motif.freq" or die "Could not open motif frequency file: $!\n";
print $motif_freq_out "Type\t$sample_prefix\n";
foreach my $i(@MER){
    if (exists $hash_motif{$i}){
            my $ratio=$hash_motif{$i}/$sum_site;
            print $motif_freq_out "$i\t$ratio\n";
                                        }
        else{
            print $motif_freq_out "$i\t0\n";
        }
}
close $motif_freq_out;

# Step 5: Calculate motif diversity score
open my $freq_in, '<', "$sample_prefix.4mer.motif.freq" or die "Could not open motif frequency file: $!\n";
open my $mds_out, '>', "$sample_prefix.4mer.motif.MDS" or die "Could not open motif frequency file for diversity: $!\n";

my $divScore = 0;
while (<$freq_in>) {
    chomp;
    next if (/Type/);
    my $line=$_;
   my @F=split/\s+/;
   if ($F[1]==0){
	$F[1]=10**-10;
   }
   my $motif_score=-$F[1]*log($F[1])/log(256);
   $divScore+=$motif_score;
}
close $freq_in;

print $mds_out "Sample\t$sample_prefix\tDiversity Score: $divScore\n";


# Helper function for reverse complement
sub revcom {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return scalar reverse($seq);
}
