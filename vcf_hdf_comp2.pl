#!/usr/bin/perl -w

# nohup ./vcf_hdf_comp2.pl /data/seq/processed-genomes/Trio/17144784_rtg/snps/snps.vcf.gz /data/gs/snp2gw_LICHDF/2011-08-12-sequenced/seqgencon.lichdf 17144784 & 
# Number missing SNPs from VCF (? by chromosome)
# Per chromosome for counts
# Count of non-passing SNPs

# Total number of SNPs in VCF file
# Percentage of SNPs in diagonal of cont-table, over whole table

# /data/seq/chhar0/rtg-working/snp/ 99591_GATK_SNPs.vcf.gz 99591_rtg_snps.vcf.gz 
# /data/seq/processed_genomes/SNPs/ bwa novoalign

# prefix option
# Marker qualities + concordance => csv file

use strict;
use warnings;
use IO::Uncompress::Gunzip;
use Getopt::Long;

#my ($vcf, $hdf, $animal, $prefix) = @ARGV;

my $vcf = undef;
my $hdf = undef;
my $animal = undef;
my $prefix = 'comparer.out';
my $minQuality = 0;
my $minDepth = 0;
my $smoothing = 300;

GetOptions(
	'vcf=s' => \$vcf,
	'hdf=s' => \$hdf,
	'animal=i' => \$animal,
	'prefix=s' => \$prefix,
	'quality=f' => \$minQuality,
	'depth=f' => \$minDepth,
	'smoothing=i' => \$smoothing
);

die "Usage: $0 --vcf=vcf_file --hdf=hdf_file --animal=animal_key [--prefix=prefix] [--quality=q] [--depth=d] [--smoothing=s]\n" 
	unless defined $vcf and defined $hdf and defined $animal and $animal =~ /^\d+$/;

print "Comparing VCF file $vcf with LICHDF file $hdf, for animal $animal...\n";

$prefix = 'comparer.out' unless defined $prefix and $prefix =~ /^[A-Za-z0-9\-_]+$/;
print "Using '$prefix' as the prefix for output files.\n";

$minQuality = 0 unless defined $minQuality and $minQuality > 0;
$minDepth = 0 unless defined $minDepth and $minDepth > 0;
$smoothing = 300 unless defined $smoothing and $smoothing > 0;

#############################################################################################

my $z = new IO::Uncompress::Gunzip $vcf, {'MultiStream' => 1} or die "Cannot open the VCF file: $!\n";

### Read the VCF file, and store the lines
print "Reading the VCF file...\n";

my %vcfLines = ();
my $countPassSNPs = 0;
my $countFailSNPs = 0;

while (my $line = <$z>) {
	chomp $line;
	next if $line =~ /^#/ or $line =~ /^\s*$/;
	my ($vcfChr, $vcfPos, undef, $a, $b, $vcfQual, $vcfPass, undef, $vcfKey, $vcfStuff) = split "\t", $line;
	#($vcfChr) = ($vcfChr =~ /Chr([A-Z0-9]+)/);
	
	my @keyParts = split ':', $vcfKey;
	my @stuffParts = split ':', $vcfStuff;
	
	my $gt = undef;
	my $dp = undef;
	for (my $i = 0; $i < @keyParts; $i++) {
		$gt = $stuffParts[$i] if $keyParts[$i] eq 'GT';
		$dp = $stuffParts[$i] if $keyParts[$i] eq 'DP';
		last if defined $gt and defined $dp;
	}
	
	$vcfChr = substr $vcfChr, 3; # chop the 'Chr' off the beginning
	
	if ($vcfPass eq 'PASS') {
		$countPassSNPs++;
	}
	else {
		$countFailSNPs++;
	}
	
	$vcfLines{"$vcfChr:$vcfPos"} = [$vcfQual, $vcfPass, $gt, $dp, $a, $b];
}
close $z;
my $numVCF = scalar keys %vcfLines;
print "$countPassSNPs pass; $countFailSNPs fail by VCF filter\n";

#############################################################################################


### Read the HDF file
print "Reading the LICHDF file...\n";
my ($rootHDFFile) = ($hdf =~ /(.*)\.lichdf$/);
die "Invalid HDF file.\n" unless defined $rootHDFFile;

my $animalFile = "$rootHDFFile.animals.lichdf";
my $markerFile = "$rootHDFFile.markers.lichdf";
my $genoFile = "$rootHDFFile.genotypes.lichdf";

#############################################################################################


# Read the animal file
my $animalLineNum = 0;
my $found = 0;
open ANIMAL, $animalFile or die "Cannot open the LICHDF animal file: $!\n";
while (my $line = <ANIMAL>) {
	chomp $line;
	$animalLineNum++;
	if ($line =~ /$animal/) {
		$found = 1;
		last;
	}
}
close ANIMAL;

die "Animal $animal is not in the LICHDF file.\n" unless $found;

#############################################################################################


my %markerLines = ();
my $markerLine = 0;
open MARKER, $markerFile or die "Cannot open the LICHDF marker file: $!\n";
while (my $line = <MARKER>) {
	chomp $line;
	my ($id, $chr, $pos, $a, $b) = split /\s+/, $line;
	
	$markerLines{"$chr:$pos:$id"} = [$markerLine, $a, $b];
	$markerLine++;
	
}
close MARKER;

my $numMarkers = scalar keys %markerLines;
die "The marker file contains no valid markers.\n" if $numMarkers < 1;

#############################################################################################

my $skip = ($animalLineNum-1)*$numMarkers;
print "Reading markers from position $skip, for $numMarkers bytes.\n";
open GENO, $genoFile or die "Cannot read the genotype file: $!\n";
binmode GENO;
seek GENO, $skip, 0 or die "The genotype file is the wrong size: $!\n";

my $genoBuffer = '';
read GENO, $genoBuffer, $numMarkers;
close GENO;
my @genotypes = split '', $genoBuffer;

die "Error!!!!!" unless @genotypes == $numMarkers;

#open GENOOUT, ">$prefix.genotype" or die "Cannot write genotype out: $!\n";
#for my $g (@genotypes) {
#	print GENOOUT chr(ord($g)+48);
#}
#close GENOOUT;

#############################################################################################


open OUTPUT, ">$prefix.found.csv" or die "Cannot open output found-marker file: $!\n";
print OUTPUT "marker,chr,pos,vcfQual,hetHDF,hetVCF,sameHet,allelesHDF,allelesVCF,markerNumber\n";

open MISSING, ">$prefix.missing.csv" or die "Cannot open missing-marker file: $!\n";
print MISSING "marker,chr,pos,hetHDF,allelesHDF,markerNumber\n";

open MARKERINFO, ">$prefix.markers.csv" or die "Cannot open moutput marker-info file: $!\n";
print MARKERINFO "marker,vcfQual,vcfPass,coverage,concordant\n";

my $numAgree = 0;
my $numDisagree = 0;
my $matchLowQual = 0;
my $numMissing = 0;

my %numAgreeByChromo = ();
my %numDisagreeByChromo = ();

my $numHomoMissing = 0;

my $numFilteredQuality = 0;
my $numFilteredDepth = 0;

my %contTable = ();

for my $key (sort byPos keys %markerLines) {

	my ($markerChr, $markerPos, $markerId) = split ':', $key;
	my ($markerLine, $aM, $bM) = @{$markerLines{$key}};
	my $genotype = $genotypes[$markerLine];
	$genotype = chr(ord($genotype)+48);

	next if $markerPos == 0 or $markerChr eq 'Y' or $markerChr eq '0';

	my $hetMarkers = 'false';
	$hetMarkers = 'true' if $genotype eq '2';

	if (defined $vcfLines{"$markerChr:$markerPos"}) {
		
		my ($vcfQual, $vcfPass, $vcfGeno, $vcfDepth, $aV, $bV) = @{$vcfLines{"$markerChr:$markerPos"}};
		
		# Check against filters
		if ($vcfQual =~ /\d+\.?\d+/ && $vcfQual < $minQuality) {
			$numFilteredQuality++;
			next;
		}
		if ($vcfDepth =~ /\d+/ && $vcfDepth < $minDepth) {
			$numFilteredDepth++;
			next;
		}
		
		# Test for concordance
		my @genoParts = split /\/|\|/, $vcfGeno;
		my $hetVCF = 'false';
		$hetVCF = 'true' unless $genoParts[0] eq $genoParts[1];
		my $hetSame = 'false';
		$hetSame = 'true' if $hetVCF eq $hetMarkers;
		if ($hetSame eq 'true') {
			$numAgreeByChromo{$markerChr}++;
			$numAgree++;
		}
		else {
			$numDisagreeByChromo{$markerChr}++;
			$numDisagree++;
		}
		
		# Test for pass/fail status (only works for RTG files)
		if ($vcfPass ne 'PASS') {
			$matchLowQual++;
		}
		
		# Calculate the contingency table
		if ($hetMarkers eq 'true' && $hetVCF eq 'true') {
			$contTable{'HETHET'}++;
		}
		elsif ($hetMarkers eq 'true') {
			$contTable{'HETM'}++;
		}
		elsif ($hetVCF eq 'true') {
			$contTable{'HETV'}++;
		}
		else {
			$contTable{'HOMHOM'}++;
		}
		
		# Replace commas with pluses to avoid messing up CSV files
		$aV =~ s/,/+/g;
		$bV =~ s/,/+/g;
		
		print OUTPUT "$markerId,$markerChr,$markerPos,$vcfQual,$hetMarkers,$hetVCF,$hetSame,$aM;$bM,$aV;$bV,$markerLine\n";
		print MARKERINFO "$markerId,$vcfQual,$vcfPass,$vcfDepth,$hetSame\n";
		
	}
	else {
		
		print MISSING "$markerId,$markerChr,$markerPos,$hetMarkers,$aM;$bM,$markerLine\n";
		$numMissing++;
		$numHomoMissing++ if $hetMarkers eq 'false';
	}
}

close OUTPUT;
close MISSING;
close MARKERINFO;

#############################################################################################

print "Done.\n";

my $ratio = $numAgree/$numDisagree;
$ratio = 0 unless defined $ratio;
$ratio = sprintf "%.2f", $ratio;

my $pctMissing = $numMissing/$numMarkers*100;
$pctMissing = 0 unless defined $pctMissing;
$pctMissing = sprintf "%.2f%%", $pctMissing;

my $numHetMissing = $numMissing - $numHomoMissing;
my $pctHomoMissing = $numHomoMissing/$numMissing*100;
$pctHomoMissing = 0 unless defined $pctHomoMissing;
$pctHomoMissing = sprintf "%.2f%%", $pctHomoMissing;
my $pctHetMissing = $numHetMissing/$numMissing*100;
$pctHetMissing = 0 unless defined $pctHetMissing;
$pctHetMissing = sprintf "%.2f%%", $pctHetMissing;

#############################################################################################

$contTable{'HETHET'} = 0 unless defined $contTable{'HETHET'};
$contTable{'HETV'} = 0 unless defined $contTable{'HETV'};
$contTable{'HETM'} = 0 unless defined $contTable{'HETM'};
$contTable{'HOMHOM'} = 0 unless defined $contTable{'HOMHOM'};

my $pctRight = ($contTable{'HETHET'} + $contTable{'HOMHOM'}) / ($contTable{'HETHET'} + $contTable{'HETV'} + $contTable{'HETM'} + $contTable{'HOMHOM'}) * 100;
$pctRight = 0 unless defined $pctRight;
$pctRight = sprintf "%.2f", $pctRight;

open CONT, ">$prefix.conttable" or die "Cannot save the contigency table: $!\n";

print CONT "Total number of SNPs in the LICHDF file: $numMarkers\n";
print CONT "Total number of variants in the VCF file: $numVCF\n";
print CONT "Number filtered for low quality (<$minQuality): $numFilteredQuality\n" if $minQuality > 0;
print CONT "Number filtered for read depth (<$minDepth): $numFilteredDepth\n" if $minDepth > 0;

print CONT "\nContingency Table:\n\n";
print CONT "        SNP Genotype\n";
print CONT "        |   Het    |   Homo   \n";
print CONT "  ------+----------+----------\n";
printf CONT "S  Het  |  %6s  |  %6s  \n", $contTable{'HETHET'}, $contTable{'HETV'};
print CONT "e ------+----------+----------\n";
printf CONT "q  Homo |  %6s  |  %6s  \n", $contTable{'HETM'}, $contTable{'HOMHOM'};
print CONT "\n\n";
print CONT "Number agree: $numAgree\t\t\tNumber disagree: $numDisagree\nNumber low quality: $matchLowQual\t\tNumber missing from VCF: $numMissing\n";
print CONT "Agree/Disagree Ratio: $ratio\tPercentage Missing: $pctMissing\n";
print CONT "Number of missing homozygotes: $numHomoMissing ($pctHomoMissing)\n";
print CONT "Percentage of concordant markers: $pctRight\n";

print CONT "\n ** Statistics with missing homozygous SNPs assumed to be correct **\n";

$numAgree += $numHomoMissing;
$numMissing = $numHetMissing;

print CONT "\nContingency Table (assuming missing homozygotes agree):\n\n";
print CONT "        SNP Genotype\n";
print CONT "        |   Het    |   Homo   \n";
print CONT "  ------+----------+----------\n";
printf CONT "S  Het  |  %6s  |  %6s  \n", $contTable{'HETHET'}, $contTable{'HETV'};
print CONT "e ------+----------+----------\n";
printf CONT "q  Homo |  %6s  |  %6s  \n", $contTable{'HETM'}, $contTable{'HOMHOM'}+$numHomoMissing;
print CONT "\n\n";

print CONT "Number agree: $numAgree\t\tNumber missing from VCF: $numMissing\n";

$ratio = $numAgree/$numDisagree;
$ratio = 0 unless defined $ratio;
$ratio = sprintf "%.2f", $ratio;

$pctMissing = $numMissing/$numMarkers*100;
$pctMissing = 0 unless defined $pctMissing;
$pctMissing = sprintf "%.2f%%", $pctMissing;

$pctRight = ($contTable{'HETHET'} + $contTable{'HOMHOM'} + $numHomoMissing) / ($contTable{'HETHET'} + $contTable{'HETV'} + $contTable{'HETM'} + $contTable{'HOMHOM'} + $numHomoMissing) * 100;
$pctRight = 0 unless defined $pctRight;
$pctRight = sprintf "%.2f", $pctRight;


print CONT "Agree/Disagree Ratio: $ratio\tPct Missing: $pctMissing\n";
print CONT "Percentage of concordant markers: $pctRight\n";

close CONT;

#############################################################################################

open CSV, ">$prefix.csv" or die "Cannot open the CSV output file: $!\n";

my %keys = ();
for my $key (keys %numAgreeByChromo) {
	$keys{$key} = 1;
}
for my $key (keys %numDisagreeByChromo) {
	$keys{$key} = 1;
}

print CSV "Chromosome,Number_agree,Number_disagree,Total,Percent_agree\n";
for my $chromo (sort byChromo keys %keys) {
	my $numAgree = $numAgreeByChromo{$chromo} if defined $numAgreeByChromo{$chromo};
	my $numDisagree = $numDisagreeByChromo{$chromo} if defined $numDisagreeByChromo{$chromo};
	$numAgree = 0 unless defined $numAgree;
	$numDisagree = 0 unless defined $numDisagree;
	my $total = $numAgree + $numDisagree;
	my $pct = $numAgree*100.0/$total if $total > 0;
	$pct = 0 unless defined $pct;
	
	$pct = sprintf "%.2f", $pct;
	
	print CSV "$chromo,$numAgree,$numDisagree,$total,$pct\n";
}

close CSV;

#############################################################################################

print "Generating density plot using R...\n";

my $rcode = "data = read.csv('$prefix.markers.csv');
noncon = data[which(data\$concordant == 'false'),];
con = data[which(data\$concordant == 'true'),];

dN = density(log10(noncon\$vcfQual), n=$smoothing);
dC = density(log10(con\$vcfQual), n=$smoothing);
nN = nrow(noncon);
nC = nrow(con);

pdf(file='$prefix.density.pdf', width=10, height=7, title='Variant Qualities from VCF file $vcf');
plot(dC, xaxt='n', col='blue', lwd=2, xlab='Variant Quality (phred scale)', ylab='Density', ylim=c(0,floor(max(c(dC\$y,dN\$y)+0.1)*10)/10),
	main='Qualities for concordant and non-concordant Variants');
axis(side=1, at=1:6, label=10^(1:6));
lines(dN, col='red', lwd=2);
abline(v=log10(300),col='green');
legend(x='topleft', legend=c(sprintf('Concordant (\%d)',nC), sprintf('Non-concordant (\%d)',nN) ), 
	lwd=2, col=c('blue','red'));
dev.off();";

open RCODE, ">$prefix.rcode.R" or die "Cannot write R code: $!\n";
print RCODE $rcode;
close RCODE;

`R CMD BATCH $prefix.rcode.R`;

#############################################################################################

sub byPos {
	my ($ac,$ap,undef) = split ':', $a;
	my ($bc,$bp,undef) = split ':', $b;
	
	if ($ac ne $bc) {
		if ($ac =~ /^\d+$/ && $bc =~ /^\d+$/) {
			return $ac <=> $bc; # numeric comparison
		}
		else {
			return $ac cmp $bc; # string comparison
		}		
	}
	else {
		return $ap <=> $bp;
	}
}

sub byChromo {
	if ($a =~ /^\d+$/ && $b =~ /^\d+$/) {
		return $a <=> $b; # numeric comparison
	}
	else {
		return $a cmp $b; # string comparison
	}
}

