#!/usr/bin/perl
########################################
#### Filters peaks by p-value and width
#### Perl script		
#######################################


use feature qw(say);


# ================== 
# Parsing Arguments
# ================== 
#initialize to NULL
my $bed_file = NULL; #bed file

#Parse the arguments
$bed_file = $ARGV[0]; #bed file
#=======================> DONE! 




# ========================================
# Parse the bed file and extend the reads
# ========================================
#open the file
open(DATA, $bed_file) || die("Can't open the bed file, probably you gave me the wrong path!");

while (<DATA>) {
    my ($chr, $start, $end, $name, $score, $strand, $signal, $pvalue, $qvalue, $summit, $minpeak, $geom) = split(/\t/,$_,12);
    
	my $size = $end - $start;
	
	if($size >= $minpeak) {
		if($signal > $geom) {
			#Now write the new line
			$summit =~ s/\015?\012?$//;
			say join "\t", $chr, $start, $end, $name, $score, $strand, $signal, $pvalue, $qvalue, $summit;
		}
	}
}
#=======================> DONE! 
