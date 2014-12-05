#!/usr/bin/perl
##############################################
#### Shifts minus strand reads by read length
#### Perl script		
#############################################



use feature qw(say);


# ================== 
# Parsing Arguments
# ================== 

#initialize to NULL
my $bed_file = NULL; #bed file
my $shift_size = NULL; #shift size
my $read_length = NULL; #read length

#Parse the arguments
$bed_file = $ARGV[0]; #bed file
$shift_size = $ARGV[1];	#shift size
$read_length = $ARGV[2]; #read length

#=======================> DONE! 



# ========================================
# Parse the bed file and extend the reads
# ========================================

#open the file
open(DATA, $bed_file) || die("Can't open the bed file, probably you gave me the wrong path!");

#loop through the rest of the file line by line (To do: look for a faster way)
while (<DATA>) {
    my ($start, $strand) = split(/\t/,$_,2);
    $strand =~ s/\015?\012?$//;


    
	#plus strand
	if ($strand eq '+') {
		#do nothing
	}
	#minus strand
	elsif ($strand eq '-') {
		$start = $start - $shift_size + $read_length; #difference between $start and $end should be equal to fragment length ($shift_size)
	}
	#bad format 
	else {
		die("It appears the bed file is not formatted properly. Specifically, I expect to find either + or - in the second column, but I found something else!");	
	}
	if ($start >= 0) {
		#Now write the new line
		say join "\t", $start, $strand;
	}
}
#=======================> DONE! 
