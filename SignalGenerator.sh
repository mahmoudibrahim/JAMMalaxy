#!/usr/bin/bash




##Finding out the path
sPath="`dirname \"$0\"`"
sPath="`( cd \"$sPath\" && pwd )`"



usage()
{
cat << EOF
Welcome to JAMM v1.0.6rev3 Signal Generator Script (GNU GPLv3). Copyright (C) 2014  Mahmoud Ibrahim.

This program comes with ABSOLUTELY NO WARRANTY; for details visit http://www.gnu.org/licenses/gpl.html. This is free software, and you are welcome to redistribute it under certain conditions; visit http://www.gnu.org/licenses/gpl.html for details.

OPTIONS:
   -s      Directory containing sample files (required)
   -g      Genome size file (required)
   -o      Output Directory (required)
   -c      directory containing input or Control files
   -r 	   file with Regions to get signal for (required)
   -b      Bin size for signal generation (default: 10)
   -f      Fragment lengths (required)
   -p	   Number of processors used by R scripts (default: 1)
 
EOF
}


# ========================= 
# Process Input parameters
# =========================

#Defaults -- Change those if you want
export LANG=C #locale defaults
export LC_ALL=C #locale defaults

#Defaults -- Do not change
sdir=""
gsize=""
out=""
signal="10"
regs=""
fraglen=""
wdir=$(mktemp -d)
ran=$RANDOM
cores="1"


while getopts "s:o:c:r:b:f:g:p:" OPTION
do
	case $OPTION in
	s) sdir=$OPTARG
	;;
	g) gsize=$OPTARG
	;;
	o) out=$OPTARG
	;;	
	c) bdir=$OPTARG
	;;
	r) regs=$OPTARG
	;;
	b) signal=$OPTARG
	;;
	f) fraglen=$OPTARG
	;;
	p) cores=$OPTARG
	;;
	?)
	usage
	exit
	;;
	esac
done


if [[ -z $fraglen ]] || [[ -z $sdir ]] || [[ -z $regs ]] || [[ -z $out ]]
then
     usage
     exit 1
fi
#=======================> DONE!



# ============================= 
# Step One: Initial Processing
# =============================
printf "\n\n============================================\nStarted JAMM Signal Generator Pipeline...Hang on!\n============================================\n\n"

if [ ! -d "$wdir" ]; then
	mkdir $wdir #make working directory
fi
if [ ! -d "$out" ]; then
	mkdir $out #make working directory
fi
mkdir $wdir/bkgd.$ran/ #directory to store background files
mkdir $wdir/sizes.$ran/ #chromosomes and sizes
mkdir $wdir/samples.$ran/ #store sample files

dupnum=$(ls -1 $sdir | wc -l) #count how many sample files


#separate chromosome sizes
printf "Loading genome size file..."
ext="$wdir/sizes.$ran/"
awk -v ext="$ext" '{ print >> ext"/size." $1 ".bed" }' $gsize
printf "Done!\n"


printf "Processing sample files..."
#load each chromosome from each sample file
for i in $sdir/*.bed; do
samplefile=$(basename $i)	
	for f in $wdir/sizes.$ran/*; do
		sizefile=$(basename $f)
		chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');
		awk -v chr="$chr" -v ext="$wdir/samples.$ran/" -v samplefile="$samplefile" -F"\t" '$1 == chr { print $2"\t"$6 >> ext"sample."chr"."samplefile }' "$i" 
	done
done
printf "Done!\n"


if [ ! -z $bdir ]; then
#concatenate all background files into one file
printf "Processing control files..."
cat $bdir/*.bed > $wdir/bkgd.$ran/ctrl.bed

for f in $wdir/sizes.$ran/*; do
	sizefile=$(basename $f)
	chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');
	awk -v chr="$chr" -v ext="$wdir/bkgd.$ran/" -F"\t" '$1 == chr { print $2"\t"$6 >> ext"bkgd."chr".ctrl.bed" }' "$wdir/bkgd.$ran/ctrl.bed"
done

rm $wdir/bkgd.$ran/ctrl.bed

printf "Done!\n"
fi


#determine average read lengths
printf "Getting average read lengths..."
	if [ ! -z $bdir ]; then
		readC=$(awk '{a=$3-$2;print a;}' "$wdir/bkgd.$ran/ctrl.bed" | perl -lane '$a+=$_;END{print $a/$.}' | awk '{print int($1)}')
	fi
	readL=""
	for s in $sdir/*.bed; do #and for each sample file
		read=$(awk '{a=$3-$2;print a;}' "$s" | perl -lane '$a+=$_;END{print $a/$.}' | awk '{print int($1)}')
		readL="$readL,$read"
	done
	readL=${readL#","}
printf "Done!\n"
#=======================> DONE!






# =========================== 
# Step Three: Getting Signal
# ===========================
mkdir $wdir/signal.$ran/
mkdir $out/signal #store signal


printf "Generating Signal...(bin size: $signal)\n"

counting=1;			
for f in $wdir/sizes.$ran/*; do #for each chromosome
		samplelist=""
		frag=""
		frag=$fraglen
		k=1
		
		sizefile=$(basename $f)
		chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');
		chrSize=$(cat $f | cut -f2);
		printf "Chromosome $chr: "
		
		#list of sample bed files and fragment lengths
		for s in $wdir/samples.$ran/*.bed; do #and for each sample file
			samplefile=$(basename $s)
			chr2=$(echo $samplefile | awk -F"." '{print $2}');
			if [ $chr == $chr2 ] #belonging to this chromosome
			then
				samplelist="$samplelist,$wdir/samples.$ran/ext.$samplefile"
				samplename=$(echo $samplefile | awk -F"." '{ print $3 }')
				samplefilename=$(echo $samplefile | cut -d'.' -f 3-)
				shift=$(echo "$frag" | cut -f "$k" -d ",")
				read=$(echo "$readL" | cut -f "$k" -d ",")	
				k=$(($k+1))
				perl "$sPath/readshifter.pl" "$wdir/samples.$ran/$samplefile" $shift $read > "$wdir/samples.$ran/ext.$samplefile"
			fi
		done
		
		#control file
		bkgdfile="None"
		if [ ! -z $bdir ]; then
			l=$(($dupnum+1))
			bshift=$(echo $frag | cut -f "$l" -d ",")
			perl "$sPath/readshifter.pl" "$wdir/bkgd.$ran/bkgd.$chr.ctrl.bed" $bshift $readC > "$wdir/bkgd.$ran/ext.bkgd.$chr.ctrl.bed"
			bkgdfile="$wdir/bkgd.$ran/ext.bkgd.$chr.ctrl.bed"
		fi
		
		#remove leading comma
		samplelist=${samplelist#","}
		frag=${frag#","}
		
		#call the peak calling R script
		Rscript "$sPath/signalmaker.r" -chromoS="$chrSize" -chromo="$chr" -bednames=$samplelist -frag=$frag -bkgd=$bkgdfile -out="$wdir/signal.$ran/" -sig="$signal" -regions="$regs" -p="$cores" -chrcount="$counting"
		cat $wdir/signal.$ran/*.bedSignal | awk -F"\t" -v j=0 '$4 > j' > "$out/signal/$chr.bedGraph"
		rm $wdir/signal.$ran/*.bedSignal
		counting=$(($counting+1));
done
counting=1;
#=======================> DONE!

rm -rf $wdir

printf "\n\n========================================\nWe're done...Congratulations!\n========================================\n\n"
