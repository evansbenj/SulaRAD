# Analysis with Kai's program

The first step to analyze each species/population with Kai's program is to generate an input file.  I wrote a script to do this (19_generates_input_for_Kai_program.pl) below.  The command lines I used are:

```
./19_Generates_input_for_Kai_program.pl nonrecal_51000plus.vcf.gz.gz.tab_with_baboon.tab_and_human_andchrXdepth.tab 0000000000000000000000000000000000000000 4_6_32_33_34_35_36_37_38_39_40 tonk_kai_newdepth.input

./19_Generates_input_for_Kai_program.pl nonrecal_51000plus.vcf.gz.gz.tab_with_baboon.tab_and_human_andchrXdepth.tab 0000000000000000000000000000000000000000 4_6_8_9_10_11_12_13 maura_kai_newdepth.input

./19_Generates_input_for_Kai_program.pl nonrecal_51000plus.vcf.gz.gz.tab_with_baboon.tab_and_human_andchrXdepth.tab 0000000000000000000000000000000000000000 4_6_2_3_4_5_6_7 hecki_kai_newdepth.input

./19_Generates_input_for_Kai_program.pl nonrecal_51000plus.vcf.gz.gz.tab_with_baboon.tab_and_human_andchrXdepth.tab 0000000000000000000000000000000000000000 4_6_18_19_20_21 borneo_kai_newdepth.input

```

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



#  This program reads in a tab delimited genotype file generated
#  by the perl program "17_adds_outgroup_to_lots_of_tab_files.pl"
#  or from vcftools vcf_to_tab
#  and generates an input file for analysis with Kai's program

# to execute type 19_Generates_input_for_Kai_program.pl 
# /home/ben/2015_SulaRADtag/good_merged_samples/tab_relative_to_genez/recal_51000plus.vcf.gz.tab_with_baboon.tab_and_human.tab 
# 1111100110000111100011100110010100000000 
# 3_6_22_23_25_26 nigra_kai_input.txt  
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 4_6_22_23_25_26 refers to (i) the column that contains the 
# outgroup nucleotide (4 in this case), (ii) the column number of the first individual in the ingroup 
# (6 in this case), and (iii) the sample number that contain the data from the individuals you want to 
# include (22, 23, 25, and 26 in this case), which are the four nigra samples itemized below.

# IMPORTANT: (i) and (ii) are columns beginning with 1 but (iii) is based on the individual samples such
# as enumerated below

# Notes: 
# ##### nigra_PM1000 is actually nigrescens_PM1000 #####

# ##### these samples have very low coverage (<10x): ####
# ##### nigrescens_PF654_sorted (7.33X) ####
# ##### maura_PM613_sorted (8.65X) ####
# ##### ochreata_PM596_sorted (9.02X)####
# ##### nigra_660_sorted (9.61X) ####
# ##### togeanus_PF549 (9.63X) ####

# Here is the order of the samples for the SulaRad project:


#   1   brunescens_PF707_stampy_sorted
#   2   hecki_PF643_stampy_sorted
#   3   hecki_PF644_stampy_sorted
#   4   hecki_PF648_stampy_sorted
#   5   hecki_PF651_stampy_sorted
#   6   hecki_PM639_stampy_sorted
#   7   hecki_PM645_stampy_sorted
#   8   maura_PF615_stampy_sorted
#   9   maura_PF713_stampy_sorted
#   10  maura_PM613_stampy_sorted
#   11  maura_PM614_stampy_sorted
#   12  maura_PM616_stampy_sorted
#   13  maura_PM618_stampy_sorted
#   14  nem_Gumgum_stampy_sorted
#   15  nem_Kedurang_stampy_sorted
#   16  nem_Malay_stampy_sorted
#   17  nem_Ngasang_stampy_sorted
#   18  nem_PM664_stampy_sorted
#   19  nem_PM665_stampy_sorted
#   20  nem_Sukai_male_stampy_sorted
#   21  nem_pagensis_stampy_sorted
#   22  nigra_PF1001_stampy_sorted
#   23  nigra_PF660_stampy_sorted
#   24  nigrescens_PM1000_stampy_sorted
#   25  nigra_PM1003_stampy_sorted
#   26  nigrescens_PF654_stampy_sorted
#   27  ochreata_PF625_stampy_sorted
#   28  ochreata_PM571_stampy_sorted
#   29  ochreata_PM596_stampy_sorted
#   30  togeanus_PF549_stampy_sorted
#   31  togeanus_PM545_stampy_sorted
#   32  tonk_PF515_stampy_sorted
#   33  tonk_PM561_stampy_sorted
#   34  tonk_PM565_stampy_sorted
#   35  tonk_PM566_stampy_sorted
#   36  tonk_PM567_stampy_sorted
#   37  tonk_PM582_stampy_sorted
#   38  tonk_PM584_stampy_sorted
#   39  tonk_PM592_stampy_sorted
#   40  tonk_PM602_stampy_sorted

# tonk
# 32_33_34_35_36_37_38_39_40 

# hecki
# 2_3_4_5_6_7  

# maura
# 8_9_10_11_12_13 

# nigra 
# 22_23_25 

# nigrescens 
# 24_26 

# ochreata 
# 27_28_29 

# togeanus 
# 30_31 

# brunn 
# 1 

# borneo 
# 18_19_20_21 

# sumatra 
# 15_17 

# pagensis 
# 21 

# malay 
# 16 


my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
my $outputfile = $ARGV[3];

unless (open DATAINPUT, $inputfile) {
    print "Can not find the input file.\n";
    exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
    print "I can\'t write to $outputfile\n";
    exit;
}

print "Creating output file: $outputfile\n";

my @sexes = split("",$ARGV[1]);
my @whotoinclude = split("_",$ARGV[2]);

my $number_of_individuals_genotyped=($#whotoinclude - 1);

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped=0;
my $y;


for ($y = 2 ; $y <= $#whotoinclude ; $y++ ) {
    if($sexes[$whotoinclude[$y]-1] == 1){
        $number_of_female_individuals_genotyped +=1;
    }   
}   


print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";


my @temp;
my @temp1;
my $string;
my $x;
my @unique;
my %ahash=();
my %xhash=();
my $type1=0;
my $type2=0;

my $w;


# initialize the ahash 
for ($y = 2 ; $y <= ($number_of_individuals_genotyped*2) ; $y+=2 ) {
    $ahash{$y}[0]=$y;
    for ($x = 1 ; $x <= ($y+1) ; $x++ ) {
        $ahash{$y}[$x]=0;
    }
}
# initialize the xhash 
for ($y = 2 ; $y <= ($number_of_individuals_genotyped*2 - ($number_of_individuals_genotyped - $number_of_female_individuals_genotyped)) ; $y++ ) {
    $xhash{$y}[0]=$y;
    for ($x = 1 ; $x <= ($y+1) ; $x++ ) {
        $xhash{$y}[$x]=0;
    }
}


# Read in datainput file
while ( my $line = <DATAINPUT>) {
    chomp($line);
    @temp=split(/[\/'\t']+/,$line);
    if($temp[0] ne '#CHROM'){
        if(($temp[0] ne "chrX")&&($temp[0] ne "chrY")&&($temp[0] ne "chrM")&&(length($temp[2]) == 1)){
            # load the autosomal data
            $string=();
            for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
                # load the first allele
                if(
                    ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '.')&&
                    ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '.')&&
                    ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '*')&&
                    ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '*') &&
                    ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '<NON_REF>')&&
                    ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '<NON_REF>')
                    )
                    {
                    $w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
                    $string=$string.$w;
                    $w = uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]];
                    $string=$string.$w;
                }   
            }
            if(defined($string)){
                @temp1=split('',$string);
                if($#temp1 >0){                 
                    # count up the number of each class of nucleotide.
                    # the first class is G or C
                    # the second class is A or T
                    # so if there is a G and a C segregating, for example, this is a non-polymorphic site
                    $type1=0;
                    $type2=0;

                    for ($y = 0 ; $y <= $#temp1 ; $y++ ) {
                        if((uc $temp1[$y] eq "G")||(uc $temp1[$y] eq "C")){
                            $type1+=1;
                        }
                        elsif((uc $temp1[$y] eq "A")||(uc $temp1[$y] eq "T")){
                            $type2+=1;
                        }
                        else{
                            print "Problemo1! ",$line,"g\n";
                        }
                    }
                    $ahash{($#temp1+1)}[$type1+1]+=1;
                }       
            }   
        }
        elsif(($temp[0] eq "chrX")&&(length($temp[2]) == 1)){
         # load the chrX data           
            $string=();
            # Need to assess whether the male genotypes have two "A/A" or one "A/" allele
            # if one, we need to modify @temp
            #print "before temp array @temp ttttt\n";
            if($#temp != (($#sexes+1)*2+$whotoinclude[1]-2)){
                # males have only one allele
                my $counterramma=0;
                foreach(@sexes){
                    if($_ == 0){
                        # insert strategically into @temp
                        splice( @temp, ($whotoinclude[1]+$counterramma), 0, "." );
                    }
                    $counterramma+=2;
                }
            }
            for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
                # load both alleles if the individual is a female
                if($sexes[$whotoinclude[$y+2]-1] eq "1"){
                    # load both alleles
                    if(
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '.')&&
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '.')&&
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '*')&&
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '*')&&
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '<NON_REF>')&&
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '<NON_REF>')
                        )
                        {
                        $w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
                        $string=$string.$w;
                        $w = uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]];
                        $string=$string.$w;
                    }   
                }
                # load one allele if the individual is a male
                elsif($sexes[$whotoinclude[$y+2]-1] eq "0"){
                    # load only the first allele
                    if(
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '.')&&
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '*')&&
                        ($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '<NON_REF>')
                        )
                        {
                        $w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
                        $string=$string.$w;
                    }   
                }
                else{
                    print "Something is wrong with figuring out what sex each individual is X ",$sexes[$whotoinclude[$y+2]-1]," ",$whotoinclude[$y+2],"\n";
                }
            } # end of cycling through each individual for chrX 
            if(defined($string)){
                @temp1=split('',$string);
                if($#temp1 >0){                 
                    # count up the number of each class of nucleotide.
                    # the first class is G or C
                    # the second class is A or T
                    # so if there is a G and a C segregating, for example, this is a non-polymorphic site
                    $type1=0;
                    $type2=0;

                    for ($y = 0 ; $y <= $#temp1 ; $y++ ) {
                        if((uc $temp1[$y] eq "G")||(uc $temp1[$y] eq "C")){
                            $type1+=1;
                        }
                        elsif((uc $temp1[$y] eq "A")||(uc $temp1[$y] eq "T")){
                            $type2+=1;
                        }
                        else{
                            print "Problemo2!\n";
                        }
                    }
                    $xhash{($#temp1+1)}[$type1+1]+=1;
                }       
            }   
        }
    } # endif to check for first line
} # end while

close DATAINPUT;


print OUTFILE "X\n";
print OUTFILE ($number_of_individuals_genotyped*2 - ($number_of_individuals_genotyped - $number_of_female_individuals_genotyped) -1),"\n";
for ($y = 2 ; $y <= ($number_of_individuals_genotyped*2 - ($number_of_individuals_genotyped - $number_of_female_individuals_genotyped)) ; $y++ ) {
    print OUTFILE $xhash{$y}[0],"\t";
    for ($x = 1 ; $x <= $y ; $x++ ) {
        print OUTFILE $xhash{$y}[$x],"\t";
    }
    print OUTFILE $xhash{$y}[$y+1],"\n";
}

print OUTFILE "A\n";
print OUTFILE $number_of_individuals_genotyped,"\n";
for ($y = 2 ; $y <= $number_of_individuals_genotyped*2 ; $y+=2 ) {
    print OUTFILE $ahash{$y}[0],"\t";
    for ($x = 1 ; $x <= $y ; $x++ ) {
        print OUTFILE $ahash{$y}[$x],"\t";
    }
    print OUTFILE $ahash{$y}[$y+1],"\n";
}


close OUTFILE;
print "Done with input file 1\n";






sub in_array {
my ($arr,$search_for) = @_;
foreach my $value (@$arr) {
    return 1 if $value eq $search_for;
}
    return 0;
}


```

We then run Kai's program using different models 100-500 times to ensure convergence on the maximum likelihood parameter values.  His program is run like this:

```
#!/bin/bash
for i in {501..1000}
do
#sqsub -r 7d -q serial -o myout ../data_2a_step_xa_m1 control_file_3_epoch $i
    /work/ben/2013.09.16_kai_program/data_2a_step_xa_m1 control_file_3_epoch $i
done
```

and uses a control file that looks like this:
```
outputFile: ../3epoch.txt
dataFile: ../../tonkeana_kai_input.txt
K: 200
nstep: 2
maxTA: 0.5 0.5
tau: 0.01
useNrSimplex: 0
nlopt_alg: NLOPT_LN_NELDERMEAD
initThetaRange:	1e-10	0.5
initGammaRange:	-50	50
initLambdaRange:	0.01	100
initRhoRange:		0.01	100
seed:
thetaOnLn: 1
lambdaOnLn: 1
rhoOnLn: 1
setBound: 1
rftol: 1e-15
maxeval: 100000
maxtime: 600
imprftol: 1e-15
nnoimp: 3
maximp: 100

```


Results can be summarized with Kai's java program 'Parsefolder' as follows:

```
java ParseFolder /work/ben/new_kai_program_2013.10.13/2015_Sularad/tonkeana/epoch3_full 1e-15 true
```

Occasionally there are problems with the outputfile because of an '-inf' value for the likelihood.  This can be solved by going to the directory and typing this:

```
sed -i -n '/inf/!p' *
```
which deletes any lines containing `inf`.

I automated this for each species using this script:

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

# this program will execute lots of commando files 
# that execute sqjobs commands.


my @directories=("equilibrium","epoch2_full","epoch2_GAMMA_A_EQ_C_0","epoch2_GAMMA_X_EQ_C_0","epoch2_GAMMA_X_EQ_LAMBDA_GAMMA_A","epoch2_lambda_eq_0.75","epoch2_THETA_01_A_EQ_THETA_10_A","epoch2_THETA_01_X_EQ_LAMBDA_THETA_01_A","epoch2_THETA_01_X_EQ_THETA_10_X","epoch2_THETA_10_X_EQ_LAMBDA_THETA_10_A","epoch3_full","epoch3_GAMMA_A_EQ_C_0","epoch3_GAMMA_X_EQ_3OVER4_GAMMA_A","epoch3_GAMMA_X_EQ_C_0","epoch3_GAMMA_X_EQ_LAMBDA_GAMMA_A","epoch3_THETA_01_A_EQ_THETA_10_A","epoch3_lambda_eq_0.75","epoch3_THETA_01_X_EQ_LAMBDA_THETA_01_A","epoch3_THETA_01_X_EQ_THETA_10_X","epoch3_THETA_10_X_EQ_LAMBDA_THETA_10_A");

my $y;
my $commandline;
my $status;
my $outfile="tonkeana";

$commandline = "java ParseFolder /work/ben/new_kai_program_2013.10.13/2015_Sularad/".$outfile."/".$directories[0]." 1e-15 true > ".$outfile.".out";
print $commandline,"\n";
$status = system($commandline);

for ($y = 1 ; $y <= $#directories ; $y++ ) {
    $commandline = "java ParseFolder /work/ben/new_kai_program_2013.10.13/2015_Sularad/".$outfile."/".$directories[$y]." 1e-15 true >> ".$outfile.".out";
    print $commandline,"\n";
    $status = system($commandline);
}

```

And I wrote a script that will parse and format the output of the above script so that it is suitable to plot in latek.  This script is called "Formats_Kais_results.pl":

```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This script reformats the output of "parse_commando.pl" which is a wrapper
# for a java program that Kai wrote to extract results from many runs

# run by typing "Formats_Kais_results.pl infile outfile"

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my %datahash;
my %weighted_parameters;
my @temp;
my $y;
my @headers;
my @models = ("equilibrium","epoch2_full","epoch2_GAMMA_A_EQ_C_0","epoch2_GAMMA_X_EQ_3OVER4_GAMMA_A","epoch2_GAMMA_X_EQ_C_0","epoch2_GAMMA_X_EQ_LAMBDA_GAMMA_A","epoch2_lambda_eq_0.75","epoch2_THETA_01_A_EQ_THETA_10_A","epoch2_THETA_01_X_EQ_LAMBDA_THETA_01_A","epoch2_THETA_01_X_EQ_THETA_10_X","epoch2_THETA_10_X_EQ_LAMBDA_THETA_10_A","epoch3_full","epoch3_GAMMA_A_EQ_C_0","epoch3_GAMMA_X_EQ_3OVER4_GAMMA_A","epoch3_GAMMA_X_EQ_C_0","epoch3_GAMMA_X_EQ_LAMBDA_GAMMA_A","epoch3_THETA_01_A_EQ_THETA_10_A","epoch3_lambda_eq_0.75","epoch3_THETA_01_X_EQ_LAMBDA_THETA_01_A","epoch3_THETA_01_X_EQ_THETA_10_X","epoch3_THETA_10_X_EQ_LAMBDA_THETA_10_A");
my @free_params = ("0","3","2","1","2","2","2","2","2","2","2","5","4","3","4","4","4","4","4","4","4");
my @paramnames = ("folder_path","nconvg","theta_01_x","theta_10_x","gamma_x","theta_01_a","theta_10_a","gamma_a","lambda","rho_1","tau_a_1","rho_2","tau_a_2");
my $switch=0;

# make a hash with the free parameter values 
my %free_params;
for ($y = 0 ; $y <= $#models; $y++ ) {
	$free_params{$models[$y]}=$free_params[$y];
}

# read in the data
# and calculate the AIC according to WAGENMAKERS and FARRELL
# Psychonomic Bulletin & Review 2004, 11 (1), 192-196

my $best_AIC=0;
while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'folder_path'){
		for ($y = 1 ; $y <= $#headers; $y++ ) {
			$datahash{$temp[0].'_'.$headers[$y]} = $temp[$y];
		}
		# calculate the AIC value
		$datahash{$temp[0].'_AIC'} = -2*$datahash{$temp[0].'_lnLike'}+2*$free_params{$temp[0]};
		# and calculate the minimum AIC value
		if($best_AIC == 0){
			$best_AIC = $datahash{$temp[0].'_AIC'};
		}	
		elsif($datahash{$temp[0].'_AIC'} < $best_AIC){
			$best_AIC = $datahash{$temp[0].'_AIC'};
		}
	}
	else{
		@headers=@temp;
	}
}
#print "best AIC is ",$best_AIC,"\n";
my $sum_e_raised_to_neg5_times_deltaAIC=0;

# now calculate deltaAIC and the numbers for the AIC weights
foreach(@models){
	if(exists($datahash{$_.'_AIC'})){
		$datahash{$_.'_deltaAIC'} = $datahash{$_.'_AIC'} - $best_AIC;
		$datahash{$_.'_e_raised_to_neg5_times_deltaAIC'} = exp(-0.5*$datahash{$_.'_deltaAIC'});
		$sum_e_raised_to_neg5_times_deltaAIC+=$datahash{$_.'_e_raised_to_neg5_times_deltaAIC'};
	}
	else{
		print "This model does not have a defined lnL ",$_,"\n";
	}
}

# now calculate the wi_deltaAIC (the AIC weights)
foreach(@models){
	if(exists($datahash{$_.'_AIC'})){
		$datahash{$_.'_wi_deltaAIC'} = $datahash{$_.'_e_raised_to_neg5_times_deltaAIC'}/$sum_e_raised_to_neg5_times_deltaAIC;
	}
	else{
		print "This model does not have a defined lnL ",$_,"\n";
	}
}	

# now print the results
# first print headers
print OUTFILE "model	thetaxa	thetaxb	gammax	thetaaa	thetaab	gammaa	lambda	rhoa	taua	rhob	taub	lnL	deltaAIC	wiAIC\n";
print "model\t";
#print OUTFILE "model\t";
for ($y = 2 ; $y <= $#paramnames; $y++ ) {
	print $paramnames[$y],"\t";
	#print OUTFILE $paramnames[$y],"\t";
}
# proviously I printed out all of the AIC stuff
#print "lnL\tAIC\tdeltaAIC\te^-0.5deltaAIC\twi_AIC\n";
#print OUTFILE "lnL\tAIC\tdeltaAIC\te^-0.5deltaAIC\twi_AIC\n";
# now I want to print out only deltaAIC and w_AIC
print "lnL\tdeltaAIC\twi_AIC\n";
#print OUTFILE "lnL\tdeltaAIC\twi_AIC\n";

print OUTFILE "1 Epoch\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\n";

# calculate the AIC weights and print out the formatted values
foreach(@models){
	if($_ eq 'epoch2_full'){
		print OUTFILE "2 Epochs\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\n";
	}
	elsif($_ eq 'epoch3_full'){
		print OUTFILE "3 Epochs\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\n";
	}
	print $_;
	print OUTFILE $_;

	if(exists($datahash{$_.'nconvg'})){
		if($datahash{$_.'nconvg'} eq 1){
			print "*\t";
			print OUTFILE "\t";
			$switch=1;
		}
	}	
	else{
		print "\t";
		print OUTFILE "\t";
	}
	for ($y = 2 ; $y <= $#paramnames; $y++ ) {
		if(exists($datahash{$_.'_'.$paramnames[$y]})){
			# print out thetas and gammas to 5 decimal points
			if(($paramnames[$y] eq 'theta_01_x')||($paramnames[$y] eq 'theta_10_x')||($paramnames[$y] eq 'gamma_x')||($paramnames[$y] eq 'theta_01_a')||($paramnames[$y] eq 'theta_10_a')||($paramnames[$y] eq 'gamma_a')){
				print sprintf("%.5f",$datahash{$_.'_'.$paramnames[$y]})," \t";
				print OUTFILE sprintf("%.5f",$datahash{$_.'_'.$paramnames[$y]})," \t";
				$weighted_parameters{$paramnames[$y]}+=$datahash{$_.'_'.$paramnames[$y]}*$datahash{$_.'_wi_deltaAIC'};
			}
			else{
				print sprintf("%.3f",$datahash{$_.'_'.$paramnames[$y]})," \t";
				print OUTFILE sprintf("%.3f",$datahash{$_.'_'.$paramnames[$y]})," \t";
				$weighted_parameters{$paramnames[$y]}+=$datahash{$_.'_'.$paramnames[$y]}*$datahash{$_.'_wi_deltaAIC'};
			}	
		}
		elsif($_ eq 'equilibrium'){
				print "- \t";
				print OUTFILE "- \t";
		}
		elsif($paramnames[$y] eq 'theta_01_x'){
			if(($_ eq 'epoch2_THETA_01_X_EQ_LAMBDA_THETA_01_A')||($_ eq 'epoch3_THETA_01_X_EQ_LAMBDA_THETA_01_A')){
				print "lq01A \t";
				print OUTFILE "lq01A \t";
				$weighted_parameters{'theta_01_x'}+=$datahash{$_.'_'.'lambda'}*$datahash{$_.'_'.'theta_01_a'}*$datahash{$_.'_wi_deltaAIC'};
			}
			if(($_ eq 'epoch2_THETA_01_X_EQ_THETA_10_X')||($_ eq 'epoch3_THETA_01_X_EQ_THETA_10_X')){
				print "q10x \t";
				print OUTFILE "q10x \t";
				$weighted_parameters{'theta_01_x'}+=$datahash{$_.'_'.'theta_10_x'}*$datahash{$_.'_wi_deltaAIC'};
			}
		}
		elsif($paramnames[$y] eq 'theta_10_x'){
			if(($_ eq 'epoch2_THETA_10_X_EQ_LAMBDA_THETA_10_A')||($_ eq 'epoch3_THETA_10_X_EQ_LAMBDA_THETA_10_A')){
				print "lq10a \t";
				print OUTFILE "lq10a \t";
				$weighted_parameters{'theta_10_x'}+=$datahash{$_.'_'.'lambda'}*$datahash{$_.'_'.'theta_10_a'}*$datahash{$_.'_wi_deltaAIC'};
			}
		}
		elsif($paramnames[$y] eq 'theta_01_a'){
			if(($_ eq 'epoch2_THETA_01_A_EQ_THETA_10_A')||($_ eq 'epoch3_THETA_01_A_EQ_THETA_10_A')){
				print "q10a \t";
				print OUTFILE "q10a \t";
				$weighted_parameters{'theta_01_a'}+=$datahash{$_.'_'.'theta_10_a'}*$datahash{$_.'_wi_deltaAIC'};
			}
		}
		elsif($paramnames[$y] eq 'rho_2'){
			print "- \t";
			print OUTFILE "- \t";
		}
		elsif($paramnames[$y] eq 'tau_a_2'){
			print "- \t";
			print OUTFILE "- \t";
		}		
		elsif($paramnames[$y] eq 'gamma_a'){
			print "0(fixed) \t";
			print OUTFILE "0(fixed) \t";
		}		
		elsif($paramnames[$y] eq 'gamma_x'){
			if(($_ eq 'epoch2_GAMMA_X_EQ_C_0')||($_ eq 'epoch3_GAMMA_X_EQ_C_0')){
				print "0(fixed) \t";
				print OUTFILE "0(fixed) \t";
				# no need to add to the weighted parameter value
			}
			elsif(($_ eq 'epoch2_GAMMA_X_EQ_3OVER4_GAMMA_A')||($_ eq 'epoch3_GAMMA_X_EQ_3OVER4_GAMMA_A')){
				print "0.75gA \t";
				print OUTFILE "0.75gA \t";
				$weighted_parameters{'gamma_x'}+=0.75*$datahash{$_.'_'.'gamma_a'}*$datahash{$_.'_wi_deltaAIC'};
			}	
			elsif(($_ eq 'epoch2_GAMMA_X_EQ_LAMBDA_GAMMA_A')||($_ eq 'epoch3_GAMMA_X_EQ_LAMBDA_GAMMA_A')){
				print "lgA \t";
				print OUTFILE "lgA \t";
				$weighted_parameters{'gamma_x'}+=$datahash{$_.'_'.'lambda'}*$datahash{$_.'_'.'gamma_a'}*$datahash{$_.'_wi_deltaAIC'};
			}	
		}
		elsif($paramnames[$y] eq 'lambda'){
				print "0.75(fixed) \t";
				print OUTFILE "0.75(fixed) \t";
				$weighted_parameters{'lambda'}+=0.75*$datahash{$_.'_wi_deltaAIC'};	
		}				
		else{
			print "fixme \t";
			print OUTFILE "fixme \t";
		}
	}	
	if(exists($datahash{$_.'_AIC'})){
		$datahash{$_.'_wi_deltaAIC'} = $datahash{$_.'_e_raised_to_neg5_times_deltaAIC'}/$sum_e_raised_to_neg5_times_deltaAIC;
		# proviously I printed out all of the AIC stuff
		# print sprintf("%.3f",$datahash{$_.'_lnLike'})," \t",sprintf("%.3f",$datahash{$_.'_AIC'})," \t",sprintf("%.3f",$datahash{$_.'_deltaAIC'})," \t ",sprintf("%.3f",$datahash{$_.'_e_raised_to_neg5_times_deltaAIC'})," \t",sprintf("%.3f",$datahash{$_.'_wi_deltaAIC'}),"\n";
		# print OUTFILE sprintf("%.3f",$datahash{$_.'_lnLike'})," \t",sprintf("%.3f",$datahash{$_.'_AIC'})," \t",sprintf("%.3f",$datahash{$_.'_deltaAIC'})," \t ",sprintf("%.3f",$datahash{$_.'_e_raised_to_neg5_times_deltaAIC'})," \t",sprintf("%.3f",$datahash{$_.'_wi_deltaAIC'}),"\n";
		# now I want to print out only deltaAIC and w_AIC
		print sprintf("%.3f",$datahash{$_.'_lnLike'})," \t",sprintf("%.2f",$datahash{$_.'_deltaAIC'})," \t ",sprintf("%.3f",$datahash{$_.'_wi_deltaAIC'}),"\n";
		print OUTFILE sprintf("%.3f",$datahash{$_.'_lnLike'})," \t",sprintf("%.2f",$datahash{$_.'_deltaAIC'})," \t ",sprintf("%.3f",$datahash{$_.'_wi_deltaAIC'}),"\n";
	}
	else{
		print "This model does not have a defined lnL ",$_,"\n";
		print OUTFILE "This model does not have a defined lnL ",$_,"\n";
	}
}

	if($switch == 1){
		print "* indicates convergence was achieved with rtol > 1e-15\n";
		print OUTFILE "* indicates convergence was achieved with rtol > 1e-15\n";
	}
close DATAINPUT;


# print out the weighted parameter values
print OUTFILE "NAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\n";

		print "Model_average\t";
		print OUTFILE "Model_average\t";

for ($y = 2 ; $y < 8; $y++ ) {
	if(exists($weighted_parameters{$paramnames[$y]})){
		print $weighted_parameters{$paramnames[$y]},"\t";
		#if($weighted_parameters{$paramnames[$y]} =~ /^[0-9,.E]+$/ ){
			print OUTFILE sprintf("%.5f",$weighted_parameters{$paramnames[$y]}),"\t";
		#}
		#else{
		#	print OUTFILE $weighted_parameters{$paramnames[$y]},"\t";
		#}	
	}
	else{
		print "0\t";
		print OUTFILE "0\t";
	}	
}
for ($y = 8 ; $y <= $#paramnames; $y++ ) {
	if(exists($weighted_parameters{$paramnames[$y]})){
		print $weighted_parameters{$paramnames[$y]},"\t";
		#if($weighted_parameters{$paramnames[$y]} =~ /^[0-9,.E]+$/ ){
			print OUTFILE sprintf("%.3f",$weighted_parameters{$paramnames[$y]}),"\t";
		#}
		#else{
		#	print OUTFILE $weighted_parameters{$paramnames[$y]},"\t";
		#}	
	}
	else{
		print "0\t";
		print OUTFILE "0\t";
	}	
}

print "NAN\tNAN\tNAN\n";
print OUTFILE "NAN\tNAN\tNAN\n";

close OUTFILE;

sub in_array {
my ($arr,$search_for) = @_;
foreach my $value (@$arr) {
	return 1 if $value eq $search_for;
}
 	return 0;
}



```

# Bootstrap analysis

First I made a directory called `epoch2_full_bootstrap`.  Then I ran this script: `Generates_lots_of_control_and_commando_files_for_bootstraps.pl`

In order to run the next of the script I needed to load the Math::Random module on sharcnet like this:
```
module load cpan/2.10
```
Then I could run this script:`Generates_bootstraps_of_tab_delimited_data_for_kai_program.pl`.  Then I modified the `directory_commando` script to make the commando files executable.  ANd then I changed it again to make it launch the runs on sharcnet.  Then (on iqaluk because it does not work on kraken) I ran the parsefolder_commando script after making sure the `Parsefolder_class` program was in the `epoch2_full_bootstrap` directory.
