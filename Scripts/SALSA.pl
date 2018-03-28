##########################
# Author: Rishi R. Masalia, rishimasalia@gmail.com, UGA, Burke Lab, April 2017
#
# Purpose: this pipeline is a direct follow up to the GWAS Pipeline by Hubner and Masalia. 
# It is designed to take significnat SNPs from a manhattan plot, and identify LD boundaries
# by looking at pairwise LD across significant SNPs on a chromosome
#
# File from where to grab the phenos
#
# Run perl SALSA.pl environments.list traits.list
########################################

$file1 = $ARGV[0]; #list of environments argument
$file2 = $ARGV[1]; #list of traits argument

open(ENVS,$file1);

while(<ENVS>){
	chomp();	
	$Env = $_;
	print "$Env\n";
	
	# get the pheno file:
	system("cp ../Outputs/Tables/SigSNPs/*\_$Env.sigsnps ../EnchiladaSuite/SALSA/TraitHopper/");
	

	open(TRAITS,$file2);
	while(<TRAITS>){
		open(CONFIG, ">Tomato.config");
		chomp();
		$comment = "\# Config File for Determining LD and SigSNP overlap";
		$Trait = $_;
		$Num_env = 1;
		$Prefix= "XRQ_213_Salsa";
		$env1_Name= "$Env";
		$env1_file = "$Trait\_$Env\.sigsnps"; 
		$env2_Name = "NA";
        	$env2_file = "NA";
		$env3_Name = "NA";
        	$env3_file = "NA";
		$env4_Name = "NA";
        	$env4_file = "NA";
		$env5_Name = "NA";
        	$env5_file = "NA";
		$thresh = 'bon'; #'bon' for bonferroni / simpleM 'qfdr' for FDR

		print CONFIG "$comment\n$Trait\n$Num_env\n$Prefix\n$env1_Name\n$env1_file\n$env2_Name\n$env2_file\n$env3_Name\n$env3_file\n$env4_Name\n$env4_file\n$env5_Name\n$env5_file\n$thresh\n";

		system("mv Tomato.config ../EnchiladaSuite/SALSA/");
		system("cp ../EnchiladaSuite/SALSA/TraitHopper/$Trait\_* ../EnchiladaSuite/SALSA/");
		system("Rscript ./SALSA_New.R ../EnchiladaSuite/SALSA/Tomato.config");	
		system("rm ../EnchiladaSuite/SALSA/*.sigsnps");
		system("rm ../EnchiladaSuite/SALSA/Ready_tmp.tped");
		system("rm ../EnchiladaSuite/SALSA/tmp.tped");
		system("rm ../EnchiladaSuite/SALSA/tmp.bins");
		
	} close TRAITS;

	#put the data into a folder:

	system("rm ../EnchiladaSuite/SALSA/TraitHopper/*");
	
} close ENVS; 

system ("cat ../EnchiladaSuite/SALSA/PICs/* > ../Outputs/Tables/Pre_Identified_Clusters.csv");

 




