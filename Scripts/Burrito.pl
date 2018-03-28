##########################
#
# Author: Sariel Hubner, sariel.hubner@botany.ubc.ca
# Ammended by: Rishi R. Masalia, rishimasalia@gmail.com, Jan 2016 - April 2017, Burke Lab, UGA
# Please don't distribute
#
# Summary: Burrito is a wrapper for the Manhattanizer (augmentation of Sariel's original Helimap pipeline)
# It treats traits and environments individually, and produces a Manhattan plot & list of significant SNPs.
# EMMAX is run with kinship and structure (PCA) to get associations
#
# Preferences are controlled in the configuration file ("Asada.config")
##########################

$file1 = $ARGV[0]; #list of environments argument
$file2 = $ARGV[1]; #list of traits argument

open(ENVS,$file1);

while(<ENVS>){
	chomp();	
	$Env = $_;
	print "$Env\n";

 	# get the pheno file:
	system("cp ../Phenos/*\_$Env.txt ../EnchiladaSuite/Burrito/TraitHopper/");

	open(TRAITS,$file2);
	while(<TRAITS>){
		print "$_\n";
		open(CONFIG, ">Asada.config");
		chomp();
		$comment = "\# Config File";
		$VC = "NA";
		$group = "XRQ_213.group";
		$VCF = "XRQ_213.vcf";
		$Trait = $_;
		$FullRun = "no";
		$Num_env = 1;
		$Num_rep= 1;
		$Prefix= "XRQ_213";
		$env1_Name= $Env;
		$env1_file = "$Trait\_$env1_Name\.txt"; 
		$env2_Name = "NA";
        	$env2_file = "NA";
		$env3_Name = "NA";
        	$env3_file = "NA";
		$env4_Name = "NA";
        	$env4_file = "NA";
		$env5_Name = "NA";
        	$env5_file = "NA";
		$kinship = "no";
		$PCA = "no";
		$fdrcutoff = 0.05;
		$color1 = "dodgerblue4"; #### Choose whatever colors you want for 1 and 2 from R colors
		$color2 = "red";

		print CONFIG "$comment\n$VC\n$group\n$VCF\n$Trait\n$FullRun\n$Num_env\n$Num_rep\n$Prefix\n$env1_Name\n$env1_file\n$env2_Name\n$env2_file\n$env3_Name\n$env3_file\n$env4_Name\n$env4_file\n$env5_Name\n$env5_file\n$kinship\n$PCA\n$fdrcutoff\n$color1\n$color2\n";

		system("mv Asada.config ../EnchiladaSuite/Burrito/");
		system("cp ../EnchiladaSuite/Burrito/TraitHopper/$Trait\_* ../EnchiladaSuite/Burrito/");
		system("Rscript ./XRQ_213Manhat.R ../EnchiladaSuite/Burrito/Asada.config");	
		system("mv ../EnchiladaSuite/Burrito/*.ps ../Outputs/Tables/PS_Files/");
		system("mv ../EnchiladaSuite/Burrito/*.sigsnps ../Outputs/Tables/SigSNPs/");
		system("mv ../EnchiladaSuite/Burrito/*.pdf ../Outputs/Plots/Manhattan_Plots/");
	} close TRAITS;

	#put the data into a folder:

	system("rm ../EnchiladaSuite/Burrito/TraitHopper/*");
	

} close ENVS; 







