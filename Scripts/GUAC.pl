##############################################
# 
# Author: Rishi R. Masalia, rishimasalia@gmail.com
#
# Purpose: this script is a direct follow up to SALSA by Masalia. 
# It is designed to take a "Pre_Identified_Cluster.csv" list and create a new list of genes underlying those clusters
#
# 
# Run perl GUAC.pl
########################################

############# Make the Config File ###########################
$buff = $ARGV[0];
open(CONFIG, ">Avocado.config");
chomp();
$comment = "\# Config File for guac: identifying candidate genes";
$MainEnv = "../../Outputs/Tables/Pre_Identified_Clusters.csv";
$genes = $buff; # +/- number of genes at the edge of every cluster 

print CONFIG "$comment\n$MainEnv\n$genes\n";

############# Run The R Scripts ###########################

system("Rscript GUAC.R Avocado.config");	


