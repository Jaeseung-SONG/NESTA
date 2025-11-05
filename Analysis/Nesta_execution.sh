Phenotype='/home/js/Thyroid_disorder/TWAS_res/'

Ex_net='/home/js/Thyroid_disorder/Coex_Net_Thyr'

CT_list=("$Ex_net"/*.rds)

# CT_list=$(ls "$Ex_net" | grep '\.rds$')

for i in "$Phenotype"/*; do

  phen_name=$(basename -s .rds "$i")
  
  # for j in "$CT_list"; do
  
  for j in "${CT_list[@]}"; do

  
    CT=$(basename -s .rds "$j")
    
    echo $j;
    
    echo $CT;
    
    Anal_name="${phen_name}_${CT}"

    Rscript NESTA.R \
      --TWAS_res "$i" \
      --Reference_net "$j" \
      --Diffuse_grid TRUE \
      --out_dir './Grid_test/' \
      --prefix "$Anal_name" \
      --Analysis_name "$CT" \
      
  done
done
