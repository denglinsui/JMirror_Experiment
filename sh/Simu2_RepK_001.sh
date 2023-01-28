# To reproduce Figure 5
#!/bin/bash
DGP_name="Rep_K1" #DataGenerating ID

#mkdir -p figs
mkdir -p logerr


for alter_prop in 0.01 #250
do
for H1_type in TwoSide
do
nohup Rscript Simu_Setting/RunDGP_Rep.R $DGP_name $H1_type $alter_prop  > logerr/simuRep-$DGP_name-$H1_type-$alter_prop.log 2>logerr/simuRep-$DGP_name-$H1_type-$alter_prop.err &
  wait
done
done

