# To reproduce Figure 3
#!/bin/bash
DGP_name="Separate_Y_Dep_pi_POW_VaryAlterProp" #DataGenerating ID

#mkdir -p figs
mkdir -p logerr


for simu_n in 250 #250
do
for simu_p in 5000 #20000 #Meier  FullRank SPAM  LowRank 
do  
for H1_type in TwoSide
do
nohup Rscript Simu_Setting/RunDGP_SepY_Dep_pi_POW_VaryAlterProp.R $DGP_name $H1_type $simu_n $simu_p  > logerr/simuSepY-POW-VaryAlterProp-$DGP_name-$H1_type-$simu_n-$simu_p.log 2>logerr/simuSepY-POW-VaryAlterProp-$DGP_name-$H1_type-$simu_n-$simu_p.err &
  wait
done
done
done


#nohup Rscript Simu_Setting/RunDGP_OC_WG_pq100.R >logerr/nohupWG100OC.out 2>logerr/nohupWG100OC.err &
#nohup Rscript Simu_Setting/RunDGP_OC_WG_pq100.R >logerr/nohupWG200OC.out 2>logerr/nohupWG200OC.err &
#nohup Rscript Simu_Setting/RunDGP_OC_NG_pq200.R >logerr/nohupNG100OC.out 2>logerr/nohupNG100OC.err &
#nohup Rscript Simu_Setting/RunDGP_OC_NG_pq200.R >logerr/nohupNG200OC.out 2>logerr/nohupNG200OC.err &
