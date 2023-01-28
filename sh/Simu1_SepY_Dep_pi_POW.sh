# To reproduce Tables S2-S4
#!/bin/bash
DGP_name="Separate_Y_Dep_pi_POW" #DataGenerating ID

#mkdir -p figs
mkdir -p logerr


for simu_n in 250 #250
do
for simu_p in 5000 #20000 #Meier  FullRank SPAM  LowRank 
do  
for H1_type in TwoSide
do
nohup Rscript Simu_Setting/RunDGP_SepY_Dep_pi_POW.R $DGP_name $H1_type $simu_n $simu_p  > logerr/simuSepY-POW-$DGP_name-$H1_type-$simu_n-$simu_p.log 2>logerr/simuSepY-POW-$DGP_name-$H1_type-$simu_n-$simu_p.err &
  wait
done
done
done

