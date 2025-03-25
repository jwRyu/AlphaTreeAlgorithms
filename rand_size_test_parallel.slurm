#!/bin/bash
#SBATCH --job-name=atree_test
#SBATCH --ntasks=1         
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --partition=regular

module load OpenCV/4.6.0-foss-2022a-contrib

./AlphaTree config/config_pm_8bit_4N_1000.txt | tee -a output/out_size_pm_8bit_4N.txt
./AlphaTree config/config_pm_8bit_4N_4000.txt | tee -a output/out_size_pm_8bit_4N.txt
./AlphaTree config/config_pm_8bit_4N_6000.txt | tee -a output/out_size_pm_8bit_4N.txt
./AlphaTree config/config_pm_8bit_4N_8000.txt | tee -a output/out_size_pm_8bit_4N.txt
./AlphaTree config/config_pm_8bit_4N_10000.txt | tee -a output/out_size_pm_8bit_4N.txt
./AlphaTree config/config_pm_8bit_4N_12000.txt | tee -a output/out_size_pm_8bit_4N.txt

./AlphaTree config/config_pm_8bit_8N_1000.txt | tee -a output/out_size_pm_8bit_8N.txt
./AlphaTree config/config_pm_8bit_8N_4000.txt | tee -a output/out_size_pm_8bit_8N.txt
./AlphaTree config/config_pm_8bit_8N_6000.txt | tee -a output/out_size_pm_8bit_8N.txt
./AlphaTree config/config_pm_8bit_8N_8000.txt | tee -a output/out_size_pm_8bit_8N.txt
./AlphaTree config/config_pm_8bit_8N_10000.txt | tee -a output/out_size_pm_8bit_8N.txt
./AlphaTree config/config_pm_8bit_8N_12000.txt | tee -a output/out_size_pm_8bit_8N.txt

./AlphaTree config/config_pm_16bit_4N_1000.txt | tee -a output/out_size_pm_16bit_4N.txt
./AlphaTree config/config_pm_16bit_4N_4000.txt | tee -a output/out_size_pm_16bit_4N.txt
./AlphaTree config/config_pm_16bit_4N_6000.txt | tee -a output/out_size_pm_16bit_4N.txt
./AlphaTree config/config_pm_16bit_4N_8000.txt | tee -a output/out_size_pm_16bit_4N.txt
./AlphaTree config/config_pm_16bit_4N_10000.txt | tee -a output/out_size_pm_16bit_4N.txt
./AlphaTree config/config_pm_16bit_4N_12000.txt | tee -a output/out_size_pm_16bit_4N.txt

./AlphaTree config/config_pm_16bit_8N_1000.txt | tee -a output/out_size_pm_16bit_8N.txt
./AlphaTree config/config_pm_16bit_8N_4000.txt | tee -a output/out_size_pm_16bit_8N.txt
./AlphaTree config/config_pm_16bit_8N_6000.txt | tee -a output/out_size_pm_16bit_8N.txt
./AlphaTree config/config_pm_16bit_8N_8000.txt | tee -a output/out_size_pm_16bit_8N.txt
./AlphaTree config/config_pm_16bit_8N_10000.txt | tee -a output/out_size_pm_16bit_8N.txt
./AlphaTree config/config_pm_16bit_8N_12000.txt | tee -a output/out_size_pm_16bit_8N.txt

./AlphaTree config/config_hybrid_8bit_4N_1000.txt | tee -a output/out_size_hybrid_8bit_4N.txt
./AlphaTree config/config_hybrid_8bit_4N_4000.txt | tee -a output/out_size_hybrid_8bit_4N.txt
./AlphaTree config/config_hybrid_8bit_4N_6000.txt | tee -a output/out_size_hybrid_8bit_4N.txt
./AlphaTree config/config_hybrid_8bit_4N_8000.txt | tee -a output/out_size_hybrid_8bit_4N.txt
./AlphaTree config/config_hybrid_8bit_4N_10000.txt | tee -a output/out_size_hybrid_8bit_4N.txt
./AlphaTree config/config_hybrid_8bit_4N_12000.txt | tee -a output/out_size_hybrid_8bit_4N.txt

./AlphaTree config/config_hybrid_8bit_8N_1000.txt | tee -a output/out_size_hybrid_8bit_8N.txt
./AlphaTree config/config_hybrid_8bit_8N_4000.txt | tee -a output/out_size_hybrid_8bit_8N.txt
./AlphaTree config/config_hybrid_8bit_8N_6000.txt | tee -a output/out_size_hybrid_8bit_8N.txt
./AlphaTree config/config_hybrid_8bit_8N_8000.txt | tee -a output/out_size_hybrid_8bit_8N.txt
./AlphaTree config/config_hybrid_8bit_8N_10000.txt | tee -a output/out_size_hybrid_8bit_8N.txt
./AlphaTree config/config_hybrid_8bit_8N_12000.txt | tee -a output/out_size_hybrid_8bit_8N.txt

./AlphaTree config/config_hybrid_16bit_4N_1000.txt | tee -a output/out_size_hybrid_16bit_4N.txt
./AlphaTree config/config_hybrid_16bit_4N_4000.txt | tee -a output/out_size_hybrid_16bit_4N.txt
./AlphaTree config/config_hybrid_16bit_4N_6000.txt | tee -a output/out_size_hybrid_16bit_4N.txt
./AlphaTree config/config_hybrid_16bit_4N_8000.txt | tee -a output/out_size_hybrid_16bit_4N.txt
./AlphaTree config/config_hybrid_16bit_4N_10000.txt | tee -a output/out_size_hybrid_16bit_4N.txt
./AlphaTree config/config_hybrid_16bit_4N_12000.txt | tee -a output/out_size_hybrid_16bit_4N.txt

./AlphaTree config/config_hybrid_16bit_8N_1000.txt | tee -a output/out_size_hybrid_16bit_8N.txt
./AlphaTree config/config_hybrid_16bit_8N_4000.txt | tee -a output/out_size_hybrid_16bit_8N.txt
./AlphaTree config/config_hybrid_16bit_8N_6000.txt | tee -a output/out_size_hybrid_16bit_8N.txt
./AlphaTree config/config_hybrid_16bit_8N_8000.txt | tee -a output/out_size_hybrid_16bit_8N.txt
./AlphaTree config/config_hybrid_16bit_8N_10000.txt | tee -a output/out_size_hybrid_16bit_8N.txt
./AlphaTree config/config_hybrid_16bit_8N_12000.txt | tee -a output/out_size_hybrid_16bit_8N.txt

./AlphaTree config/config_hybrid_32bit_4N_1000.txt | tee -a output/out_size_hybrid_32bit_4N.txt
./AlphaTree config/config_hybrid_32bit_4N_4000.txt | tee -a output/out_size_hybrid_32bit_4N.txt
./AlphaTree config/config_hybrid_32bit_4N_6000.txt | tee -a output/out_size_hybrid_32bit_4N.txt
./AlphaTree config/config_hybrid_32bit_4N_8000.txt | tee -a output/out_size_hybrid_32bit_4N.txt
./AlphaTree config/config_hybrid_32bit_4N_10000.txt | tee -a output/out_size_hybrid_32bit_4N.txt
./AlphaTree config/config_hybrid_32bit_4N_12000.txt | tee -a output/out_size_hybrid_32bit_4N.txt

./AlphaTree config/config_hybrid_32bit_8N_1000.txt | tee -a output/out_size_hybrid_32bit_8N.txt
./AlphaTree config/config_hybrid_32bit_8N_4000.txt | tee -a output/out_size_hybrid_32bit_8N.txt
./AlphaTree config/config_hybrid_32bit_8N_6000.txt | tee -a output/out_size_hybrid_32bit_8N.txt
./AlphaTree config/config_hybrid_32bit_8N_8000.txt | tee -a output/out_size_hybrid_32bit_8N.txt
./AlphaTree config/config_hybrid_32bit_8N_10000.txt | tee -a output/out_size_hybrid_32bit_8N.txt
./AlphaTree config/config_hybrid_32bit_8N_12000.txt | tee -a output/out_size_hybrid_32bit_8N.txt

./AlphaTree config/config_hybrid_64bit_4N_1000.txt | tee -a output/out_size_hybrid_64bit_4N.txt
./AlphaTree config/config_hybrid_64bit_4N_4000.txt | tee -a output/out_size_hybrid_64bit_4N.txt
./AlphaTree config/config_hybrid_64bit_4N_6000.txt | tee -a output/out_size_hybrid_64bit_4N.txt
./AlphaTree config/config_hybrid_64bit_4N_8000.txt | tee -a output/out_size_hybrid_64bit_4N.txt
./AlphaTree config/config_hybrid_64bit_4N_10000.txt | tee -a output/out_size_hybrid_64bit_4N.txt
./AlphaTree config/config_hybrid_64bit_4N_12000.txt | tee -a output/out_size_hybrid_64bit_4N.txt

./AlphaTree config/config_hybrid_64bit_8N_1000.txt | tee -a output/out_size_hybrid_64bit_8N.txt
./AlphaTree config/config_hybrid_64bit_8N_4000.txt | tee -a output/out_size_hybrid_64bit_8N.txt
./AlphaTree config/config_hybrid_64bit_8N_6000.txt | tee -a output/out_size_hybrid_64bit_8N.txt
./AlphaTree config/config_hybrid_64bit_8N_8000.txt | tee -a output/out_size_hybrid_64bit_8N.txt
./AlphaTree config/config_hybrid_64bit_8N_10000.txt | tee -a output/out_size_hybrid_64bit_8N.txt
./AlphaTree config/config_hybrid_64bit_8N_12000.txt | tee -a output/out_size_hybrid_64bit_8N.txt

./AlphaTree config/config_pm_8bit_4N_20000.txt
./AlphaTree config/config_pm_8bit_4N_25000.txt
./AlphaTree config/config_pm_8bit_4N_31640.txt
./AlphaTree config/config_pm_8bit_8N_20000.txt
./AlphaTree config/config_pm_8bit_8N_25000.txt
./AlphaTree config/config_pm_8bit_8N_31640.txt