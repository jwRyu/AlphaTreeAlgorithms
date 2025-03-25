#!/bin/bash
#SBATCH --job-name=atree_test
#SBATCH --ntasks=1         
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --partition=regular

module load OpenCV/4.6.0-foss-2022a-contrib

./AlphaTree config/config_hq_8bit_4N_1000.txt | tee -a output/out_size_hq_8bit_4N.txt
./AlphaTree config/config_hq_8bit_4N_2000.txt | tee -a output/out_size_hq_8bit_4N.txt
./AlphaTree config/config_hq_8bit_4N_4000.txt | tee -a output/out_size_hq_8bit_4N.txt
./AlphaTree config/config_hq_8bit_4N_6000.txt | tee -a output/out_size_hq_8bit_4N.txt
./AlphaTree config/config_hq_8bit_4N_8000.txt | tee -a output/out_size_hq_8bit_4N.txt
./AlphaTree config/config_hq_8bit_4N_10000.txt | tee -a output/out_size_hq_8bit_4N.txt

./AlphaTree config/config_hq_8bit_8N_1000.txt | tee -a output/out_size_hq_8bit_8N.txt
./AlphaTree config/config_hq_8bit_8N_2000.txt | tee -a output/out_size_hq_8bit_8N.txt
./AlphaTree config/config_hq_8bit_8N_4000.txt | tee -a output/out_size_hq_8bit_8N.txt
./AlphaTree config/config_hq_8bit_8N_6000.txt | tee -a output/out_size_hq_8bit_8N.txt
./AlphaTree config/config_hq_8bit_8N_8000.txt | tee -a output/out_size_hq_8bit_8N.txt
./AlphaTree config/config_hq_8bit_8N_10000.txt | tee -a output/out_size_hq_8bit_8N.txt

./AlphaTree config/config_hq_16bit_4N_1000.txt | tee -a output/out_size_hq_16bit_4N.txt
./AlphaTree config/config_hq_16bit_4N_2000.txt | tee -a output/out_size_hq_16bit_4N.txt
./AlphaTree config/config_hq_16bit_4N_4000.txt | tee -a output/out_size_hq_16bit_4N.txt
./AlphaTree config/config_hq_16bit_4N_6000.txt | tee -a output/out_size_hq_16bit_4N.txt
./AlphaTree config/config_hq_16bit_4N_8000.txt | tee -a output/out_size_hq_16bit_4N.txt
./AlphaTree config/config_hq_16bit_4N_10000.txt | tee -a output/out_size_hq_16bit_4N.txt

./AlphaTree config/config_hq_16bit_8N_1000.txt | tee -a output/out_size_hq_16bit_8N.txt
./AlphaTree config/config_hq_16bit_8N_2000.txt | tee -a output/out_size_hq_16bit_8N.txt
./AlphaTree config/config_hq_16bit_8N_4000.txt | tee -a output/out_size_hq_16bit_8N.txt
./AlphaTree config/config_hq_16bit_8N_6000.txt | tee -a output/out_size_hq_16bit_8N.txt
./AlphaTree config/config_hq_16bit_8N_8000.txt | tee -a output/out_size_hq_16bit_8N.txt
./AlphaTree config/config_hq_16bit_8N_10000.txt | tee -a output/out_size_hq_16bit_8N.txt

./AlphaTree config/config_hhq_8bit_4N_1000.txt | tee -a output/out_size_hhq_8bit_4N.txt
./AlphaTree config/config_hhq_8bit_4N_2000.txt | tee -a output/out_size_hhq_8bit_4N.txt
./AlphaTree config/config_hhq_8bit_4N_4000.txt | tee -a output/out_size_hhq_8bit_4N.txt
./AlphaTree config/config_hhq_8bit_4N_6000.txt | tee -a output/out_size_hhq_8bit_4N.txt
./AlphaTree config/config_hhq_8bit_4N_8000.txt | tee -a output/out_size_hhq_8bit_4N.txt
./AlphaTree config/config_hhq_8bit_4N_10000.txt | tee -a output/out_size_hhq_8bit_4N.txt

./AlphaTree config/config_hhq_8bit_8N_1000.txt | tee -a output/out_size_hhq_8bit_8N.txt
./AlphaTree config/config_hhq_8bit_8N_2000.txt | tee -a output/out_size_hhq_8bit_8N.txt
./AlphaTree config/config_hhq_8bit_8N_4000.txt | tee -a output/out_size_hhq_8bit_8N.txt
./AlphaTree config/config_hhq_8bit_8N_6000.txt | tee -a output/out_size_hhq_8bit_8N.txt
./AlphaTree config/config_hhq_8bit_8N_8000.txt | tee -a output/out_size_hhq_8bit_8N.txt
./AlphaTree config/config_hhq_8bit_8N_10000.txt | tee -a output/out_size_hhq_8bit_8N.txt

./AlphaTree config/config_hhq_16bit_4N_1000.txt | tee -a output/out_size_hhq_16bit_4N.txt
./AlphaTree config/config_hhq_16bit_4N_2000.txt | tee -a output/out_size_hhq_16bit_4N.txt
./AlphaTree config/config_hhq_16bit_4N_4000.txt | tee -a output/out_size_hhq_16bit_4N.txt
./AlphaTree config/config_hhq_16bit_4N_6000.txt | tee -a output/out_size_hhq_16bit_4N.txt
./AlphaTree config/config_hhq_16bit_4N_8000.txt | tee -a output/out_size_hhq_16bit_4N.txt
./AlphaTree config/config_hhq_16bit_4N_10000.txt | tee -a output/out_size_hhq_16bit_4N.txt

./AlphaTree config/config_hhq_16bit_8N_1000.txt | tee -a output/out_size_hhq_16bit_8N.txt
./AlphaTree config/config_hhq_16bit_8N_2000.txt | tee -a output/out_size_hhq_16bit_8N.txt
./AlphaTree config/config_hhq_16bit_8N_4000.txt | tee -a output/out_size_hhq_16bit_8N.txt
./AlphaTree config/config_hhq_16bit_8N_6000.txt | tee -a output/out_size_hhq_16bit_8N.txt
./AlphaTree config/config_hhq_16bit_8N_8000.txt | tee -a output/out_size_hhq_16bit_8N.txt
./AlphaTree config/config_hhq_16bit_8N_10000.txt | tee -a output/out_size_hhq_16bit_8N.txt

./AlphaTree config/config_hhq_32bit_4N_1000.txt | tee -a output/out_size_hhq_32bit_4N.txt
./AlphaTree config/config_hhq_32bit_4N_2000.txt | tee -a output/out_size_hhq_32bit_4N.txt
./AlphaTree config/config_hhq_32bit_4N_4000.txt | tee -a output/out_size_hhq_32bit_4N.txt
./AlphaTree config/config_hhq_32bit_4N_6000.txt | tee -a output/out_size_hhq_32bit_4N.txt
./AlphaTree config/config_hhq_32bit_4N_8000.txt | tee -a output/out_size_hhq_32bit_4N.txt
./AlphaTree config/config_hhq_32bit_4N_10000.txt | tee -a output/out_size_hhq_32bit_4N.txt

./AlphaTree config/config_hhq_32bit_8N_1000.txt | tee -a output/out_size_hhq_32bit_8N.txt
./AlphaTree config/config_hhq_32bit_8N_2000.txt | tee -a output/out_size_hhq_32bit_8N.txt
./AlphaTree config/config_hhq_32bit_8N_4000.txt | tee -a output/out_size_hhq_32bit_8N.txt
./AlphaTree config/config_hhq_32bit_8N_6000.txt | tee -a output/out_size_hhq_32bit_8N.txt
./AlphaTree config/config_hhq_32bit_8N_8000.txt | tee -a output/out_size_hhq_32bit_8N.txt
./AlphaTree config/config_hhq_32bit_8N_10000.txt | tee -a output/out_size_hhq_32bit_8N.txt

./AlphaTree config/config_hhq_64bit_4N_1000.txt | tee -a output/out_size_hhq_64bit_4N.txt
./AlphaTree config/config_hhq_64bit_4N_2000.txt | tee -a output/out_size_hhq_64bit_4N.txt
./AlphaTree config/config_hhq_64bit_4N_4000.txt | tee -a output/out_size_hhq_64bit_4N.txt
./AlphaTree config/config_hhq_64bit_4N_6000.txt | tee -a output/out_size_hhq_6  4bit_4N.txt
./AlphaTree config/config_hhq_64bit_4N_8000.txt | tee -a output/out_size_hhq_64bit_4N.txt
./AlphaTree config/config_hhq_64bit_4N_10000.txt | tee -a output/out_size_hhq_64bit_4N.txt

./AlphaTree config/config_hhq_64bit_8N_1000.txt | tee -a output/out_size_hhq_64bit_8N.txt
./AlphaTree config/config_hhq_64bit_8N_2000.txt | tee -a output/out_size_hhq_64bit_8N.txt
./AlphaTree config/config_hhq_64bit_8N_4000.txt | tee -a output/out_size_hhq_64bit_8N.txt
./AlphaTree config/config_hhq_64bit_8N_6000.txt | tee -a output/out_size_hhq_64bit_8N.txt
./AlphaTree config/config_hhq_64bit_8N_8000.txt | tee -a output/out_size_hhq_64bit_8N.txt
./AlphaTree config/config_hhq_64bit_8N_10000.txt | tee -a output/out_size_hhq_64bit_8N.txt


./AlphaTree config/config_hq_8bit_4N_12000.txt
./AlphaTree config/config_hq_8bit_8N_12000.txt
./AlphaTree config/config_hq_16bit_4N_12000.txt
./AlphaTree config/config_hq_16bit_8N_12000.txt
./AlphaTree config/config_hhq_8bit_4N_12000.txt
./AlphaTree config/config_hhq_8bit_8N_12000.txt
./AlphaTree config/config_hhq_16bit_4N_12000.txt
./AlphaTree config/config_hhq_16bit_8N_12000.txt
./AlphaTree config/config_hhq_32bit_4N_12000.txt
./AlphaTree config/config_hhq_32bit_8N_12000.txt
./AlphaTree config/config_hhq_64bit_4N_12000.txt
./AlphaTree config/config_hhq_64bit_8N_12000.txt
