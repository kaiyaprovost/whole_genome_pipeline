#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=9999999:00:00
#PBS -N zip_old
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

cd /vz-nas2-archive/SLiM_Pipeline-Archive

for i in */*.*; do
if [ "${i: -2}" != "gz" ]; then
echo $i
#gzip -f $i
fi
done
