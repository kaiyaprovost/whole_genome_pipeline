#!/bin/bash

qsub -v bamlist=/home/kprovost/nas3/ANGSD_pipeline/Parrots.bamlist /home/kprovost/nas3/ANGSD_pipeline/A6.ANGSD_GL_Parrots_DXY.job 

qsub -v bamlist=/home/kprovost/nas3/ANGSD_pipeline/Parrots.bamlist /home/kprovost/nas3/ANGSD_pipeline/A6.ANGSD_GL_Parrots_DXY_parus.job 

qsub -v bamlist=/home/kprovost/nas3/ANGSD_pipeline/Parrots.bamlist /home/kprovost/nas3/ANGSD_pipeline/A6.ANGSD_GL_Parrots.job 

qsub -v bamlist=/home/kprovost/nas3/ANGSD_pipeline/Parrots.bamlist /home/kprovost/nas3/ANGSD_pipeline/A6.ANGSD_GL_Parrots_parus.job 
