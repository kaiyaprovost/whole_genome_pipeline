!/bin/bash

# STEP04-GATK_SNP_GVCF_combine_bel111_BOTH.job
# declare -a arr=("WA07" "WB07" "WB11" "WD02" "WE08" "WF02" "WF05" "WG05" "WH08")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/bellii/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_belBOTH.job;
# done; echo ""; done;
# 
# STEP04-GATK_SNP_GVCF_combine_bel111_BOTH.job
# declare -a arr=("WA06" "WC01" "WD07" "WD11" "WF03")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/bellii/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_belBOTH.job;
# done; echo ""; done;

# STEP04-GATK_SNP_GVCF_combine_cri111_BOTH.job
# declare -a arr=("245111_P002_WE10" "245111_P002_WA09" "245111_P002_WH05")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/crissale/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_criBOTH.job;
# done; echo ""; done;

# STEP04-GATK_SNP_GVCF_combine_cri111_BOTH.job
# declare -a arr=("245111_P002_WB06" "245111_P002_WC06")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/crissale/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_criBOTH.job;
# done; echo ""; done;
# 
# STEP04-GATK_SNP_GVCF_combine_fla111_BOTH.job
# declare -a arr=("245111_P002_WA08" "245111_P002_WB05" "245111_P002_WC11" "245111_P002_WD06" "245111_P002_WG01" "245111_P002_WH01") 
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/flaviceps/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_flaBOTH.job;
# done; echo ""; done;
# 
# STEP04-GATK_SNP_GVCF_combine_fla111_BOTH.job
# declare -a arr=("245111_P002_WA02" "245111_P002_WB02" "245111_P002_WE02" "245111_P002_WE03" "245111_P002_WB01" "245111_P002_WD01" "245111_P002_WF01")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/flaviceps/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_flaBOTH.job;
# done; echo ""; done;

# STEP04-GATK_SNP_GVCF_combine_nit111_BOTH.job
# declare -a arr=("245111_P002_WC05" "245111_P002_WD05" "245111_P002_WE07" "245111_P002_WA04" "245111_P002_WH03")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/nitens/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_nitBOTH.job;
# done; echo ""; done;

# STEP04-GATK_SNP_GVCF_combine_bil107-109-111_BOTH.job
# declare -a arr=("245107_P01_WD01" "245107_P01_WE01")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/bilineata/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bilBOTH.job;
# done; echo ""; done;

# declare -a arr=("245109_P002_WA07" "245109_P002_WB05" "245109_P002_WE01" "245109_P002_WG01")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/bilineata/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bilBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WA03" "245111_P002_WE05" "245111_P002_WG09")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/bilineata/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bilBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_bil107-109-111_BOTH.job
# declare -a arr=("245107_P01_WA01" "245107_P01_WB01" "245107_P01_WC01")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/bilineata/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bilBOTH.job;
# done; echo ""; done;

# declare -a arr=("245109_P002_WD03" "245109_P002_WD05" "245109_P002_WE02" "245109_P002_WF01" "245109_P002_WF02" "245109_P002_WG06")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/bilineata/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bilBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WA01")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/bilineata/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bilBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_cur109-111_BOTH.job
# declare -a arr=("245109_P002_WB04" "245109_P002_WF06" "245109_P002_WG02")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/curvirostre/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_curBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WC02" "245111_P002_WE06" "245111_P002_WF10" "245111_P002_WG10" "245111_P002_WH10")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/curvirostre/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_curBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_cur109-111_BOTH.job
# declare -a arr=("245109_P002_WC06" "245109_P002_WD02" "245109_P002_WH06")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/curvirostre/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_curBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WA12" "245111_P002_WD10" "245111_P002_WG03")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/curvirostre/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_curBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_mel109-111_BOTH.job
# declare -a arr=("245109_P002_WA04" "245109_P002_WB03" "245109_P002_WC05" "245109_P002_WE04")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/melanura/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_melBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WB10" "245111_P002_WC10" "245111_P002_WE01" "245111_P002_WH09")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/melanura/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_melBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_mel109-111_BOTH.job
# declare -a arr=("245109_P002_WA01" "245109_P002_WA02" "245109_P002_WB02" "245109_P002_WD04" "245109_P002_WG03" "245109_P002_WH03")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/melanura/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_melBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WB03" "245111_P002_WC03")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/melanura/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_melBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_fus107-109-111_BOTH.job
# declare -a arr=("245107_P01_WG02")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/fusca/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_fusBOTH.job;
# done; echo ""; done;

# declare -a arr=("245109_P002_WB01" "245109_P002_WC01" "245109_P002_WG04" "245109_P002_WH04")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/fusca/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_fusBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WA05" "245111_P002_WF06" "245111_P002_WG06")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/fusca/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_fusBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_fus107-109-111_BOTH.job
# declare -a arr=("245107_P01_WC02" "245107_P01_WD02" "245107_P01_WE02" "245107_P01_WF02")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/fusca/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_fusBOTH.job;
# done; echo ""; done;

# declare -a arr=("245109_P002_WA03" "245109_P002_WB06" "245109_P002_WC04" "245109_P002_WH02")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/fusca/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_fusBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WA11" "245111_P002_WF08")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/fusca/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_fusBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_bru107-109-111_BOTH.job
# declare -a arr=("245107_P01_WB02")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/brunneicapillus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bruBOTH.job;
# done; echo ""; done;

# declare -a arr=("245109_P002_WE06")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/brunneicapillus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bruBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WC08" "245111_P002_WC09" "245111_P002_WD03" "245111_P002_WD09" "245111_P002_WF09" "245111_P002_WG02" "245111_P002_WH02")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/brunneicapillus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bruBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_bru107-109-111_BOTH.job
# declare -a arr=("245107_P01_WA02" "245107_P01_WF01" "245107_P01_WG01" "245107_P01_WH01")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/brunneicapillus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bruBOTH.job;
# done; echo ""; done;

# declare -a arr=("245109_P002_WE03" "245109_P002_WF04")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/brunneicapillus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bruBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WB08" "245111_P002_WD08" "245111_P002_WG08")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/brunneicapillus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_bruBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_sin107-109-111_BOTH.job
# declare -a arr=("245107_P01_WD03")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/sinuatus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_sinBOTH.job;
# done; echo ""; done;

# declare -a arr=("245107_P01_WC03")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/sinuatus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_sinBOTH.job;
# done; echo ""; done;

# declare -a arr=("245109_P002_WA05" "245109_P002_WA06" "245109_P002_WC02" "245109_P002_WC03" "245109_P002_WD01" "245109_P002_WD06" "245109_P002_WH01")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/sinuatus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_sinBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WA10" "245111_P002_WB09")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/sinuatus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_sinBOTH.job;
# done; echo ""; done;

#STEP04-GATK_SNP_GVCF_combine_sin107-109-111_BOTH.job
# declare -a arr=("245107_P01_WA03" "245107_P01_WB03" "245107_P01_WH02")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245107/sinuatus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_sinBOTH.job;
# done; echo ""; done;

# declare -a arr=("245109_P002_WF03")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245109/Plate2/sinuatus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_sinBOTH.job;
# done; echo ""; done;

# declare -a arr=("245111_P002_WC07" "245111_P002_WE09" "245111_P002_WH06")
# cd /vz-nas1-active/ProcessedGenomicReads/AMN_245111/Plate2/sinuatus/December_Rerun/BAM_files/not_merged/
# for i in "${arr[@]}"; do echo "XXXXXXXXXX"; echo "$i"; for filename in ${PWD}/*${i}*dedup.reheadered.bam; do echo $filename;
# prefix=${filename#$PWD/}; 
# prefix=$(echo $prefix | cut -f 1 -d '.')
# qsub -v bam=$filename,name=$prefix /home/kprovost/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP06-GATK_SNP_BQSR_sinBOTH.job;
# done; echo "XXXXXXXXXX"; done;

