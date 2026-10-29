# Author: Yuxuan Lan
# Date: 2026-10-28
# This is a tidy version of /ei/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Analysis/QC/QC_Low-C.v1.sh

# ==== 0. RUN BEFORE ALL CODES =====
# export PATH=/ei/cb/common/Scripts:$PATH 
sampleid=HTA15
# sampleid=HTP16a
# sampleid=CCA3
# sampleid=CCP6

root=/ei/scratch/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C
analysis_dir=$root/Analysis/QC
cd $analysis_dir
# samtools used here 
source package c92263ec-95e5-43eb-a527-8f1496d56f1a

# real reads resite under 
rawread_root=/ei/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/data/raw/Audrea_Munsterberg/Low-C/$sampleid
# ls -l $rawread_root| egrep '^d'|wc -l
symlink_dir=$root/Reads

# ==== 1. download chicken reference
if true; then
    cd $root/References/
    wget https://ftp.ensembl.org/pub/release-115/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz \
        
    wget https://ftp.ensembl.org/pub/release-115/gff3/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gff3.gz

    # make bwa index
    sbatch -p ei-short -c 1 -J bwa_index -o bwa_index.%j.%N.log --mem 10G \
        --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
        --wrap="cd $root/References/ && source bwa-0.7.17 && \
                bwa index Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"

    # make for bwa-mem2
    source bwa-mem2-2.2.1_CBG
    cd $root/References/bwa-mem2
    ln -vs $root/References/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa ./
    sbatch -p ei-short -c 1 -J bwa2_index -o bwa2_index.%j.%N.log --mem 100G \
        --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
        --wrap="cd $root/References/bwa-mem2/ && source bwa-mem2-2.2.1_CBG&& \
                bwa-mem2 index Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
    # 18688815
    
    # index fa
    cd $root/References/bwa-mem2
    sbatch -p ei-short -c 1 -J bwa2_index -o bwa2_index.%j.%N.log --mem 20G -t 0:20:00 \
        --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
        --wrap="cd $root/References/bwa-mem2/ && source package c92263ec-95e5-43eb-a527-8f1496d56f1a && \
            samtools faidx $root/References/bwa-mem2/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
fi
# ==== 2. symlink raw reads ====
if true; then
    gzfile=$(ls $rawread_root/*_R1.fastq.gz)
    ln -vs $gzfile $symlink_dir/${sampleid}_R1.fastq.gz
    gzfile=$(ls $rawread_root/*_R2.fastq.gz)
    ln -vs $gzfile $symlink_dir/${sampleid}_R2.fastq.gz
fi

# ==== 3. trim-galore ====
if true; then
    cd $analysis_dir/trim_galore-0.6.10
    
    sbatch -p ei-short -c 4 -J trimgalore_$sampleid -o trimgalore_$sampleid.%j.%N.log --mem 20G \
        --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
        --wrap="source package 04b61fb6-8090-486d-bc13-1529cd1fb791 \
                && cd $analysis_dir/trim_galore-0.6.10/ && \
                trim_galore --cores 4 --fastqc --output_dir ${analysis_dir}/trim_galore-0.6.10/ \
                --paired ${symlink_dir}/${sampleid}_R1.fastq.gz ${symlink_dir}/${sampleid}_R2.fastq.gz"
    # 18607252
    #     JobID        AllocNodes        NodeList                        JobName      User    Account    Elapsed  Timelimit     MaxRSS     ReqMem  MaxVMSize      NCPUS ExitCode
    # ------------ ---------- --------------- ------------------------------ --------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- --------
    # 18607252              1         e512n21               trimgalore_HTA15      lany         ei   05:25:22   12:00:00                  300G                     4      0:0
    # 18607252.ba+          1         e512n21                          batch                   ei   05:25:22               525628K            158719972K          4      0:0
    # 18607252.ex+          1         e512n21                         extern                   ei   05:25:23                   44K                   44K          4      0:0

    # stats on trim-galore: NO. THESE CAN BE RETRIEVED FROM fastqc report AND trimming_report.txt
    if false; then
        sbatch -p ei-short -c 1 -J trimgalore_stat_$sampleid -o trimgalore_stat_$sampleid.%j.%N.log --mem 20G \
        --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
        --wrap="cd $analysis_dir/trim_galore-0.6.10 && \
            zcat $symlink_dir/${sampleid}_R1.fastq.gz | wc -l && \
            zcat ${analysis_dir}/trim_galore-0.6.10/${sampleid}_R1_val_1.fq.gz | wc -l "

        sbatch -p ei-short -c 1 -J trimgalore_stat_$sampleid -o trimgalore_stat_$sampleid.%j.%N.log --mem 20G \
        --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
        --wrap="cd $analysis_dir/trim_galore-0.6.10 && \
            zcat $symlink_dir/${sampleid}_R2.fastq.gz | wc -l && \
            zcat ${analysis_dir}/trim_galore-0.6.10/${sampleid}_R2_val_2.fq.gz | wc -l"

    fi
fi

#