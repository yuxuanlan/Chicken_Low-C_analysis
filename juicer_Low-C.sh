# This is a clean version of 
# /ei/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Analysis.todelete/Juicer-1.6/Low-C/juicer_Low-C.v2.sh

root=/ei/scratch/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C
analysis_dir=$root/Analysis/Juicer-1.6/Low-C
QC_dir=/ei/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Analysis/QC
cd $analysis_dir

fi
# samtools used here 
# source package c92263ec-95e5-43eb-a527-8f1496d56f1a

# ==== 0. make bwa-mem2 index to use in Juicer ====
if true; then
    source bwa-mem2-2.2.1_CBG
    cd $root/References/Gallus_gallus/GRCg6a/

    wget https://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
    wget https://ftp.ensembl.org/pub/release-106/gtf/gallus_gallus/Gallus_gallus.GRCg6a.106.gtf.gz

    cd $root/References/Gallus_gallus/GRCg6a/bwa-mem2
    ln -vs $root/References/Gallus_gallus/GRCg6a/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz ./
    sbatch -p ei-short -c 1 -J bwa2_index -o bwa2_index.%j.%N.log --mem 100G -t 3:00:00 \
        --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
        --wrap="cd $root/References/Gallus_gallus/GRCg6a/bwa-mem2 && source bwa-mem2-2.2.1_CBG && \
                bwa-mem2 index Gallus_gallus.GRCg6a.dna.toplevel.fa"
fi
# ==== 1. prep ====
if true; then
    # make chromesome size file
    if true; then
        cd $root/References/Gallus_gallus/GRCg6a/bwa-mem2
        sbatch -p ei-short -c 1 -J faidx -o faidx.%j.%N.log --mem 20G -t 0:20:00 \
            --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
            --wrap="cd $root/References/Gallus_gallus/GRCg6a/bwa-mem2 && source package c92263ec-95e5-43eb-a527-8f1496d56f1a && \
                samtools faidx Gallus_gallus.GRCg6a.dna.toplevel.fa"

        cut -f1,2 Gallus_gallus.GRCg6a.dna.toplevel.fa.fai > Gallus_gallus.chrom_size
    fi

    # make restriction site
    if true; then
        python $root/Software/Juicer-1.6/misc/generate_site_positions.py MboI Gallus_gallus \
            $root/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.GRCg6a.dna.toplevel.fa 
        cp $root/References/Gallus_gallus_MboI.txt $analysis_dir/restriction_sites/

    fi
fi

# ==== 2. run juicer on Low-C data: didn't run arrowhead to make TADs ====
if true; then
    if true; then
        # sampleid=HTA15
        # sampleid=CCA3
        # sampleid=CCP6
        # sampleid=HTP16a
        # sampleid=AN 
        # sampleid=Pos
        sampleid=All
        if [[ $sampleid == Pos ]]; then
            if [ ! -d $analysis_dir/$sampleid.run2 ]; then mkdir $analysis_dir/$sampleid.run2; fi
            cd $analysis_dir/$sampleid.run2/fastq/
            ln -svf $QC_dir/trim_galore-0.6.10/HTP16a_R1_val_1.fq.gz ./HTP16a_R1.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/HTP16a_R2_val_2.fq.gz ./HTP16a_R2.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/CCP6_R1_val_1.fq.gz ./CCP6_R1.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/CCP6_R2_val_2.fq.gz ./CCP6_R2.fastq.gz 
        fi
        if [[ $sampleid == AN ]]; then
            if [ ! -d $analysis_dir/$sampleid.run2 ]; then mkdir $analysis_dir/$sampleid.run2; fi
            cd $analysis_dir/$sampleid.run2/fastq/
            ln -svf $QC_dir/trim_galore-0.6.10/HTA15_R1_val_1.fq.gz ./HTA15_R1.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/HTA15_R2_val_2.fq.gz ./HTA15_R2.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/CCA3_R1_val_1.fq.gz ./CCA3_R1.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/CCA3_R2_val_2.fq.gz ./CCA3_R2.fastq.gz 
        fi
        if [[ $sampleid == All ]]; then
            if [ ! -d $analysis_dir/$sampleid.run2/fastq ]; then mkdir -p $analysis_dir/$sampleid.run2/fastq; fi
            cd $analysis_dir/$sampleid.run2/fastq/
            ln -svf $QC_dir/trim_galore-0.6.10/HTA15_R1_val_1.fq.gz ./HTA15_R1.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/HTA15_R2_val_2.fq.gz ./HTA15_R2.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/CCA3_R1_val_1.fq.gz ./CCA3_R1.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/CCA3_R2_val_2.fq.gz ./CCA3_R2.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/HTP16a_R1_val_1.fq.gz ./HTP16a_R1.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/HTP16a_R2_val_2.fq.gz ./HTP16a_R2.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/CCP6_R1_val_1.fq.gz ./CCP6_R1.fastq.gz 
            ln -svf $QC_dir/trim_galore-0.6.10/CCP6_R2_val_2.fq.gz ./CCP6_R2.fastq.gz
        fi
        
        cd $analysis_dir/$sampleid.run2 # && rm -rf aligned/ debug/ HIC_tmp/ splits/
        # sbatch will submit a string of jobs and finish in no time. 

        # -t: bwa threads
        # -D Juicer scripts parent dir; -d topDir; -y site file
        # for single sample: -Q 6:00:00 -L 12:00:00
        # for AN/Pos: -Q 1-00:00:00 -L 3-12:00:00

        # run all stages
        if true; then
            cd $analysis_dir/$sampleid.run2 && rm -rf aligned/ debug/ HIC_tmp/ splits/
            sbatch -p ei-short -c 1 -J juicer.$sampleid -o juicer.$sampleid.%j.%N.log --mem 10G -t 0:10:00 \
                --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
                --wrap="cd $analysis_dir/$sampleid.run2 && \
                    $root/Software/Juicer-1.6/scripts/juicer_chr19.sh -t 20 -q ei-compute -l ei-compute -Q 2-00:00:00 -L 15-12:00:00 \
                    -y $root/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus_MboI.txt \
                    -D $root/Software/Juicer-1.6 \
                    -d $analysis_dir/$sampleid.run2 \
                    -z $root/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.GRCg6a.dna.toplevel.fa \
                    -p $root/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.chrom_size"
            # sq |grep a1776419231
        fi
        
    fi

fi

# ===== 4. post analysis: contact pair visualisation - anterior vs posterior. NO FILTERING  ====
if true; then
    # 4.1 BAM coverage visualisation on UCSC browser
    #   NOT USED: THE merged_dedup.sort.bam includes all reads, dup reads were tags only.
    if false; then
        cd $analysis_dir/time_analysis
        # merge anterior bam files
        if true; then
            # sort bams
            sampleid=CCA3
            # sampleid=CCP6
            # sampleid=HTP16a
            if [ ! -f $analysis_dir/$sampleid.run2/aligned/merged_dedup.sort.bam ]; then
                cd $analysis_dir/time_analysis
                sbatch -p ei-compute -c 4 -J samtools-sort.$sampleid -o samtools-sort.$sampleid.%j.%N.log --mem 100G -t 3:00:00 \
                    --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
                    --wrap="cd $analysis_dir/time_analysis && source samtools-1.11_CBG && \
                    samtools sort -@4 -o $analysis_dir/$sampleid.run2/aligned/merged_dedup.sort.bam \
                        $analysis_dir/$sampleid.run2/aligned/merged_dedup.bam && \
                    samtools index $analysis_dir/$sampleid.run2/aligned/merged_dedup.sort.bam"
                # CCA3: 20526095 done
                # CCP6: 20526113 done
                # HTP16a: 20526119 done
            fi

            # how many reads in merged_dedup.sort.bam?
            sbatch -p ei-compute -c 4 -J samtools-view -o samtools-view.%j.%N.log --mem 50G -t 12:00:00 \
                --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk  \
                --wrap="cd /ei/scratch/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Analysis/Juicer-1.6/Low-C/HTA15.run1/aligned \
                    && source samtools-1.11_CBG && samtools view -@4 merged_dedup.sort.bam|cut -f1|sort|uniq|wc -l"
            # 20598380

            # merge bams:
            sbatch -p ei-compute -c 4 -J samtools-merge.AN -o samtools-merge.AN.%j.%N.log --mem 100G -t 3:00:00 \
                    --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
                    --wrap="cd $analysis_dir/time_analysis && source samtools-1.11_CBG && \
                            samtools merge -@4 -f AN.merged_dedup.bam \
                                $analysis_dir/HTA15.run2/aligned/merged_dedup.sort.bam \
                                $analysis_dir/CCA3.run2/aligned/merged_dedup.sort.bam "
            
            # index bam:
            timetype=Pos 
            cd $analysis_dir/time_analysis
            sbatch -p ei-compute -c 4 -J samtools-index.$timetype -o samtools-merge.$timetype.%j.%N.log --mem 100G -t 3:00:00 \
                    --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
                    --wrap="cd $analysis_dir/time_analysis && source samtools-1.11_CBG && \
                            samtools index $timetype.merged_dedup.bam "

            # AN: 20529845 done
            timetype=AN 
            cd $analysis_dir/time_analysis

            # bedcov.bed file:
            sbatch -p ei-compute -c 1 -J genomecov_$timetype -o genomecov_$timetype.%j.%N.log --mem 50G -t 5:00:00 \
                --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
                --wrap="source bedtools-2.31.1_CBG && cd $analysis_dir/time_analysis && \
                    bedtools genomecov -d  -ibam $analysis_dir/time_analysis/AN.merged_dedup.bam \
                    -g $root/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.chrom_size > \
                    $analysis_dir/time_analysis/$timetype.merged_dedup.bedcov.bed" 
            head $timetype.merged_dedup.bedcov.bed
            # MT  1   24826
            # MT  2   24867

            # instead run these steps in sbatch
            if false; then
                # convert to proper bed format
                cat $timetype.merged_dedup.bedcov.bed | \
                    awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > $analysis_dir/time_analysis/$timetype.merged_dedup.bedcov.bedgraph

                # interset with enhancer.bed
                source bedtools-2.31.1_CBG
                bedtools intersect -wa -a $analysis_dir/time_analysis/$timetype.merged_dedup.bedcov.bedgraph \
                    -b $analysis_dir/enhancer.1.flank1k.2.bed > $analysis_dir/time_analysis/$timetype.merged_dedup.intersect.enhancer.bedgraph
            fi 
            timetype=Pos
            sbatch -p ei-compute -c 1 -J bedtools-intersect.$timetype -o bedtools-intersect.$timetype.%j.%N.log --mem 30G -t 03:00:00 \
                    --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
                    --wrap="cd $analysis_dir/time_analysis && source bedtools-2.31.1_CBG && \
                        cat $timetype.merged_dedup.bedcov.bed | \
                        awk '{print \$1\"\t\"\$2-1\"\t\"\$2\"\t\"\$3}' > $analysis_dir/time_analysis/$timetype.merged_dedup.bedcov.bedgraph && \
                        bedtools intersect -wa -a $analysis_dir/time_analysis/$timetype.merged_dedup.bedcov.bedgraph \
                            -b $analysis_dir/enhancer.1.flank1k.2.bed > \
                            $analysis_dir/time_analysis/$timetype.merged_dedup.intersect.enhancer.bedgraph"
            # 20555266

            /Users/lany/mount_drives/HPC_Projects_Scratch/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Software/wigToBigWig_macos/wigToBigWig \
                $timetype.merged_dedup.intersect.enhancer.bedgraph \
                /Users/lany/mount_drives/HPC_Projects_Scratch/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.chrom_size \
                $timetype.merged_dedup.intersect.enhancer.bw

            # double check to confirm the coverage
            if true; then
                # view coverage of an area in AN vs PoS
                less -N AN.merged_dedup.bedcov.bedgraph
                # 524028413 3       87378167        87378168        0
                # 524028414 3       87378168        87378169        0
                # 524028415 3       87378169        87378170        7
                # 524028416 3       87378170        87378171        7
                
                source samtools-1.11_CBG
                samtools index AN.merged_dedup.bam
                samtools tview  -p 3:87378167 AN.merged_dedup.bam

                less -N AN.merged_dedup.bedcov.bedgraph
                # 524028413 3       87378167        87378168        28
                # 524028414 3       87378168        87378169        28
                # 524028415 3       87378169        87378170        28
                # 524028416 3       87378170        87378171        28
            fi
        fi

        # bedgraph to bw
        if true; then
            /Users/lany/mount_drives/HPC_Projects_Scratch/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Software/wigToBigWig_macos/wigToBigWig \
                merged_dedup.bedcov.bedgraph \
                /Users/lany/mount_drives/HPC_Projects_Scratch/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.chrom_size \
                merged_dedup.bedcov.bw

        fi
    fi

    # 4.2 contact pair coverage visualisation using UCSC browser
    if true; then
        cd $analysis_dir/time_analysis.merged30

        timetype=AN 
        if [[ $timetype == Pos ]]; then
            cat $analysis_dir/HTP16a.run2/aligned/merged30.txt $analysis_dir/CCP6.run2/aligned/merged30.txt | \
                cut -d' ' -f2,3,6,7 > $timetype.merged30.txt
        fi
        if [[ $timetype == AN ]]; then
            cat $analysis_dir/HTA15.run2/aligned/merged30.txt $analysis_dir/CCA3.run2/aligned/merged30.txt | \
                cut -d' ' -f2,3,6,7 > $timetype.merged30.txt
        fi
        # head $timetype.merged30.txt
        # pair 1;   pair 2;
        # 10 562716 10 3743033
        # 10 743546 10 3665128

        # this include some error entries, to be removed (not many):
        # 0 MT -129 0 16 MT 1 1
        # 0 MT -2 0 16 MT 177 1

        # separate pair 1 and pair2, make into bed format
        awk -F" " '{if($2>=0 && $4>=0) {print $1"\t"$2-1"\t"$2}}' $timetype.merged30.txt > tmp1
        awk -F" " '{if($2>=0 && $4>=0) {print $3"\t"$4-1"\t"$4}}' $timetype.merged30.txt > tmp2

        cat tmp1 tmp2 |grep -E -f grep.chr.txt | sort -k 1,1 -k2,2n|uniq | awk '$2>=0'> $timetype.merged30.bed
        # rm tmp1 tmp2 $timetype.merged30.txt
        head $timetype.merged30.bed
        # 1   216 217
        # 1   216 217
        # cat $timetype.merged30.bed | sed -i 's/ /\t/g' > tmp.bed
        source bedtools-2.31.1_CBG
        bedtools genomecov -i $timetype.merged30.bed \
            -g $root/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.chrom_size -bga \
            > $timetype.merged30.bed_genomecov.bedgraph

        # convert to bw - THIS IS THE BW FILE FINALLY UPLOADED TO BROWSER!
        /Users/lany/mount_drives/HPC_Projects_Scratch/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Software/wigToBigWig_macos/wigToBigWig \
                $timetype.merged30.bed_genomecov.bedgraph \
                /Users/lany/mount_drives/HPC_Projects_Scratch/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.chrom_size \
                $timetype.merged30.bw
    fi
fi

# ===== 7. gene associated with UNFILTERED contact pairs (cp) - UNFILTERED =====
#  contact pairs: inter contig, +-5kb enhancer region only 
if true; then 
    # 7.1 prepare flanking region, use enhancer.sort.2.bed, where enhancer name of TSHZ2 is changed to BMP7
    if true; then
        cd $analysis_dir/time_analysis.merged30/nofilter

        # flanking over (include) enhancer regions
        bedtools flank -b 5000 -i $analysis_dir/enhancer.sort.2.bed -g $root/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.chrom_size \
            > ./enhancer.1.flank5k.bed 
        cat $analysis_dir/enhancer.sort.2.bed ./enhancer.1.flank5k.bed |sort -k1,1 -k2,2n | bedtools merge > ./enhancer.2.flank5k.bed 
        # added gene name to bed: 
        head enhancer.3.flank5k.bed
        # 10  1568138 1578639 RBPMS2
        # 15  12698285    12709579    TBX5
    fi

    # 7.2 modify cp, remove: 
    # 1. inter contig pairs; 
    if true; then
        source bedtools-2.31.1_CBG

        cd $analysis_dir/time_analysis.merged30/nofilter/map_to_gene
        timetype=AN
        # timetype=Pos
        # head ../$timetype.merged30.txt
        # 10 562716 10 3743033
        # 10 743546 10 3665128

        echo "intra_cp,in_enhancer,with_gene_names" > $timetype.process_stat.txt 
        # remove inter contig cp and short distant cp
        awk '$1==$3 && $2>=0 && $4>=0' ../../$timetype.merged30.txt > $timetype.intra_cp.txt
        # 10 562716 10 3743033
        # 10 743546 10 3665128
        n1=$(cat $timetype.intra_cp.txt|wc -l)
        #  AN: 296,053

        # split intra_cp.txt into mate1 and mate2, move the mate pair to last 2 columns
        awk -F' ' '$2-1>0 {print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' $timetype.intra_cp.txt > $timetype.intra_cp.M1.bed
        awk -F' ' '$4-1>0 {print $3"\t"$4-1"\t"$4"\t"$1"\t"$2}' $timetype.intra_cp.txt > $timetype.intra_cp.M2.bed
        
        # working on M1, search which gene its mate pair touches
        if true; then
            # mate 1 needs to be inside enhancer region.
            bedtools intersect -wa -wb -a $timetype.intra_cp.M1.bed -b ../enhancer.3.flank5k.bed | \
                cut -f1-5,9 | sort -k1,1 -k2,2n > $timetype.intra_cp.M1.gene.bed
            n2=$(cat $timetype.intra_cp.M1.gene.bed|wc -l)
            # AN: 47
            # head AN.intra_cp.M1.gene.bed
            # M1                  M1_mate         M1_gene
            # 10  1570747 1570748 10  14263230    RBPMS2
            # 10  1572571 1572572 10  1644181 RBPMS2

            # find genes where its mate touches
            # make its mate bedfile:

            # awk '{print $4"\t"$5-1"\t"$5"\t"$1"\t"$2"\t"$3"\t"$6}' $timetype.intra_cp.mate1.gene.bed|head
            # mate                        M1                      M1_gene
            # 10  14263229    14263230    10  1570747 1570748 10  RBPMS2
            awk '{print $4"\t"$5-1"\t"$5"\t"$1"\t"$2"\t"$3"\t"$6}' $timetype.intra_cp.M1.gene.bed |\
                bedtools intersect -wa -wb -b ../../Gallus_gallus.GRCg6a.106.gff3 -a - |\
                cut -f1-8,16|sed 's/ID=.*;geneID=\(.*\)/geneID=\1/' |\
                grep "geneID="|sort -k1,1 -k2,2n > $timetype.tmp.M1.bed 
            # head $timetype.tmp.M1.bed 
            # mate                M1                  M1_gene M2_gene
            # 10  1644180 1644181 10  1572571 1572572 RBPMS2  10  geneID=ENSGALG00000030672;gene_name=TRIP4
            # 10  1644180 1644181 10  1572571 1572572 RBPMS2  10  geneID=ENSGALG00000040087;gene_name=ZNF609

            # tidy up
            awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3"\t"$9}' $timetype.tmp.M1.bed | sort | uniq > $timetype.intra_cp.M1_mate.bed 
            # note the repeat below as two genes in same coord
            # 10  1572571 1572572 RBPMS2  10  1644180 1644181 geneID=ENSGALG00000030672;gene_name=TRIP4
            # 10  1572571 1572572 RBPMS2  10  1644180 1644181 geneID=ENSGALG00000040087;gene_name=ZNF609
            n3=$(cut -f1-4 $timetype.intra_cp.M1_mate.bed|sort|uniq|wc -l)
            echo "M1,$n1,$n2,$n3" >> $timetype.process_stat.txt

        fi

        # working on M2, search which gene its mate pair touches
        if true; then
            # mate 1 needs to be inside enhancer region.
            bedtools intersect -wa -wb -a $timetype.intra_cp.M2.bed -b ../enhancer.3.flank5k.bed | \
                cut -f1-5,9 | sort -k1,1 -k2,2n > $timetype.intra_cp.M2.gene.bed
            n2=$(cat $timetype.intra_cp.M2.gene.bed|wc -l)
            
            awk '{print $4"\t"$5-1"\t"$5"\t"$1"\t"$2"\t"$3"\t"$6}' $timetype.intra_cp.M2.gene.bed |\
                bedtools intersect -wa -wb -b ../../Gallus_gallus.GRCg6a.106.gff3 -a - |\
                cut -f1-8,16|sed 's/ID=.*;geneID=\(.*\)/geneID=\1/' |\
                grep "geneID="|sort -k1,1 -k2,2n > $timetype.tmp.M2.bed 

            # tidy up
            awk '{print $4"\t"$5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3"\t"$9}' $timetype.tmp.M2.bed | sort | uniq > $timetype.intra_cp.M2_mate.bed 
            # note the repeat below as two genes in same coord
            # 10  1572571 1572572 RBPMS2  10  1644180 1644181 geneID=ENSGALG00000030672;gene_name=TRIP4
            # 10  1572571 1572572 RBPMS2  10  1644180 1644181 geneID=ENSGALG00000040087;gene_name=ZNF609
            n3=$(cut -f1-4 $timetype.intra_cp.M2_mate.bed|sort|uniq|wc -l)
            echo "M2,$n1,$n2,$n3" >> $timetype.process_stat.txt
        fi

        # merge M1 and M2 outputs
        cat $timetype.intra_cp.M1_mate.bed $timetype.intra_cp.M2_mate.bed | sort -k1,1 -k2,2n > $timetype.intra_cp.M1M2gene.bed

    fi
fi

