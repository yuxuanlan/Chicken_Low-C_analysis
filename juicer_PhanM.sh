# this is a tidy version of 
# /ei/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Analysis/Juicer-1.6/PhanM/juicer_PhanM.v2.sh
# version history: 
# v1: run Arrowhead with -r 5000 for high resolution HIC data
# v2: -r 10000 for medium resolution

root=/ei/scratch/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C
analysis_dir=$root/Analysis/Juicer-1.6/PhanM
cd $analysis_dir
software_script_dir=/ei/projects/5/5c5d0df3-b131-4ef6-8c6d-e7c45e3358ce/results/Audrea_Munsterberg/Low-C/Software/Juicer-1.6/scripts
# sample=HH22
# download hic data
if true; then
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263755/suppl/GSE263755%5FHiC%5Fheart%5FHH22%5FgalGal6%5FWT%5Fmerged%2Ehic
fi

# ==== 1. call TAD, -r 5k (high resolution) output to TADs.run2/ ====
if true; then
    
    cd $analysis_dir/TADs.run2
    sbatch -p ei-short -c 1 -J arrowhead-PhanM -o arrowhead-PhanM.%j.%N.log --mem 40G -t 2:00:00 \
            --mail-type=ALL --mail-user=yuxuan.lan@earlham.ac.uk \
            --wrap="cd $analysis_dir/TADs.run2 && source jre-8u92 && \
                $software_script_dir/juicer_tools arrowhead -r 10000 \
                $analysis_dir/TADs.run2/GSE263755_HiC_heart_HH22_galGal6_WT_merged.hic \
                $analysis_dir/TADs.run2/"
    # 24415465
fi
# ==== 2. intersect TADs with enhancer regions, try different flanking sizes ====

# ==== 2.3 enhancer+-2Mb regions overlapping gene (ensembl-106 annotation) ====
if true; then
    cd $analysis_dir/TADs.run2/domain_analysis.ensembl106
    source bedtools-2.31.1_CBG

    # make 2M flanking on enhancers: done in 2.2
    if false; then
        bedtools flank -b 2000000 -i Putative_enhancers_galGal6.bed -g $root/References/Gallus_gallus/GRCg6a/bwa-mem2/Gallus_gallus.chrom_size \
            > ./enhancer.1.flank2M.bed 
        cat Putative_enhancers_galGal6.bed ./enhancer.1.flank2M.bed |sort -k1,1n -k2,2n | bedtools merge > ./enhancer.2.flank2M.bed 
        # added gene name to bed manually: 
        cat enhancer.3.flank2M.bed
        # 3   85379126    89379959    BMP5_3-870  0   +
        # 10  0   3573639 RBPMS2_10-588   0   +
        # 15  10703285    13062184    TBX5_15-768-1407    0   +
        # 18  0   3825973 MYOCD1_18-768   0   +
        # 20  9699877 13700378        BMP7_20-660 0   +
        # 26  2966342 6055710 APOBEC2_26-743  0   +
    fi
    
    # overlap TADs with 7 genes
    if true; then
        # convert bedpe to BED format
        if [ ! -f ../10000_blocks.bed ]; then 
            cut -f1-3,8 ../10000_blocks.bedpe | sort -k1,1 -k2,2n > ../10000_blocks.bed
        fi
        bedtools intersect -wa -wb -b ./gene_regions_ensembl106.bed \
            -a ../10000_blocks.bed | sort -k1,1n -k2,2n > TADs_in_gene.bed
        cat TADs_in_gene.bed
        # 3   87340000      87540000    1.4302139926634538  3   87360676    87415579    BMP5_gene
        # 10  1460000       1760000 0.9146852753055239  10  1554284 1575874 RBPMS2_gene
        # 15  12530000      12730000    1.2259513917204536  15  12598201    12638563    TBX5_gene
        # 15  12530000      12950000    1.420341346461756   15  12598201    12638563    TBX5_gene
        # 18  700000        1040000 1.418649198186509   18  808100  1027820 MYOCD_gene
        # 20  11640000      11970000    1.0715921250402447  20  11682314    11723878    BMP7_gene

        # from 5k TADs:
        # cat TADs_in_gene.bed
        # 15      12535000        12735000        1.441539498976025       15      12598201        12638563        TBX5_gene
        # 15      12540000        12610000        0.9372958534364728      15      12598201        12638563        TBX5_gene
        # 15      12620000        12720000        1.03606711942831        15      12598201        12638563        TBX5_gene
        # 20      11650000        11745000        1.2035528767196726      20      11682314        11723878        BMP7_gene
    fi

    # overlap TADs to enhancer
    if true; then
        # convert bedpe to BED format
        if [ ! -f ../10000_blocks.bed ]; then 
            cut -f1-3,8 ../10000_blocks.bedpe | sort -k1,1 -k2,2n > ../10000_blocks.bed
            head ../10000_blocks.bedpe | less -S
            # x1=y1, x2=y2, chr1=chr2
            # chr1  x1          x2          chr2    y1          y2          name    score               strand1 strand2 color       score               uVarScore           lVarScore           upSign              loSign
            # 10    10595000    11135000    10      10595000    11135000    .       0.42272225863106483 .       .       255,255,0   0.42272225863106483 0.3539977823229645  0.34434262744442734 0.4107744107744108  0.5016835016835017

        fi        
        
        # Find TADs overlap with enhancer
        bedtools intersect -wa -wb -b enhancer.3.flank2M.bed \
            -a ../10000_blocks.bedpe | cut -f1-3,8,17-20,22 | sort -k1,1n -k2,2n > TADs_in_enhancer.flank2M.bed 
        
        # one enhancer can overlap multiple TADs
        head TADs_in_enhancer.flank2M.bed
        # TADs                                            enhancer
        # 3   85550000    85730000    1.3036238402237244  3   85379126    89379959    BMP5_3-870  +
        # 3   85740000    86070000    0.847350645110547   3   85379126    89379959    BMP5_3-870  +
        # 3   86710000    87010000    1.2100747168279211  3   85379126    89379959    BMP5_3-870  +
    fi

    # TADs contain 7 genes and its enhancer: 
    if true; then
        echo "# col 1-9: TAD with enhancer; col 10-17: TAD with gene" >  TADs_in_enhancer_gene.bed
        bedtools intersect -wa -wb -a TADs_in_enhancer.flank2M.bed -b TADs_in_gene.bed >> TADs_in_enhancer_gene.bed

        cat TADs_in_enhancer_gene.bed
        # col 1-9: TAD with enhancer; col 10-17: TAD with gene
        # 3   87340000      87540000    1.4302139926634538  3   85379126    89379959    BMP5_3-870          +   3   87340000    87540000    1.4302139926634538  3   87360676    87415579    BMP5_gene
        # 10  1460000       1760000     0.9146852753055239  10  0           3573639     RBPMS2_10-588       +   10  1460000     1760000     0.9146852753055239  10  1554284     1575874     RBPMS2_gene
        # 10  1610000       1760000     0.8720830808111032  10  0           3573639     RBPMS2_10-588       +   10  1460000     1760000     0.9146852753055239  10  1554284     1575874     RBPMS2_gene
        # 15  12530000      12730000    1.2259513917204536  15  10703285    13062184    TBX5_15-768-1407    +   15  12530000    12730000    1.2259513917204536  15  12598201    12638563    TBX5_gene
        # 15  12530000      12730000    1.2259513917204536  15  10703285    13062184    TBX5_15-768-1407    +   15  12530000    12950000    1.420341346461756   15  12598201    12638563    TBX5_gene
        # 15  12530000      12950000    1.420341346461756   15  10703285    13062184    TBX5_15-768-1407    +   15  12530000    12730000    1.2259513917204536  15  12598201    12638563    TBX5_gene
        # 15  12530000      12950000    1.420341346461756   15  10703285    13062184    TBX5_15-768-1407    +   15  12530000    12950000    1.420341346461756   15  12598201    12638563    TBX5_gene
        # 18  700000        1040000     1.418649198186509   18  0           3825973     MYOCD1_18-768       +   18  700000      1040000     1.418649198186509   18  808100      1027820     MYOCD_gene
        # 20  11640000      11970000    1.0715921250402447  20  9699877     13700378    BMP7_20-660         +   20  11640000    11970000    1.0715921250402447  20  11682314    11723878    BMP7_gene
        # 20  11770000      11960000    1.200501330986765   20  9699877     13700378    BMP7_20-660         +   20  11640000    11970000    1.0715921250402447  20  11682314    11723878    BMP7_gene
        # TADs in 5k res
        # cat TADs_in_enhancer_gene.bed
        # 15  12535000    12735000    1.441539498976025   15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 15  12535000    12735000    1.441539498976025   15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 15  12535000    12735000    1.441539498976025   15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 15  12540000    12610000    0.9372958534364728  15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 15  12540000    12610000    0.9372958534364728  15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 15  12540000    12610000    0.9372958534364728  15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 15  12620000    12720000    1.03606711942831    15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 15  12620000    12720000    1.03606711942831    15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 15  12620000    12720000    1.03606711942831    15  10703285    13062184    TBX5_15-768-1407    +   15  12598201    12638563    TBX5_gene
        # 20  11650000    11745000    1.2035528767196726  20  9699877 13700378    BMP7_20-660 +   20  11682314    11723878    BMP7_gene

    fi
fi

