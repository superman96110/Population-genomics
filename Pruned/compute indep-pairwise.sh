#确保bim文件的第二列，SNP名称不能为"."，应为chr:position，或者rs号
#(base) [supeng@jianglin gwas]$ head 664_sheep_filter_maf005_geno01_mind01.bim
#1       1:3786  0       3786    A       G
#1       1:3800  0       3800    G       A
#1       1:3810  0       3810    T       A
#1       1:3820  0       3820    C       T
#1       1:3832  0       3832    G       A


plink --bfile 664_sheep_filter_maf005_geno01_mind01 --indep-pairwise 50 10 0.2 --chr-set 26 --keep-allele-order --out 664_sheep_filter_maf005_geno01_mind01_pruned

plink --bfile 664_sheep_filter_maf005_geno01_mind01 --extract 664_sheep_filter_maf005_geno01_mind01_pruned.prune.in --make-bed --chr-set 26 --keep-allele-order --out 664_sheep_filter_maf005_geno01_mind01_pruned501002
