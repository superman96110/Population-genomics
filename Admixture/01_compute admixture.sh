#计算admixture，需要pruned的bfile文件
#!/bin/bash
for k in {2..14}; do
    /home/jianglin/software/admixture_linux-1.3.0/admixture --cv sheep_snpname_filter_gtex_WGS_RNA_pruned_remove.bed $k | tee log${k}.out
done
