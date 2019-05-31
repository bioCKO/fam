# 基因家族全基因组分析

以梨树NAC基因家族为例，演示基因家族全基因分析内容和方法。

## 系统要求
本文的代码均在Ubuntu18.04操作系统下测试，需要安装Python(3.7.1以上版本)和R(3.4.4以上版本)

### Perl
```
sudo apt install bioperl libtext-trim-perl libsvg-perl
```

### Muscle
`sudo apt install muscle fasttree emboss`

### Interproscan
从ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/下载最新版的InterProScan (64-bit Linux)，并在某一目录下解压缩，例如/opt

```shell
$ wget -c ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.34-73.0/interproscan-5.34-73.0-64-bit.tar.gz
$ cd /opt
$ sudo tar xzvf interproscan-5.34-73.0-64-bit.tar.gz
# 在shell路径添加执行程序所在目录
export PATH=$PATH:/opt/interproscan-5.34-73.0
```

```sh
# 下载BLAST
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz
# 安装
tar xzvf ncbi-blast-2.7.1+-x64-linux.tar.gz
export PATH=$PATH:/opt/ncbi-blast-2.7.1+/bin
```

```bash
wget http://meme-suite.org/meme-software/4.12.0/meme_4.12.0.tar.gz
tar xzvf meme_4.12.0.tar.gz
sudo apt-get install libxslt1-dev libxml2-dev zlib1g-dev openmpi-bin libopenmpi-dev csh 
./configure --prefix=/home/bc/opt/meme --with-url="http://meme.nbcr.net/meme" 
make && make install
export PATH=$PATH:~/opt/meme/bin
```

## 基因鉴定
### 方法：
梨树(_Pyrus bretschneideri_)全基因组（包括基因组序列、蛋白序列、CDS序列和GFF）由 http://peargenome.njau.edu.cn/ 下载。
拟南芥 (_Arabidopsis thaliana_)全基因组由[Phytozome](http://phytozome.jgi.doe.gov/pz/portal.html)下载。 
对于CDS序列和蛋白质序列，由于其包含了不同的剪切本，选取其中一条最长的剪切本。
由[Pfam数据库](http://pfam.xfam.org/family/)下载NAC结构域[PF02365](http://pfam.xfam.org/family/PF02365)。
利用InterPro[1]搜索梨树全基因组蛋白的结构域，提取含有NAC隐马可夫模型（HMM）的基因。

进入~/nac/search目录，并将拟南芥和梨树的蛋白序列拷贝至pep目录下
```sh
mkdir -p ~/nac/search/pep ~/nac/search/scan ~/nac/search/name
```

进入~/nac/search目录，执行interproscan
```sh
cd ~/nac/search
interproscan.sh -i pep/Athaliana.pep -dp -goterms -iprlookup -f tsv -appl Pfam -o scan/Athaliana.scan
interproscan.sh -i pep/Pbretschneideri.pep -dp -goterms -iprlookup -f tsv -appl Pfam -o scan/Pbretschneideri.scan 
search_final.pl pep scan scan_out
```

根据gene_aliases_20130831.txt中拟南芥基因的名称列表，为拟南芥的NAC命名。

```bash
cd ~/nac/search/name
# 根据gene_aliases_20130831.txt中拟南芥基因的名称列表，为拟南芥的NAC命名
wget -c ftp://ftp.arabidopsis.org/home/tair/Genes/gene_aliases_20130831.txt
gene_name_ath.pl gene_aliases_20130831.txt ../scan_out/Athaliana.id ath.id0
# 调整ath.id0中第2列的基因名称，获得ath.id
cut -f 1,2 ath.id0 |sort -k 2 >ath.id
# 提取拟南芥NAC蛋白序列
ex_map_fa.pl ath.id ../scan_out/Athaliana.pep ath.pep
# 以梨树的NAC基因为查询序列，BLAST比对ath.pep
makeblastdb -in ath.pep -dbtype prot
blastp -query ../scan_out/Pbretschneideri.pep -db ath.pep -num_threads 4 -outfmt 6 -max_target_seqs 1 -max_hsps 1 -out blast.out
# 根据BLAST的结果，对梨树NAC基因命名
name_fam_gene.pl blast.out pbr.id pbr_name.info
# 合并ath.id和pbr.id，以备后用
cat ath.id pbr.id >ath_pbr.id
```

```bash
cd ~/nac/search
# 把scan_out目录中所有文件内的基因ID，转换成基因名
map_id_dir.pl name/ath_pbr.id scan_out out
```

## 进化树
### 方法：
提取基因中保守的结构域氨基酸序列，用于多序列联配，所用软件为MUSCLE，参数使用默认值。根据多序列联配的结果，使用FastTree软件构建系统发生树，所用算法为最大使然法（Maximum Likelihood)。同时用贝斯法构建进化树。所用软件为MrBayes 3.2。

MrBayes参数:

lset rates=gamma;  #碱基替换率使用gamma参数
ngen=1000000 #循环次数为为1000000代
samplefreq=100 #每100次循环取样
nchains=4 # 进行4次隐马氏链计算
stopval=0.01 stoprule=yes; #p-value < 0.01后停止运行。
sumt burnin=100; #计算一致树，前100次结果被舍去。

```bash
mkdir -p ~/nac/pbr/tree && cd ~/nac/pbr/tree
muscle -in ~/nac/search/out/Pbretschneideri.domain -out domain.aln
# 多序列比的结果中取出gap数超过50%的列
remove_alig_gap.pl 0.5 domain.aln domain.fasta
# 把domain.fasta转换成Phylip的输入文件，作NJ树outtree
fasta2phy.pl domain.fasta 1.phy
# 转换outtree为NWK个格式
parse_phylip_nwk.pl outtree nj.nwk

# 作BI树
fasta2mrbayes.pl domain.fasta domain.nex
mpirun -np 4 mb domain.nex
nex2nwk.pl domain.nex.con.tre bi.nwk
ex_nwk_id.pl bi.nwk bi.id

# 作ML树
fasttree -wag -gamma -out ml.nwk domain.fasta
```
### 结果
梨树183个NAC基因分为33个亚组（subgroup），亚组根据其成员基因的名称进行的命名 （Figure 2）。

![](img/fig.2.png)
Figure 2: Consensus phylogenetic tree of NAC genes in P. bretschneideri constructed from amino acid multiple sequence alignment using MrBayes software (left). Schematic representations of the conserved motifs of PbNAC genes (middle). Different motifs are highlighted with different colored boxes with numbers 1 to 20. Lines represent protein regions without detected motif. Schematic representations of the exon–intron compositions of PbNAC genes. Exons, represented by boxes. Dashed lines connecting two exons represent an intron.

#### 用R语言的ggtree包作进化树图

```R
library("ggtree")
library("colorspace")
library("Cairo")
	
df <- read.table("group.txt", header = F, sep = "\t")
c2gene <- split(as.character(df[, 2]), as.character(df[, 1]))
c2gene[[1]] <- NULL
	
tree <- read.tree("bi.nwk")
tree <- groupOTU(tree, c2gene)
c2node <- lapply(c2gene, function(x) MRCA(tree, tip=x) )
	
cairo_pdf("ggtree_bi.pdf")
p <- ggtree(tree, aes(color=group), layout="circular", branch.length="none") +  
geom_tiplab(size=2, aes(angle=angle)) + 
geom_text2(aes(label=sprintf("%.2f", as.numeric(label)), subset = !is.na(as.numeric(label)) & as.numeric(label) > 0.88), size=1, hjust=1, vjust=-0.4) + 
xlim(NA, 20) + 
scale_color_manual(values=c("black", rainbow_hcl(length(c2gene))))
	
s <- "p"
cs <- names(c2node)
for (c in cs) {
	s1 <- sprintf("geom_cladelabel(node=%s, label=\"%s\", hjust='center', offset.text=1, offset=4, fontsize=2.5)", c2node[c], c)
	s <- paste(s, s1, sep=" + ")
}
eval(parse(text=s))
dev.off()
```

## 结构域特征

### 方法

```bash
mkdir -p ~/nac/pbr/info && cd ~/nac/pbr/info
meme ~/nac/search/out/Pbretschneideri.pep -p 4 -mod anr -nmotifs 20 -minw 6 -maxw 100 -maxsize 1000000 -text >meme.txt
draw_meme.pl ~/nac/pbr/tree/bi.id ~/nac/search/out/Pbretschneideri.pep meme.txt >meme.svg

map_id.pl ~/nac/search/name/pbr.id ~/nac/search/scan/Pbretschneideri.scan pbr.scan
draw_meme_scan.pl ~/nac/pbr/tree/bi.id 
~/nac/search/out/Pbretschneideri.pep meme.txt pbr.scan >meme_scan.svg

map_id.pl ~/nac/search/name/pbr.id ~/nac/search/scan/Pbretschneideri.scan pbr.scan
map_id.pl ~/nac/search/name/pbr.id pbr.gff3 pbr.gff
draw_exon_phase.pl ~/nac/pbr/tree/bi.id pbr.scan  pbr.gff >exon.svg
```
### 结果
PbNAC的motif结构见图2。基因exon-intron见下图。
## 

```bash
map_id.pl ~/nac/search/name/pbr.id drought.rpkm drought.rpkm.txt
```

```R
library(pheatmap)
data_all <- as.matrix(read.table("drought.rpkm.txt", header=T))
nac_gene <- scan("/home/bc/nac/search/name/pbr.id", what="")
gene_ids <- intersect(nac_gene, rownames(data_all))
expr <- data_all[gene_ids,]
write.table(expr, "drought_expr.txt", quote=F)
```

```bash
awk '{print $1"\t"$1"("$2")"}' ~/nac/pbr/tree/group.txt  >group_gene.id
map_id.pl group_gene.id drought_expr.txt drought_expr_group.txt
```

```R
library(pheatmap)
expr <- read.table("drought_expr_group.txt", header=T)
pheatmap(expr, scale="row", cluster_rows=T, cellwidth = 20, cellheight = 10, fontsize = 8, filename = "drought.svg")
```



```R
	
library("ggtree")
library("colorspace")
library("Cairo")
	
df <- read.table("group.txt", header=T, sep="\t")
c2gene <- split(as.character(df[, 2]), as.character(df[, 1]))
c2gene[[1]] <- NULL
	
tree <- read.tree("ml.nwk")
tree <- groupOTU(tree, c2gene)
c2node <- lapply(c2gene, function(x) MRCA(tree, tip=x) )
	
cairo_pdf("ggtree_ml.pdf", height=150)
p <- ggtree(tree, aes(color=group), branch.length="none") + 
	geom_tiplab() + 
	geom_text2(aes(label=sprintf("%.2f", as.numeric(label)), subset = !is.na(as.numeric(label)) & as.numeric(label) > 0.88), size=2, hjust=1, vjust=-0.4) + 
	xlim(NA, 40) + 
	scale_color_manual(values=c("black", rainbow_hcl(length(c2gene))))
s <- "p"
cs <- names(c2node)
for (c in cs) {
	s1 <- sprintf("geom_cladelabel(node=%s, label=\"%s\",  align=TRUE, offset=6)", c2node[c], c)
	s <- paste(s, s1, sep=" + ")
}
eval(parse(text=s))
dev.off()
```

```
library(ggplot2)

df <- read.table("syn_ks.txt", header=F)
ix <- df$V5 <= 3
df <- df[ix,]

cairo_pdf("ks.pdf")
p <- ggplot(df) + geom_bar(aes(x=V5)) + stat_bin(aes(x=V5), binwidth=0.05, boundary=0, color= "red", fill = "white") + scale_x_continuous("Ks", breaks=seq(0, 3, by=0.5))
p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(),axis.line=element_line(colour="black"))
dev.off()


from pandas import Series, DataFrame
import pandas as pd
df1 = pd.read_table('map.id', header=None)
df2 = pd.read_table('name.info', header=None)
df3 = pd.read_table('chr.info', header=None)
df = pd.merge(df1, df2, on=0)
df = pd.merge(df, df3, on=0)
df.to_csv('gene.info', index=False, header=False, sep='\t')
```



## 共线性检测

多倍化引起的染色体数目的改变是基因组中基因产生多个拷贝的方式之一。 由基因组或基因组片段复制导致的同源基因，不仅表现为基因序列具有相似性，而且染色体上相邻同源基因的排列顺序具有很高的一致性。 这种基因排列的保守性称为共线性（synteny）。考察同一物种内同源基因染色体区段内的共线性同源基因对，包括同源基因对数目，同义置换率（ks），可估算这两个基因间的复制方式和时间。
方法：

假设对基因A和B进行共线性分析，分别提取基因A和B 所在染色体上的一定数量的相邻基因，得到染色体片段A和B。 在这2个染色体片段间的同源基因对的排列如果符合一定的线性关系（由linear regression的相关系数决定），则认为这2个片段间存在共线性关系，继而推测其复制方式为WGD或segmental duplication。 共线性区段中的同源基因对的平均Ks值反映了复制的大体时间。

### 方法
利用软件包 [Zebra](https://github.com/caibinperl/zebra.git)检测基因家族成员间的共线性。
参数：
BLAST E-value:  1e-5   #决定同源基因的BLAST的E值
Number of adjacent genes: 30   #提取NAC基因左右2侧各15个基因
Minimum number of homologous gene pairs in synteny: 3  #在共线性区段内至少存在3对同源基因
E-value: 0.01 #为共线性区段存在的E值，提供统计学上的判断。
Q-value: 0.9  #linear regression的相关系数。

### 结果

以PbNAC73a和PbNAC73b为例，这2个基因间的共线性如图。

![](img/lm.png)

更多关余基因家族共线性的检测结果详见：[Zebra: A tool for detection of microsynteny in a gene family](https://github.com/caibinperl/zebra)