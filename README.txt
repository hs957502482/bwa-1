生成索引
bwa index  -S 5 -C 7  ref.fasta           -S 5代表sa_interval=2^5  -C 7代表occ_interval=2^7
ALN步骤 
bwa aln -C 7 -t 8 ref.fastq read.fastq > result.sai     -C 7代表occ_interval=2^7  -t 8表示8线程
SAMSE步骤
bwa samse -C 7 -t 8 ef.fastq  result.sai  read.fastq > result.sam   -C 7代表occ_interval=2^7  -t 8表示8线程

相关时间信息会记录在数据所在文件夹的上级文件夹中的bwa33.txt文件中。
ALN步骤和SAMSE步骤的occ_interval要与索引的一致，sa_interval只需生成索引时指定，后续比对不需指定。

其他使用方式与BWA相同，注意的是比对时SA_INTV和OCC_INTV要与生成索引时的一致。