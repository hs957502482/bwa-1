SAMSE步骤加入了多线程，需手动在bwase.c中的 bwa_sai2sam_se_core函数中更改opt->n_threads。
在aln和samse步骤中加入了-o 参数，为 OCC_INTV_SHIFT的值，OCC_INTV=2^ OCC_INTV_SHIFT。
SA_INTV的值需手动更改，在bwtindex.c中的bwa_idx_build中更改bwt_cal_sa(bwt, 32)，sa_intv默认为32。

其他使用方式与BWA相同，注意的是比对时SA_INTV和OCC_INTV要与生成索引时的一致。