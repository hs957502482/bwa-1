#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "bwase.h"
#include "bwtaln.h"
#include "bntseq.h"
#include "utils.h"
#include "kstring.h"
#include "bwa.h"
#include "ksw.h"
#include <pthread.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int g_log_n[256];






void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi)
{
	// 如果没有比对结果
	int i, cnt, best;
	if (n_aln == 0) {
		//fprintf(stderr,"n_aln:%d",s->sa);

		s->type = BWA_TYPE_NO_MATCH;
		s->c1 = s->c2 = 0;
		return;
	}

	if (set_main) {// 如果需要设置主比对结果
		best = aln[0].score; // 设置最佳得分为第一个比对结果的得分

		for (i = cnt = 0; i < n_aln; ++i) {
			const bwt_aln1_t *p = aln + i;

			if (p->score > best) break;
			  // 根据一定概率选择比对结果中的一个作为主要结果
			if (drand48() * (p->l - p->k + 1 + cnt) > (double)cnt) {
				s->n_mm = p->n_mm; s->n_gapo = p->n_gapo; s->n_gape = p->n_gape;// 设置错配、gap open和gap extend的数目
				s->ref_shift = (int)p->n_del - (int)p->n_ins;
				s->score = p->score;
				s->sa = p->k + (bwtint_t)((p->l - p->k + 1) * drand48()); // 随机选择sa区间的一个位置

			}
			cnt += p->l - p->k + 1;//更新比对结果的长度
		}
		s->c1 = cnt;//设置主比对结果的长度
		for (; i < n_aln; ++i) cnt += aln[i].l - aln[i].k + 1; // 计算剩余比对结果的总长度
		s->c2 = cnt - s->c1;// 设置剩余比对结果的总长度
		s->type = s->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE; // 如果主比对结果的长度大于1，则设置序列类型为重复，否则设置为唯一
	}

	if (n_multi) {// 如果需要设置多重比对结果
		int k, rest, n_occ, z = 0;
		for (k = n_occ = 0; k < n_aln; ++k) {
			const bwt_aln1_t *q = aln + k;
			n_occ += q->l - q->k + 1;// 计算所有比对结果的总长度
		}
		if (s->multi) free(s->multi);
		if (n_occ > n_multi + 1) { // if there are too many hits, generate none of them
			s->multi = 0; s->n_multi = 0;// 如果比对结果数量超过给定的最大数量加1，则不生成多重比对结果
			return;
		}
		/* The following code is more flexible than what is required
		 * here. In principle, due to the requirement above, we can
		 * simply output all hits, but the following samples "rest"
		 * number of random hits. */
		rest = n_occ > n_multi + 1? n_multi + 1 : n_occ; // find one additional for ->sa // 计算要生成的多重比对结果数量
		s->multi = calloc(rest, sizeof(bwt_multi1_t));
		for (k = 0; k < n_aln; ++k) {
			const bwt_aln1_t *q = aln + k;
			if (q->l - q->k + 1 <= rest) {// 如果比对结果的长度不超过需要生成的多重比对结果数量
                bwtint_t l;
				for (l = q->k; l <= q->l; ++l) {
					s->multi[z].pos = l;
					s->multi[z].gap = q->n_gapo + q->n_gape;
					s->multi[z].ref_shift = (int)q->n_del - (int)q->n_ins;
					s->multi[z++].mm = q->n_mm;
				}
				rest -= q->l - q->k + 1;
			} else { // Random sampling (http://code.activestate.com/recipes/272884/). In fact, we never come here. 
				int j, i;
				for (j = rest, i = q->l - q->k + 1; j > 0; --j) {
					double p = 1.0, x = drand48();
					while (x < p) p -= p * j / (i--);
					s->multi[z].pos = q->l - i;
					s->multi[z].gap = q->n_gapo + q->n_gape;
					s->multi[z].ref_shift = (int)q->n_del - (int)q->n_ins;
					s->multi[z++].mm = q->n_mm;
				}
				rest = 0;
				break;
			}
		}
		s->n_multi = z;
	}
}

void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s)
{
	bwa_aln2seq_core(n_aln, aln, s, 1, 0);
}

int bwa_approx_mapQ(const bwa_seq_t *p, int mm)
{
	int n;
	if (p->c1 == 0) return 23;
	if (p->c1 > 1) return 0;
	if (p->n_mm == mm) return 25;
	if (p->c2 == 0) return 37;
	n = (p->c2 >= 255)? 255 : p->c2;
	return (23 < g_log_n[n])? 0 : 23 - g_log_n[n];
}

bwtint_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int ref_len, int *strand)
{
	bwtint_t pos_f;
	int is_rev;
	*strand = 0; // initialise strand to 0 otherwise we could return without setting it
	pos_f = bwt_sa(bwt, sapos); // position on the forward-reverse coordinate
	if (pos_f < bns->l_pac && bns->l_pac < pos_f + ref_len) return (bwtint_t)-1;
	pos_f = bns_depos(bns, pos_f, &is_rev); // position on the forward strand; this may be the first base or the last base
	*strand = !is_rev;
	if (is_rev) pos_f = pos_f + 1 < ref_len? 0 : pos_f - ref_len + 1; // position of the first base
	return pos_f; // FIXME: it is possible that pos_f < bns->anns[ref_id].offset
}

/*
 * Derive the actual position in the read from the given suffix array
 * coordinates. Note that the position will be approximate based on
 * whether indels appear in the read and whether calculations are
 * performed from the start or end of the read.
 */
void bwa_cal_pac_pos_core(const bntseq_t *bns, const bwt_t *bwt, bwa_seq_t *seq, const int max_mm, const float fnr)
{
	int max_diff, strand;
	if (seq->type != BWA_TYPE_UNIQUE && seq->type != BWA_TYPE_REPEAT) return;    // 如果序列类型不是唯一比对或重复比对，则直接返回
    // 根据给定的测序错误率和最大允许差异计算近似的比对质量值
	max_diff = fnr > 0.0? bwa_cal_maxdiff(seq->len, BWA_AVG_ERR, fnr) : max_mm;
	seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
	//fprintf(stderr, "%d\n", seq->ref_shift);
	// 根据 SA 值和参考基因组信息计算序列的参考基因组位置
	seq->pos = bwa_sa2pos(bns, bwt, seq->sa, seq->len + seq->ref_shift, &strand);
	seq->strand = strand;
	seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff); // 更新比对质量值
	if (seq->pos == (bwtint_t)-1) seq->type = BWA_TYPE_NO_MATCH; // 如果无法找到参考基因组位置，则将序列类型设置为无匹配
}

void bwa_cal_pac_pos(bwt_t *bwt,int tid, int n_threads,const bntseq_t *bns, const char *prefix, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr)
{
	double time_intv = 0;
	int i, j, strand, n_multi;
	char str[1024];
	//bwt_t *bwt;
	// load forward SA
	// double t_real = realtime();
	// strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
	// strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
	// time_intv = realtime() - t_real;

	for (i = 0; i != n_seqs; ++i) {
		
		bwa_seq_t *p = &seqs[i];
		#ifdef HAVE_PTHREAD
			if (i % n_threads != tid)  
			{
		
			continue;
			}
		#endif
		bwa_cal_pac_pos_core(bns, bwt, p, max_mm, fnr); // 调用核心函数计算序列的参考基因组位置
		for (j = n_multi = 0; j < p->n_multi; ++j) {  // 计算多重比对结果的参考基因组位置
			bwt_multi1_t *q = p->multi + j;
			q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len + q->ref_shift, &strand); // 调用函数将 SA（Suffix Array）值转换为参考基因组位置
			q->strand = strand;
			if (q->pos != p->pos && q->pos != (bwtint_t)-1)
				p->multi[n_multi++] = *q;
		}
		p->n_multi = n_multi;
		//fprintf(stderr,"p->pos:%d\n",p->pos);
	}

	//bwt_destroy(bwt);
	//return time_intv;
}

#define SW_BW 50

bwa_cigar_t *bwa_refine_gapped_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, ubyte_t *seq, int ref_shift, bwtint_t *_rb, int *n_cigar)
{
	bwa_cigar_t *cigar = 0;
	uint32_t *cigar32 = 0;
	ubyte_t *rseq;//参考序列
	int64_t k, rb, re, rlen;// 定义变量用于存储序列的索引范围和长度
	int8_t mat[25];// 存储比对打分矩阵
	int w;//比对窗口大小

	bwa_fill_scmat(1, 3, mat);// 填充比对打分矩阵
	rb = *_rb; re = rb + len + ref_shift;// 获取序列起始和终止位置
	assert(re <= l_pac);// 确保序列位置不超出参考序列的范围

	rseq = bns_get_seq(l_pac, pacseq, rb, re, &rlen);  // 获取参考序列片段

	assert(re - rb == rlen);// 确保参考序列长度与期望长度相同
	w = abs((int)rlen - len) * 1.5; // 计算比对窗口大小，考虑到序列长度的差异

	ksw_global(len, seq, rlen, rseq, 5, mat, 5, 1, SW_BW > w? SW_BW : w, n_cigar, &cigar32); // 进行全局比对
	assert(*n_cigar > 0);// 确保至少有一个 CIGAR 操作
	if ((cigar32[*n_cigar - 1]&0xf) == 1) cigar32[*n_cigar - 1] = (cigar32[*n_cigar - 1]>>4<<4) | 3; // change endding ins to soft clipping
	if ((cigar32[0]&0xf) == 1) cigar32[0] = (cigar32[0]>>4<<4) | 3; // change beginning ins to soft clipping
	if ((cigar32[*n_cigar - 1]&0xf) == 2) --*n_cigar; // delete endding del
	if ((cigar32[0]&0xf) == 2) { // delete beginning del
		*_rb += cigar32[0]>>4;// 更新序列的起始位置
		--*n_cigar;
		memmove(cigar32, cigar32+1, (*n_cigar) * 4);
	}
	cigar = (bwa_cigar_t*)cigar32;// 将结果转换为 bwa_cigar_t 结构，并释放参考序列内存
	for (k = 0; k < *n_cigar; ++k)
		cigar[k] = __cigar_create((cigar32[k]&0xf), (cigar32[k]>>4));
	free(rseq);
	return cigar;// 返回 CIGAR 
}

char *bwa_cal_md1(int n_cigar, bwa_cigar_t *cigar, int len, bwtint_t pos, ubyte_t *seq,
                  bwtint_t l_pac, ubyte_t *pacseq, kstring_t *str, int *_nm)
{
    bwtint_t x, y; // 序列和参考序列上的位置索引
    int z, u, c, nm = 0; // 临时变量和不匹配计数
    str->l = 0; // 重置字符串长度为0
    x = pos; y = 0; // 初始化序列和参考序列的位置索引

    if (cigar) { // 如果存在 CIGAR 操作
        int k, l;
        for (k = u = 0; k < n_cigar; ++k) {
            l = __cigar_len(cigar[k]); // 获取当前 CIGAR 操作的长度
            if (__cigar_op(cigar[k]) == FROM_M) { // 如果是匹配操作
                for (z = 0; z < l && x+z < l_pac; ++z) {
                    c = pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3; // 获取参考碱基
                    if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) { // 如果参考碱基与序列碱基不匹配
                        ksprintf(str, "%d", u); // 添加匹配长度到字符串中
                        kputc("ACGTN"[c], str); // 添加参考碱基到字符串中
                        ++nm; // 不匹配计数加1
                        u = 0; // 重置匹配计数
                    } else ++u; // 如果匹配，匹配计数加1
                }
                x += l; y += l; // 更新序列和参考序列的位置索引
            } else if (__cigar_op(cigar[k]) == FROM_I || __cigar_op(cigar[k]) == FROM_S) { // 如果是插入或软裁剪操作
                y += l; // 更新序列位置索引
                if (__cigar_op(cigar[k]) == FROM_I) nm += l; // 如果是插入操作，不匹配计数加上插入的碱基数
            } else if (__cigar_op(cigar[k]) == FROM_D) { // 如果是删除操作
                ksprintf(str, "%d", u); // 添加匹配长度到字符串中
                kputc('^', str); // 添加删除标记到字符串中
                for (z = 0; z < l && x+z < l_pac; ++z) // 添加删除的碱基到字符串中
                    kputc("ACGT"[pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3], str);
                u = 0; // 重置匹配计数
                x += l; nm += l; // 更新序列位置索引和不匹配计数
            }
        }
    } else { // 如果不存在 CIGAR 操作
        for (z = u = 0; z < (bwtint_t)len && x+z < l_pac; ++z) { // 遍历序列中的碱基
            c = pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3; // 获取参考碱基
            if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) { // 如果参考碱基与序列碱基不匹配
                ksprintf(str, "%d", u); // 添加匹配长度到字符串中
                kputc("ACGTN"[c], str); // 添加参考碱基到字符串中
                ++nm; // 不匹配计数加1
                u = 0; // 重置匹配计数
            } else ++u; // 如果匹配，匹配计数加1
        }
    }
    ksprintf(str, "%d", u); // 添加最后一个匹配长度到字符串中
    *_nm = nm; // 更新不匹配计数
    return strdup(str->s); // 返回 MD 字符串的副本
}

void bwa_correct_trimmed(bwa_seq_t *s)
{
	if (s->len == s->full_len) return;
	if (s->strand == 0) { // forward
		if (s->cigar && __cigar_op(s->cigar[s->n_cigar-1]) == FROM_S) { // the last is S
			s->cigar[s->n_cigar-1] += s->full_len - s->len;
		} else {
			if (s->cigar == 0) {
				s->n_cigar = 2;
				s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
				s->cigar[0] = __cigar_create(0, s->len);
			} else {
				++s->n_cigar;
				s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
			}
			s->cigar[s->n_cigar-1] = __cigar_create(3, (s->full_len - s->len));
		}
	} else { // reverse
		if (s->cigar && __cigar_op(s->cigar[0]) == FROM_S) { // the first is S
			s->cigar[0] += s->full_len - s->len;
		} else {
			if (s->cigar == 0) {
				s->n_cigar = 2;
				s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
				s->cigar[1] = __cigar_create(0, s->len);
			} else {
				++s->n_cigar;
				s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
				memmove(s->cigar + 1, s->cigar, (s->n_cigar-1) * sizeof(bwa_cigar_t));
			}
			s->cigar[0] = __cigar_create(3, (s->full_len - s->len));
		}
	}
	s->len = s->full_len;
}

void bwa_refine_gapped(ubyte_t *pacseq,int tid, int n_threads,const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq)
{
	//ubyte_t *pacseq;
	int i, j, k;
	kstring_t *str;

	// if (!_pacseq) {
	// 	pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
	// 	err_rewind(bns->fp_pac);
	// 	err_fread_noeof(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
	// } else pacseq = _pacseq;

	for (int i = 0; i != n_seqs; ++i) {

		bwa_seq_t *s = seqs + i;
		#ifdef HAVE_PTHREAD
			if (i % n_threads != tid)  
			{
				continue;
			}
		#endif

		seq_reverse(s->len, s->seq, 0); // IMPORTANT: s->seq is reversed here!!! // 反转序列

		for (j = k = 0; j < s->n_multi; ++j) {// 处理多重比对的每一个结果
			bwt_multi1_t *q = s->multi + j;
			int n_cigar;
			if (q->gap) { // gapped alignment
				q->cigar = bwa_refine_gapped_core(bns->l_pac, pacseq, s->len, q->strand? s->rseq : s->seq, q->ref_shift, &q->pos, &n_cigar);
				q->n_cigar = n_cigar;
				if (q->cigar) s->multi[k++] = *q;// 如果生成了 CIGAR，则保留
			} else s->multi[k++] = *q; // 否则保留原结果
		}

		s->n_multi = k; // this squeezes out gapped alignments which failed the CIGAR generation // 压缩掉没有生成 CIGAR 的间隙比对结果
		if (s->type == BWA_TYPE_NO_MATCH || s->type == BWA_TYPE_MATESW || s->n_gapo == 0) continue;
 		 // 处理单一比对结果的间隙
		s->cigar = bwa_refine_gapped_core(bns->l_pac, pacseq, s->len, s->strand? s->rseq : s->seq, s->ref_shift, &s->pos, &s->n_cigar);

		if (s->cigar == 0) s->type = BWA_TYPE_NO_MATCH; // 如果没有生成 CIGAR，则标记为无匹配

	}

	// generate MD tag
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	for (int i = 0; i != n_seqs; ++i) {

		bwa_seq_t *s = seqs + i;
		#ifdef HAVE_PTHREAD
			if (i % n_threads != tid)  
			{
				continue;
			}
		#endif
		if (s->type != BWA_TYPE_NO_MATCH) {
			int nm;
			s->md = bwa_cal_md1(s->n_cigar, s->cigar, s->len, s->pos, s->strand? s->rseq : s->seq, bns->l_pac, pacseq, str, &nm);
			s->nm = nm;
		}
	}
	free(str->s); free(str);

	// correct for trimmed reads
	for (int i = 0; i < n_seqs; ++i) 
	{
		#ifdef HAVE_PTHREAD
			if (i % n_threads != tid)  
			{
				continue;
			}
		#endif
		bwa_correct_trimmed(seqs + i);
	}
	//if (!_pacseq) free(pacseq);
}

int64_t pos_end(const bwa_seq_t *p)
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + p->len;
}

int64_t pos_end_multi(const bwt_multi1_t *p, int len) // analogy to pos_end()
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + len;
}

static int64_t pos_5(const bwa_seq_t *p)
{
	if (p->type != BWA_TYPE_NO_MATCH)
		return p->strand? pos_end(p) : p->pos;
	return -1;
}

void bwa_print_seq(FILE *stream, bwa_seq_t *seq) {
	char buffer[4096];
	const int bsz = sizeof(buffer);
	int i, j, l;
	
	if (seq->strand == 0) {
		for (i = 0; i < seq->full_len; i += bsz) {
			l = seq->full_len - i > bsz ? bsz : seq->full_len - i;
			for (j = 0; j < l; j++) buffer[j] = "ACGTN"[seq->seq[i + j]];
			err_fwrite(buffer, 1, l, stream);
		}
	} else {
		for (i = seq->full_len - 1; i >= 0; i -= bsz) {
			l = i + 1 > bsz ? bsz : i + 1;
			for (j = 0; j < l; j++) buffer[j] = "TGCAN"[seq->seq[i - j]];
			err_fwrite(buffer, 1, l, stream);
		}
	}
}

void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{
	int j;
	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
		int seqid, nn, am = 0, flag = p->extra_flag;
		char XT;

		if (p->type == BWA_TYPE_NO_MATCH) {
			p->pos = mate->pos;
			p->strand = mate->strand;
			flag |= SAM_FSU;
			j = 1;
		} else j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

		// get seqid
		nn = bns_cnt_ambi(bns, p->pos, j, &seqid);
		if (p->type != BWA_TYPE_NO_MATCH && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
			flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences

		// update flag and print it
		if (p->strand) flag |= SAM_FSR;
		if (mate) {
			if (mate->type != BWA_TYPE_NO_MATCH) {
				if (mate->strand) flag |= SAM_FMR;
			} else flag |= SAM_FMU;
		}
		err_printf("%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
		err_printf("%d\t%d\t", (int)(p->pos - bns->anns[seqid].offset + 1), p->mapQ);

		// print CIGAR
		if (p->cigar) {
			for (j = 0; j != p->n_cigar; ++j)
				err_printf("%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
		} else if (p->type == BWA_TYPE_NO_MATCH) err_printf("*");
		else err_printf("%dM", p->len);

		// print mate coordinate
		if (mate && mate->type != BWA_TYPE_NO_MATCH) {
			int m_seqid;
			long long isize;
			am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
			bns_cnt_ambi(bns, mate->pos, mate->len, &m_seqid);
			err_printf("\t%s\t", (seqid == m_seqid)? "=" : bns->anns[m_seqid].name);
			isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
			if (p->type == BWA_TYPE_NO_MATCH) isize = 0;
			err_printf("%d\t%lld\t", (int)(mate->pos - bns->anns[m_seqid].offset + 1), isize);
		} else if (mate) err_printf("\t=\t%d\t0\t", (int)(p->pos - bns->anns[seqid].offset + 1));
		else err_printf("\t*\t0\t0\t");

		// print sequence and quality
		bwa_print_seq(stdout, p);
		err_putchar('\t');
		if (p->qual) {
			if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
			err_printf("%s", p->qual);
		} else err_printf("*");

		if (bwa_rg_id[0]) err_printf("\tRG:Z:%s", bwa_rg_id);
		if (p->bc[0]) err_printf("\tBC:Z:%s", p->bc);
		if (p->clip_len < p->full_len) err_printf("\tXC:i:%d", p->clip_len);
		if (p->type != BWA_TYPE_NO_MATCH) {
			int i;
			// calculate XT tag
			XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';
			// print tags
			err_printf("\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD)? "NM" : "CM", p->nm);
			if (nn) err_printf("\tXN:i:%d", nn);
			if (mate) err_printf("\tSM:i:%d\tAM:i:%d", p->seQ, am);
			if (p->type != BWA_TYPE_MATESW) { // X0 and X1 are not available for this type of alignment
				err_printf("\tX0:i:%d", p->c1);
				if (p->c1 <= max_top2) err_printf("\tX1:i:%d", p->c2);
			}
			err_printf("\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
			if (p->md) err_printf("\tMD:Z:%s", p->md);
			// print multiple hits
			if (p->n_multi) {
				err_printf("\tXA:Z:");
				for (i = 0; i < p->n_multi; ++i) {
					bwt_multi1_t *q = p->multi + i;
					int k;
					j = pos_end_multi(q, p->len) - q->pos;
					nn = bns_cnt_ambi(bns, q->pos, j, &seqid);
					err_printf("%s,%c%d,", bns->anns[seqid].name, q->strand? '-' : '+',
						   (int)(q->pos - bns->anns[seqid].offset + 1));
					if (q->cigar) {
						for (k = 0; k < q->n_cigar; ++k)
							err_printf("%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
					} else err_printf("%dM", p->len);
					err_printf(",%d;", q->gap + q->mm);
				}
			}
		}
		err_putchar('\n');
	} else { // this read has no match
		//ubyte_t *s = p->strand? p->rseq : p->seq;
		int flag = p->extra_flag | SAM_FSU;
		if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;
		err_printf("%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
		//Why did this work differently to the version above??
		//for (j = 0; j != p->len; ++j) putchar("ACGTN"[(int)s[j]]);
		bwa_print_seq(stdout, p);
		err_putchar('\t');
		if (p->qual) {
			if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
			err_printf("%s", p->qual);
		} else err_printf("*");
		if (bwa_rg_id[0]) err_printf("\tRG:Z:%s", bwa_rg_id);
		if (p->bc[0]) err_printf("\tBC:Z:%s", p->bc);
		if (p->clip_len < p->full_len) err_printf("\tXC:i:%d", p->clip_len);
		err_putchar('\n');
	}
}

void bwase_initialize() 
{
	int i;
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
}

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt;
	//bwt_t *bwt2;
	int n_seqs,n_occ;
	bwa_seq_t *seqs,*seqs2;
	bntseq_t *bns;
	const char *prefix;
	const gap_opt_t *opt;
	ubyte_t *pacseq;
	//float fnr;
} thread_aux_t;


static void *worker_se(void *data)
{
	// double t1,t2,t3;
	// double t_aln2seq,t_calpac,t_refine;
	thread_aux_t *d = (thread_aux_t*)data;
	// t1 = realtime();


	// for(int i=0; i < d->n_seqs; ++i)
	// {
	// 	bwa_seq_t *p = (d->seqs) + i;
	// 	#ifdef HAVE_PTHREAD
	// 		if (i % d->opt->n_threads != d->tid)  
	// 		{
		
	// 		continue;
	// 		}
	// 	#endif
		
	// 		bwa_aln2seq_core(p->n_aln, p->aln, p, 1, d->n_occ);
	// }

	bwa_cal_pac_pos(d->bwt,d->tid,d->opt->n_threads,d->bns, d->prefix, d->n_seqs, d->seqs, d->opt->max_diff, d->opt->fnr); // forward bwt will be destroyed here

	bwa_refine_gapped(d->pacseq,d->tid,d->opt->n_threads,d->bns, d->n_seqs, d->seqs ,0);

	
	return 0;
}
#endif



void bwa_sai2sam_se_core(const char *prefix, const char *fn_sa, const char *fn_fa, int n_occ, const char *rg_line)
{
	double ta,tb,tc,td,te,tf=0;

	double t_real1;
	double t_real2;
	double t_real3 = realtime();
	double t_real4;
	double samse_time = 0;
	double read_seq_time = 0;//
	double time_intv = 0;//
	double  restore_bwt_sa_time = 0;//
	double destroy_time = 0;//
	double init_time = 0;
	int times = 0;
	double aln2seq_time = 0;
	extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
	int i, n_seqs, m_aln;
	long long tot_seqs = 0;
	bwt_aln1_t *aln = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bntseq_t *bns;
	FILE *fp_sa;
	gap_opt_t opt1;
	char magic[4];
	// initialization
	t_real1 = realtime();
	bwase_initialize();
	bns = bns_restore(prefix);
	srand48(bns->seed);
	fp_sa = xopen(fn_sa, "r");

	m_aln = 0;
	err_fread_noeof(magic, 1, 4, fp_sa);
	if (strncmp(magic, SAI_MAGIC, 4) != 0) {
		fprintf(stderr, "[E::%s] Unmatched SAI magic. Please re-run `aln' with the same version of bwa.\n", __func__);
		exit(1);
	}

	err_fread_noeof(&opt1, sizeof(gap_opt_t), 1, fp_sa);	

	bwa_print_sam_hdr(bns, rg_line);

	// set ks

	ks = bwa_open_reads(opt1.mode, fn_fa);
	gap_opt_t *opt = &opt1;
	
	
	bwt_t *bwt;
	
	{ // load BWT&SA
		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
		free(str);
		
	}
	opt->n_threads = 4;//线程数为4
	int sa_intv = bwt->sa_intv;
	ubyte_t *pacseq;////存储ref的信息
	pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
	err_rewind(bns->fp_pac);
	err_fread_noeof(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
	// core loop
	init_time = realtime() - t_real1;
	t_real1 = realtime();

	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
		++times;
		t_real2 = realtime();
		tot_seqs += n_seqs;
		t = clock();
		// read alignment
#ifdef HAVE_PTHREAD
	if (opt->n_threads <= 1){
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			int n_aln;
			err_fread_noeof(&n_aln, 4, 1, fp_sa);
			if (n_aln > m_aln) {
				m_aln = n_aln;
				aln = (bwt_aln1_t*)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
			}
			err_fread_noeof(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);
			bwa_aln2seq_core(n_aln, aln, p, 1, n_occ);//处理多重比对等信息
		}
		

		fprintf(stderr, "[bwa_aln_core] convert to sequence coordinate... ");
		bwa_cal_pac_pos(bwt,0,opt->n_threads,bns, prefix, n_seqs, seqs, opt->max_diff, opt->fnr);//计算序列的参考基因组位置 // forward bwt will be destroyed here
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_aln_core] refine gapped alignments... ");
		bwa_refine_gapped(pacseq,0,opt->n_threads,bns, n_seqs, seqs, 0);//利用 ksw_global 函数进行全局比对，然后对比对结果进行处理，最终返回细化后的 CIGAR 结构。
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_aln_core] print alignments... ");
		for (i = 0; i < n_seqs; ++i)
			bwa_print_sam1(bns, seqs + i, 0, opt->mode, opt->max_top2);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}else{//多线程

			
			pthread_t *tid;
			pthread_attr_t attr;
			thread_aux_t *data;
			t_real4 = realtime();
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			int n_aln;
			err_fread_noeof(&n_aln, 4, 1, fp_sa);
			if (n_aln > m_aln) {
				m_aln = n_aln;
				aln = (bwt_aln1_t*)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
			}
			err_fread_noeof(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);
			bwa_aln2seq_core(n_aln, aln, p, 1, n_occ);
		}
			aln2seq_time += realtime() - t_real4; 
			int i,j;
			ta = realtime();
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
			tb = realtime() - ta + tb;
			te = realtime();
			for (j = 0; j < opt->n_threads; ++j) {
				data[j].tid = j; data[j].bwt = bwt;
				data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
				data[j].bns = bns;data[j].pacseq = pacseq;
				data[j].prefix = prefix; data[j].n_occ = n_occ; 
				pthread_create(&tid[j], NULL, worker_se, data + j);
			}
			tf = realtime() - te + tf;

			//bwa_refine_gapped(0,opt->n_threads,bns, n_seqs, seqs2 ,pacseq);
			//fprintf(stderr,"before_join\n");
			tc=realtime();
			for (j = 0; j < opt->n_threads; ++j) 
			{
				pthread_join(tid[j], 0);
			}
			td = realtime()-tc+td;

			for (i = 0; i < n_seqs; ++i) {
				//fprintf(stderr,"bwa_print_sam1:%d\n",tid);
				bwa_print_sam1(bns, seqs + i, 0, opt->mode, opt->max_top2);
			}
			
			free(data);
			free(tid);

	}
#else

		int i;
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			bwa_aln2seq_core(p->n_aln, p->aln, p, 1, n_occ);
		}
		
		bwa_cal_pac_pos(bwt,0,opt->n_threads,bns, prefix, n_seqs, seqs, opt->max_diff, opt->fnr); // forward bwt will be destroyed here
		
		bwa_refine_gapped(0,opt->n_threads,bns,n_seqs, seqs ,0);


#endif
		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[bwa_aln_core] %lld sequences have been processed.\n", tot_seqs);
		time_intv += realtime() - t_real2;
	}
	read_seq_time = realtime() - t_real1 - time_intv;//read_seq
	// destroy
	t_real1 = realtime();
	bwa_seq_close(ks);
	bns_destroy(bns);
	err_fclose(fp_sa);
	free(aln);
	free(pacseq);
	destroy_time = realtime() - t_real1;
	bwt_destroy(bwt);
	fprintf(stderr,"线程创建时间:%f,pthead_join时间:%.3f,pthead_creat时间:%.3f\n",tb,td,tf);
	samse_time =   realtime() - t_real3;
	freopen("../bwa33.txt","a",stderr);
	fprintf(stderr,"occ_intv:%d,sa_intv:%d,read_name:%s,samse_time:%.3f , 串行init_time:%.3f , restore_bwt_sa_time:%.3f , read_seq_time: %.3f ,while内部其他运算:%.3f ,destroy:%.3f,aln2seq:%.3f,循环内串行:%.3f\n",OCC_INTERVAL,sa_intv,fn_fa,samse_time,init_time,restore_bwt_sa_time,read_seq_time ,time_intv , destroy_time,aln2seq_time,read_seq_time+aln2seq_time);
	fprintf(stderr,"线程创建时间:%f,pthead_join时间:%.3f,pthead_creat时间:%.3f\n",tb,td,tf);

	fclose(stderr);

}

int bwa_sai2sam_se(int argc, char *argv[])
{
	int c, n_occ = 3;
	char *prefix, *rg_line = 0;
	while ((c = getopt(argc, argv, "hn:f:r:o:")) >= 0) {
		switch (c) {
		case 'h': break;
		case 'o': OCC_INTV_SHIFT = atoi(optarg);break;
		case 'r':
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
			break;
		case 'n': n_occ = atoi(optarg); break;
		case 'f': xreopen(optarg, "w", stdout); break;
		//case 't': opt->n_threads = atoi(optarg); break;
		default: return 1;
		}
	}
	double tempocc = OCC_INTV_SHIFT;
	OCC_INTERVAL  = (int)pow(2,tempocc);
	OCC_INTV_MASK = OCC_INTERVAL - 1;
	if (optind + 3 > argc) {
		fprintf(stderr, "Usage: bwa samse [-n max_occ] [-f out.sam] [-r RG_line] <prefix> <in.sai> <in.fq>\n");
		return 1;
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[%s] fail to locate the index\n", __func__);
		return 1;
	}
	bwa_sai2sam_se_core(prefix, argv[optind+1], argv[optind+2], n_occ, rg_line);
	free(prefix);
	return 0;
}
