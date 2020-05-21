#include <zlib.h>
#include "kseq.h"
#include "bamcat.h"

#include <cassert>
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>       /* ceil */

KSEQ_INIT(gzFile, gzread)

static
int32_t
Nucl2Int(char nucl) {
  switch (nucl) {
    case 'a':
    case 'A':
      return 0;
    case 'c':
    case 'C':
      return 1;
    case 'g':
    case 'G':
      return 2;
    case 't':
    case 'T':
      return 3;
    default:
      std::cerr << "Unexpected symbol " << nucl << std::endl;
      assert(false);
  }
}

//If prefix_len is negative -- go back in the string
static
int32_t
Convert2Int(const char* seq, int32_t prefix_len) {
  assert(prefix_len != 0 && std::abs(prefix_len) < 16);
  int32_t ans = 0;
  if (prefix_len > 0) {
    for (int32_t i = 0; i < prefix_len; ++i) {
      assert(seq[i] != '\0');
      ans = (ans << 2) | Nucl2Int(seq[i]);
    }
  } else {
    for (int32_t i = 0; i < -prefix_len; ++i) {
      ans |= Nucl2Int(*(seq - i)) << (2 * i);
    }
  }
  assert(ans < (1 << (2 * prefix_len)));
  return ans;
}

//Collects kmer stats for the region of length |reg_len|
//If reg_len is negative -- go back in the string
static
void
CollectKmerStat(const char* seq, int32_t reg_len, int32_t kmer_len, int32_t *stats) {
  assert(kmer_len > 0);
  assert(reg_len != 0);
  memset(stats, 0, sizeof(int32_t) * (1 << (2 * kmer_len)));
  if (reg_len > 0) {
    for (int32_t i = 0; (i + kmer_len) <= reg_len; i = i + kmer_len) {
      stats[Convert2Int(seq + i, kmer_len)]++;
    }
  } else {
    for (int32_t i = 0; (i + kmer_len) <= -reg_len; i = i + kmer_len) {
      stats[Convert2Int(seq - i, -kmer_len)]++;
    }
  }
}

bool enable_repeat_check = false;

bool
CheckTrivialDNA(const char* const seq, const uint32_t remaining, const uint32_t offset) {
  if (!enable_repeat_check)
    return false;

  //TODO configure trivial DNA analysis
  static const int32_t SIZE_FACTOR = 6;
  static const int32_t REPEAT_NUM = 5;
  static const int32_t MIN_K = 2;
  static const int32_t MAX_K = 5;

  int32_t stats_buff[1 << (2 * MAX_K)];
  for (int32_t k = MIN_K; k <= MAX_K; ++k) {
    const int32_t possible_kmer_cnt = 1 << (2 * k);
    const int32_t reg_len = k * SIZE_FACTOR;

    //exploring sequence to the right
    for (int32_t shift = 0; shift < k; ++shift) {
      if (reg_len + shift > remaining)
        break;

      CollectKmerStat(seq + shift, reg_len, k, stats_buff);
      if (*std::max_element(stats_buff, stats_buff + possible_kmer_cnt) >= REPEAT_NUM) {
        //char subbuff[reg_len + 1];
        //memcpy(subbuff, seq + shift, reg_len);
        //subbuff[reg_len] = '\0';
        //fprintf(stderr, "Trivial DNA (k=%d) upstream\n", k);
        //fprintf(stderr, "%s\n", subbuff);
        return true;
      }
    }

    //exploring sequence to the left
    for (int32_t shift = 0; shift < k; ++shift) {
      if (reg_len + shift > offset)
        break;
      CollectKmerStat(seq - shift - 1, -reg_len, k, stats_buff);
      if (*std::max_element(stats_buff, stats_buff + possible_kmer_cnt) >= REPEAT_NUM) {
        //char subbuff[reg_len + 1];
        //memcpy(subbuff, seq - shift - reg_len, reg_len);
        //subbuff[reg_len] = '\0';
        //fprintf(stderr, "Trivial DNA (k=%d) downstream\n", k);
        //fprintf(stderr, "%s\n", subbuff);
        return true;
      }
    }
  }
  return false;
}

size_t MaskedPositions(const std::string &s) {
    const char *c_str = s.c_str();
    size_t total_masked = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (CheckTrivialDNA(c_str + i, s.size() - i, i)) {
            ++total_masked;
        }
    }
    std::cout.precision(2);
    std::cout << "Masked total " << total_masked << " (" << double(total_masked / s.size()) << "%) out of " << s.size() << std::endl;
    return total_masked;
}

using namespace std;

/*
  basemap[] works by storing a very small array that maps a base to
  its complement, by dereferencing the array with the ASCII char's
  decimal value as the index
  (int) 'A' = 65;
  (int) 'C' = 67;
  (int) 'G' = 71;
  (int) 'T' = 84;
  (int) 'a' = 97;
  (int) 'c' = 99;
  (int) 'g' = 103;
  (int) 't' = 116;
  (int) 'N' = 78;
  (int) 'U' = 85;
  (int) 'u' = 117;
  for example: basemap['A'] => basemap[65] => 'T' etc.
*/

static const char basemap[255] =
    {
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
        '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
        '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
        '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', /*  90 -  99 */
        '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
        '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a', '\0', '\0', /* 110 - 119 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
        '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
    };

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core) {

    if (core->flag & BAM_FUNMAP) {
        return false;
    }
    if (core->flag & BAM_FSECONDARY) {
        return false;
    }

    return true;
}

char *rc(char *input, int len) {
   char *seq = new char[len+1];
   input+=len-1;
   for (int i = 0; i < len; i++) {
      seq[i]=basemap[*input];
      input--;
   }
   return seq;
}

samfile_t * open_alignment_file(std::string path) {
    samfile_t * fp = NULL;
    std::string flag = "r";
    if (path.substr(path.size() - 3).compare("bam") == 0) {
        //BAM file!
        flag += "b";
    }
    if ((fp = samopen(path.c_str(), flag.c_str(), 0)) == 0) {
        fprintf(stderr, "qaCompute: Failed to open file %s\n", path.c_str());
    }
    return fp;
}

int
main (int argc, char **argv) {
    if (argc < 3 || (argc > 3 && std::string(argv[3]) != "--ignore-microsatellites")) {
        fprintf(stderr, "Usage: %s\n <alignment file (SAM/BAM)> <reference in FASTA> [--ignore-microsatellites]", argv[0]);
    }

    if (argc > 3) {
        enable_repeat_check = true;
        cerr << "Microsatellite repeats around alignment differences will be checked" << endl;
    } else {
        cerr << "No microsatellite repeats checked" << endl;
    }

    map<string, string> reference;
    // loadFasta
    FILE *in = fopen(argv[2], "r");
    gzFile f = gzdopen(fileno(in), "r");
    kseq_t *seq = kseq_init(f);
    int l = 0;
    while ((l=kseq_read(seq)) >= 0) {
        reference[seq->name.s] = seq->seq.s;
        //MaskedPositions(std::string(seq->seq.s));
        const string name = seq->name.s;
        reference[name + "rc"] = rc(seq->seq.s, strlen(seq->seq.s));
    }
    kseq_destroy(seq);
    gzclose(f);
    fclose(in);

    bam1_t *b = bam_init1();
    samfile_t * fp = open_alignment_file(argv[1]);
    if (fp == NULL) {
        exit(1);
    }
    bam_header_t* head = fp->header; // sam header

    //cout << id << "\t" << ref << "\t" << int(ceil(-1*idy/100*len)) << "\t" << idy << "\t0\t" << (isFwd == true ? seqLow : seqLen-seqHi) << "\t" << (isFwd == true ? seqHi : seqLen-seqLow) << "\t" << seqLen << "\t" <<  (isFwd == true ? 0 : 1) << "\t" << (isFwd == true ? refLo : refLen-refHigh) << "\t" << (isFwd == true ? refHigh : refLen-refLo) << "\t" << refLen << "\t" << len << "\t" << (seqHi-seqLow) << "\t" << (refHigh-refLo) << "\t" << matches << "\t" << errors << "\t" << idyMismatches << "\t" << indels << endl;
    cout << "read_id\tref_id\t???\tidentity\t0!\tquery_strand_start\tquery_strand_end\tquery_len\treverse_strand\tref_strand_start\tref_strand_end\tref_len\talignment_len\tquery_span\tref_span\tmatches\terrors\tmismatch_identity\tindels" << endl;

    while (samread(fp, b) >= 0) {
        //Get bam core.
        const bam1_core_t *core = &b->core;
        if (core == NULL) {
            printf("Input file is corrupt!");
            exit(1);
        }
        if (!is_mapped(core)) {
            continue;
        }

        const string id = bam1_qname(b);
        //cerr << "Processing " << id << endl;
        const string ref = head->target_name[core->tid];
        const uint32_t* cigar = bam1_cigar(b);
        const uint32_t refLen = head->target_len[core->tid];
        const uint32_t refLo = core->pos + 1;
        const uint32_t refHigh = bam_calend(core, cigar);
        assert(refLo <= refHigh && refHigh <= refLen);
        uint32_t seqLow = 0;
        uint32_t seqHi = core->l_qseq;
        uint32_t seqLen = core->l_qseq;
        const bool isFwd = !(core->flag & BAM_FREVERSE);

        // now parse, we need both sequences for this
        if (strlen(reference[ref].c_str()) <= 0) {
            fprintf(stderr, "Error: unknown reference sequence %s", ref.c_str());
            exit(1);
        }
        assert(refLen == reference[ref].size());
        const char * const refSeqStart = reference[ref].c_str();
        const char *refSeq = refSeqStart+refLo-1;
        char *qrySeq = new char[core->l_qseq+1];
        char const *qrySeqCleanup = qrySeq;
        for (int i = 0; i < core->l_qseq; i++) {
            qrySeq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)];
        }
        qrySeq[core->l_qseq]='\0';
        int errors = 0;
        int matches = 0;
        int indels = 0;
        uint32_t len = 0;
        for (int k = 0; k < core->n_cigar; ++k) {
            assert(refSeq >= refSeqStart);
            uint32_t refPos(refSeq - refSeqStart);
            assert(refPos <= refLen);
            assert(refPos <= refHigh);

            int cop = cigar[k] & BAM_CIGAR_MASK; // operation
            int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
            switch (cop) {
            // we don't care about clipping, not part of length

            case BAM_CSOFT_CLIP:
                if (k == 0) seqLow+=cl;
                else seqHi-=cl;
                qrySeq+=cl;
                break;

            case BAM_CHARD_CLIP:
                if (k == 0) seqLow+=cl;
                else seqHi-=cl;
                seqLen+=cl;
                seqHi+=cl;
                break;

            // we don't care about matches either, not an error
            case BAM_CMATCH:

            case BAM_CEQUAL:

            case BAM_CDIFF:
                for (int i = 0; i < cl; i++) {
                    assert(refSeq < refSeqStart + refLen);
                    refPos = refSeq - refSeqStart;
                    if (toupper(*refSeq) == toupper(*qrySeq)) {
                        matches++;
                    } else {
                        if (!CheckTrivialDNA(refSeq, refLen - refPos, refPos)) {
                            errors++;
                        }
                    }
                    refSeq++;
                    qrySeq++;
                }
                len+=cl;
                break;

            case BAM_CINS:
                if (!CheckTrivialDNA(refSeq, refLen - refPos, refPos)) {
                    errors+=cl;
                    indels+=cl;
                }
                qrySeq+=cl;
                len+=cl;
                break;

            case BAM_CDEL:
                if (!CheckTrivialDNA(refSeq, refLen - refPos, refPos)) {
                    errors+=cl;
                    indels+=cl;
                }
                refSeq+=cl;
                len+=cl;
                break;

            case BAM_CREF_SKIP:
                refSeq+=cl;
                len+=cl;
                break;
            default:
                cerr << "Error: unknown base " << cop << " of len " << cl << endl;
                exit(1);
                break;
            }
        }

        double idy = (1-((double) errors / len /*(seqHi-seqLow)*/ /*(refHigh - refLo + 1)*/)) * 100;
        double idyMismatches =  (1-((double) (errors-indels) / len /*(seqHi-seqLow)*/ /*(refHigh - refLo + 1)*/)) * 100;
        cout << id << "\t" << ref << "\t" << int(ceil(-1*idy/100*len)) << "\t" << idy << "\t0\t" << (isFwd == true ? seqLow : seqLen-seqHi) << "\t" << (isFwd == true ? seqHi : seqLen-seqLow) << "\t" << seqLen << "\t" <<  (isFwd == true ? 0 : 1) << "\t" << (isFwd == true ? refLo : refLen-refHigh) << "\t" << (isFwd == true ? refHigh : refLen-refLo) << "\t" << refLen << "\t" << len << "\t" << (seqHi-seqLow) << "\t" << (refHigh-refLo) << "\t" << matches << "\t" << errors << "\t" << idyMismatches << "\t" << indels << endl;
        delete[] qrySeqCleanup;
    }
    samclose(fp);
    return 0;
}
