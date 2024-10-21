#include "make_read_features.hpp"
#include "../parameter/parameters.hpp"

#include <cmath>

using namespace std;

static double kBaseOneHotEncodeTable[][4] = {
    { 1.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0, 0.0 },
    { 0.0, 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 0.0, 1.0 }
};

void build_phred33_to_prob_table(double table[])
{
    for (int i = 0; i < 256; ++i) {
        double p = i;
        p -= 33.0;
        p /= -10.0;
        p = 1.0 - pow(10.0, p);
        table[i] = p;
    }
}

/*
def phred_quality_score_to_prob(phred_quality_scores):
    # QUAL = (-10logP) + 33
    phred_quality_scores -= 33.0
    phred_quality_scores /= -10.0
    return 1.0 - np.power(10.0, phred_quality_scores)
*/
static void
s_build_fastq_quality_probs(const SAM_Parser::fq_qual_type* qual, const int seq_size, vector<double>& probs)
{
    probs.clear();
    for (int i = 0; i < seq_size; ++i) {
        double p = qual[i];
        p -= 33.0;
        p /= -10.0;
        p = 1.0 - pow(10.0, p);
        probs.push_back(p);
    }
}

#if 0
static void
s_build_kinetic_features(const int* codev1, const SAM_Parser::kinetics_signal_value_type* v, const int seq_size, vector<double>& features)
{
    features.clear();
    for (int i = 0; i < seq_size; ++i) {
        double f = codev1[ v[i] ];
        hbn_assert (f >= 0 && f <= 952);
        f /= 952.0;
        features.push_back(f);
    }
}
#else
static void
s_build_kinetic_features(const int* codev1, const SAM_Parser::kinetics_signal_value_type* v, const int seq_size, vector<double>& features)
{
    features.clear();
    for (int i = 0; i < seq_size; ++i) {
        double f = v[i];
        f /= 255.0;
        features.push_back(f);
    }
}
#endif 

void build_read_base_features(SAM_Parser& readinfo,
    const double* Phred33ToProbTable,
    std::vector<double>& fwd_qual_probs,
    std::vector<double>& rev_qual_probs,
    std::vector<double>& fwd_ipd_features,
    std::vector<double>& rev_ipd_features,
    std::vector<double>& fwd_pw_features,
    std::vector<double>& rev_pw_features)
{
    const int seq_size = readinfo.query_size();

    fwd_qual_probs.clear();
    const SAM_Parser::fq_qual_type* fwd_qual = readinfo.fwd_qual();
    for (int i = 0; i < seq_size; ++i) {
        fwd_qual_probs.push_back(Phred33ToProbTable[fwd_qual[i]]);
    }
    
    rev_qual_probs.clear();
    const SAM_Parser::fq_qual_type* rev_qual = readinfo.rev_qual();
    for (int i = 0; i < seq_size; ++i) {
        rev_qual_probs.push_back(Phred33ToProbTable[rev_qual[i]]);
    }

    s_build_kinetic_features(readinfo.codev1_table(), readinfo.fwd_ipd(), seq_size, fwd_ipd_features);
    s_build_kinetic_features(readinfo.codev1_table(), readinfo.rev_ipd(), seq_size, rev_ipd_features);

    s_build_kinetic_features(readinfo.codev1_table(), readinfo.fwd_pw(), seq_size, fwd_pw_features);
    s_build_kinetic_features(readinfo.codev1_table(), readinfo.rev_pw(), seq_size, rev_pw_features);
}

void dump_kmer_features(const eval_feature_value_type* features, const int kmer_size)
{
    for (int i = 0; i < kmer_size; ++i) {
        fprintf(stderr, "%d\t", i);
        const eval_feature_value_type* f = features + i * EVAL_KMER_BASE_FEATURES;
        for (int j = 0; j < EVAL_KMER_BASE_FEATURES; ++j) {
            fprintf(stderr, "\t%g", f[j]);
        }
        fprintf(stderr, "\n");
    }
}

void build_kmer_features(const char* fwd_seq,
    const char* rev_seq,
    std::vector<double>& fwd_qual_probs,
    std::vector<double>& rev_qual_probs,
    std::vector<double>& fwd_ipd_features,
    std::vector<double>& rev_ipd_features,
    std::vector<double>& fwd_pw_features,
    std::vector<double>& rev_pw_features,
    const int fn,
    const int rn,
    const int strand,
    const int offset,
    const int seq_size,
    const int flanking_bases,
    std::vector<float>& features)
{
    const char* x_seq;
    const char* y_seq;
    const double* x_qual;
    const double* y_qual;
    const double* x_ipd;
    const double* y_ipd;
    const double* x_pw;
    const double* y_pw;
    int xn, yn;
    if (strand == FWD) {
        x_seq = fwd_seq;
        y_seq = rev_seq;
        x_qual = fwd_qual_probs.data();
        y_qual = rev_qual_probs.data();
        x_ipd = fwd_ipd_features.data();
        y_ipd = rev_ipd_features.data();
        x_pw = fwd_pw_features.data();
        y_pw = rev_pw_features.data();
        xn = fn;
        yn = rn;
    } else {
        y_seq = fwd_seq;
        x_seq = rev_seq;
        y_qual = fwd_qual_probs.data();
        x_qual = rev_qual_probs.data();
        y_ipd = fwd_ipd_features.data();
        x_ipd = rev_ipd_features.data();
        y_pw = fwd_pw_features.data();
        x_pw = rev_pw_features.data();
        xn = rn;
        yn = fn;
    }

    vector<int> pos_list;
    int left_padding = (flanking_bases > offset) ? flanking_bases - offset : 0;
    for (int i = 0; i < left_padding; ++i) pos_list.push_back(0);
    int from = max(0, offset - flanking_bases);
    int to = min(seq_size, offset + 1 + flanking_bases);
    for (int i = from; i < to; ++i) pos_list.push_back(i);
    int right_padding = (offset + 1 + flanking_bases > seq_size) ? (offset + 1 + flanking_bases - seq_size) : 0;
    for (int i = 0; i < right_padding; ++i) pos_list.push_back(seq_size - 1);

    const vector<int>::size_type kmer_size = flanking_bases * 2 + 1;
    hbn_assert(pos_list.size() == kmer_size);
    hbn_assert(x_seq[offset] == REV_MOD_BASE_Z);

    size_t last_feature_offset = features.size();
    for (vector<int>::size_type _i = 0; _i < kmer_size; ++_i) {
        int p = pos_list[_i];
        int fi = 0;

        int c = x_seq[p];
        c = nst_nt4_table[c];
        for (int k = 0; k < 4; ++k, ++fi) features.push_back(kBaseOneHotEncodeTable[c][k]);

        features.push_back(x_ipd[p]);
        ++fi;

        features.push_back(x_pw[p]);
        ++fi;

        hbn_assert(fi == EVAL_KMER_BASE_FEATURES);
    }
}
