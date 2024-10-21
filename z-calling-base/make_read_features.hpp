#ifndef __MAKE_READ_FEATURES_HPP
#define __MAKE_READ_FEATURES_HPP

#include "../corelib/sam_parser.hpp"

#include <torch/script.h>
#include <vector>

typedef double eval_feature_value_type;
constexpr const int EVAL_KMER_BASE_FEATURES = 6;

struct MolModCallz
{
    int id;
    int strand;
    int offset;
    double prob;
};

void build_phred33_to_prob_table(double table[]);

void build_read_base_features(SAM_Parser& readinfo,
    const double* Phred33ToProbTable,
    std::vector<double>& fwd_qual_probs,
    std::vector<double>& rev_qual_probs,
    std::vector<double>& fwd_ipd_features,
    std::vector<double>& rev_ipd_features,
    std::vector<double>& fwd_pw_features,
    std::vector<double>& rev_pw_features);

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
    std::vector<float>& features);

void dump_kmer_features(const eval_feature_value_type* features, const int kmer_size);

#endif // __MAKE_READ_FEATURES_HPP
