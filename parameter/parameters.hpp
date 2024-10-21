#ifndef __Z_PARAMETERS_HPP
#define __Z_PARAMETERS_HPP

#include "../corelib/sam_parser.hpp"

constexpr const char* const MOTIF_Z = "A";
constexpr const int MOTIF_Z_LEN = 1;
constexpr const char MOD_BASE_Z = 'A';
constexpr const int MOD_CODE_Z = 0;
constexpr const char REV_MOD_BASE_Z = 'T';
constexpr const int REV_MOD_CODE_Z = 3;

static inline bool motif_match_z(const char* seq, const int pos) { return seq[pos] == MOD_BASE_Z; }

static inline int reverse_strand_motif_pos(int offset, int seq_size, int motif_size) {
    return (seq_size - 1) - (offset + motif_size - 1);
}

constexpr const int MIN_READ_SIZE = 1000;
constexpr const int MIN_MAPQ = 10;
constexpr const double MIN_IDENTITY = 98.0;
constexpr const int MIN_STRAND_PASS = 3;
constexpr const double MIN_EC = 6;
constexpr const int KMER_SIZE = 3;
constexpr const bool SIGNAL_IS_ENCODED = true;
constexpr const int NUM_THREADS = 1;
constexpr const double SAMPLE_FRAC = 1.0;

constexpr const int READ_BATCH_SIZE = 10000;
constexpr const int EVAL_SAMPLE_BATCH_SIZE = 10000;

constexpr const int KMER_BASE_FEATURES = 8;

typedef SAM_Parser::fq_qual_type fq_qual_type;
typedef SAM_Parser::kinetics_signal_value_type kinetics_signal_value_type;
typedef kinetics_signal_value_type feature_value_type;
 
#endif // __Z_PARAMETERS_HPP
