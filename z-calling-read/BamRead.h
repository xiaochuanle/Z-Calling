#include <iostream>
#include <filesystem>
#include "../3rdparty/htslib/include/sam.h"
#include <vector>

struct BamRead {

    BamRead() = delete;

    BamRead(bam1_t *aln, int read_label, int bucketCount);

    bool is_forward;
    bool is_mapped;

    // Basuc alignment section
    std::string query_name;
    std::vector<double> probability_proportions;
    int32_t p_num;
    int label;
    bool is_valid;

    int32_t flag;
    std::string reference_name;
    int32_t mapping_quality;
    int32_t reference_start;
    int32_t reference_end;
    int32_t reference_length;
    int32_t query_alignment_start;
    int32_t query_alignment_end;
    int32_t query_alignment_length;
    std::string cigar_string;
    std::string query_sequence;
    std::string query_qualities;
    std::string reference_seq;
};

std::string get_aux_tag_str(const bam1_t *b, const char tag[2]);

std::vector<double> splitStringToDoubles(const std::string &input);

std::vector<double> calculateProportionInBuckets(const std::vector<double> &numbers, int bucketCount);
