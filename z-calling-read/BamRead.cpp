#include "BamRead.h"
#include <sstream>
#include <vector>
#include <string>
#include <ctime>


BamRead::BamRead(bam1_t *aln, int read_label, int bucketCount)
        : reference_start(-1), reference_end(-1), reference_length(-1), query_alignment_start(-1),
          query_alignment_end(-1), query_alignment_length(-1) {

    query_name = bam_get_qname(aln);
    std::string zp = get_aux_tag_str(aln, "ZP");
    if (zp.size() < 8) {
        is_valid = false;
    } else {
        is_valid = true;
        zp = zp.substr(7);
        std::vector<double> numbers = splitStringToDoubles(zp);
        probability_proportions = calculateProportionInBuckets(numbers, bucketCount);
        p_num = numbers.size();
        label = read_label;
    }
}


std::vector<double> splitStringToDoubles(const std::string &input) {
    std::vector<double> result;
    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ',')) {
        result.push_back(std::stod(token));
    }

    return result;
}


std::vector<double> calculateProportionInBuckets(const std::vector<double> &numbers, int bucketCount) {
    std::vector<int> buckets(bucketCount, 0);
    double bucketWidth = 1.0 / bucketCount;

    for (double num: numbers) {
        if (num >= 0 && num <= 1.0) {
            int bucketIndex = static_cast<int>(num / bucketWidth);
            if (bucketIndex == bucketCount) {
                bucketIndex = bucketCount - 1;
            }
            buckets[bucketIndex]++;
        }
    }

    std::vector<double> proportions(bucketCount, 0.0);
    int totalNumbers = numbers.size();
    for (int i = 0; i < bucketCount; ++i) {
        proportions[i] = static_cast<double>(buckets[i]) / totalNumbers;
    }
    return proportions;
}


std::string get_aux_tag_str(const bam1_t *b, const char tag[2]) {
    kstring_t res = KS_INITIALIZE;
    if (bam_aux_get_str(b, tag, &res) == 1)
    {
        int len = ks_len(&res);
        char *ks_s = ks_str(&res);
        std::string s(ks_s, ks_s + len);
        ks_free(&res);
        return s;
    } else {
        ks_free(&res);
        return "";
    }
}

