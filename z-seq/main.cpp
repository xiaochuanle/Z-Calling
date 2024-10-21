#include <iostream>
#include <filesystem>
#include "../3rdparty/argparse/argparse.h"
#include "../3rdparty/htslib/include/sam.h"
#include <fstream>

char fourbits2base(uint8_t val);

std::string get_aux_tag_str(const bam1_t *b, const char tag[2]);

std::vector<double> splitStringToDoubles(const std::string &input);


int main(int argc, char **argv) {
    argparse::ArgumentParser program("z-seq", "0.1.0");
    program.add_argument("input_bam_path");
    program.add_argument("output_file_path");
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    std::filesystem::path input_bam_path = program.get<std::string>("input_bam_path");
    std::filesystem::path output_bam_path = program.get<std::string>("output_file_path");
    samFile *bam_in = sam_open(input_bam_path.c_str(), "r");
    bam_hdr_t *bam_header = sam_hdr_read(bam_in);
    bam1_t *aln = bam_init1();
    std::ofstream out;
    out.open(output_bam_path);
    while (sam_read1(bam_in, bam_header, aln) >= 0) {
        std::string query_name = bam_get_qname(aln);
        std::string zp = get_aux_tag_str(aln, "ZP");
        if (zp.size() < 8) {
            continue;
        }
        out << ">" << query_name << std::endl;
        zp = zp.substr(7);
        std::vector<double> zp_nums = splitStringToDoubles(zp);
        uint8_t *data = bam_get_seq(aln);
        int len = aln->core.l_qseq;
        std::string s(len, '\0');
        int j = 0;
        for (int i = 0; i < len; i++) {
            char base;
            if (i % 2 == 1)
                base = fourbits2base(data[i / 2] & 0xF);
            else
                base = fourbits2base((data[i / 2] >> 4) & 0xF);
            if (base == 'G' || base == 'C') {
                s[i] = base;
            } else if (base == 'A') {
                double prob = zp_nums[j];
                if (prob >= 0.5) {
                    s[i] = 'Z';
                } else {
                    s[i] = 'A';
                }
                j++;
            } else if (base == 'T') {
                double prob = zp_nums[j];
                if (prob >= 0.5) {
                    s[i] = 'O';
                } else {
                    s[i] = 'T';
                }
                j++;
            }
        }
        out << s << std::endl;
    }
    out.close();
    return 0;
}

char fourbits2base(uint8_t val) {
    switch (val) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            std::cerr << "ERROR: Wrong base with value " << (int) val << std::endl;
            return 'N';
    }
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

std::vector<double> splitStringToDoubles(const std::string &input) {
    std::vector<double> result;
    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ',')) {
        result.push_back(std::stod(token));
    }

    return result;
}
