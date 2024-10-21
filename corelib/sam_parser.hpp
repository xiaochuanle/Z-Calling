#ifndef __SAM_PARSER_HPP
#define __SAM_PARSER_HPP

#include "../3rdparty/htslib/include/sam.h"
#include "hbn_aux.h"
#include "seq_name2id_map.hpp"
#include "unpacked_seqdb.hpp"

#include <string>
#include <vector>

class SAM_Parser 
{
public:
    typedef u8      fq_qual_type;
    typedef u8      kinetics_signal_value_type;

public:
    SAM_Parser();

    ~SAM_Parser() {
        ks_free(&M_bam_aux);
    }

    bool parse(sam_hdr_t* hdr, bam1_t* bam, const int read_id, const int min_query_size, HbnUnpackedDatabase* ref = nullptr, SeqName2IdMap* ref_name2id = nullptr);

public:
    int query_id() const {
        return M_query_id;
    }

    const char* query_name() const {
        return M_query_name.c_str();
    }

    int query_strand() const {
        return M_query_strand;
    }

    int query_size() const {
        return M_query_size;
    }

    const char* fwd_seq() const {
        return M_fwd_seq.c_str();
    }

    const char* rev_seq() const {
        return M_rev_seq.c_str();
    }

    const fq_qual_type* fwd_qual() const {
        return M_fwd_qual.data();
    }

    const fq_qual_type* rev_qual() const {
        return M_rev_qual.data();
    }

    const kinetics_signal_value_type* fwd_ipd() const {
        return M_fwd_ipd.data();
    }

    const kinetics_signal_value_type* rev_ipd() const {
        return M_rev_ipd.data();
    }

    const kinetics_signal_value_type* fwd_pw() const {
        return M_fwd_pw.data();
    }

    const kinetics_signal_value_type* rev_pw() const {
        return M_rev_pw.data();
    }

    int fwd_pass() const {
        return M_fwd_pass;
    }

    int rev_pass() const {
        return M_rev_pass;
    }

    double effective_coverage() const {
        return M_effective_coverage;
    }

    bool is_mapped() const {
        return M_is_mapped;
    }

    int qb() const {
        return M_qb;
    }

    int qe() const {
        return M_qe;
    }

    int sid() const {
        return M_sid;
    }

    const char* sname() const {
        return M_sname.c_str();
    }

    int sb() const {
        return M_sb;
    }

    int se() const {
        return M_se;
    }

    int ssize() const {
        return M_ss;
    }

    int mapQ() const {
        return M_mapQ;
    }

    double identity() const {
        return M_identity;
    }

    const char* qas() const {
        return M_qas.c_str();
    }

    const char* sas() const {
        return M_sas.c_str();
    }

    int as_size() const {
        return M_as_size;
    }

    const int* qas_pos_list() const {
        return M_qas_pos_list.data();
    }

    const int* sas_pos_list() const {
        return M_sas_pos_list.data();
    }

    void dump_map_info() const {
        if (!M_is_mapped) return;
        fprintf(stderr, "[%d, %d, %d, %d, %d]", M_query_id, M_query_strand, M_qb, M_qe, M_query_size);
        fprintf(stderr, " x ");
        fprintf(stderr, "[%d, %d, %d, %d]", M_sid, M_sb, M_se, M_ss);
        fprintf(stderr, ", %d, %g", M_mapQ, M_identity);
        fprintf(stderr, ", %s", M_query_name.c_str());
        fprintf(stderr, ", %s", M_sname.c_str());
        fprintf(stderr, "\n");
    }

    void dump() const {
        fprintf(stderr, "query name:     %s\n", M_query_name.c_str());
        fprintf(stderr, "query strand:   %d\n", M_query_strand);
        fprintf(stderr, "query size:     %d\n", M_query_size);
        fprintf(stderr, "fwd-pass:       %d\n", M_fwd_pass);
        fprintf(stderr, "rev-pass:       %d\n", M_rev_pass);

        const char* query_seq = (M_query_strand == FWD) ? fwd_seq() : rev_seq();
        fprintf(stderr, "query:          ");
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%c", query_seq[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, "%c", query_seq[i]);
        fprintf(stderr, "\n");

        const char* rev_query_seq = (M_query_strand == FWD) ? rev_seq() : fwd_seq();
        fprintf(stderr, "query:          ");
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%c", rev_query_seq[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, "%c", rev_query_seq[i]);
        fprintf(stderr, "\n");

        fprintf(stderr, "quality:        ");
        const fq_qual_type* qual = (M_query_strand == FWD) ? fwd_qual() : rev_qual();
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%c", qual[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, "%c", qual[i]);
        fprintf(stderr, "\n");

        const kinetics_signal_value_type* fwd_ipds = fwd_ipd();
        fprintf(stderr, "fwd-ipd:        ");
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%d,", fwd_ipds[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, ",%d", fwd_ipds[i]);
        fprintf(stderr, "\n");

        const kinetics_signal_value_type* rev_ipds = rev_ipd();
        fprintf(stderr, "rev-ipd:        ");
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%d,", rev_ipds[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, ",%d", rev_ipds[i]);
        fprintf(stderr, "\n");

        const kinetics_signal_value_type* fwd_pws = fwd_pw();
        fprintf(stderr, "fwd-pw:         ");
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%d,", fwd_pws[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, ",%d", fwd_pws[i]);
        fprintf(stderr, "\n");

        const kinetics_signal_value_type* rev_pws = rev_pw();
        fprintf(stderr, "rev-pw:         ");
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%d,", rev_pws[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, ",%d", rev_pws[i]);
        fprintf(stderr, "\n");

        dump_map_info();
    }

    const int* codev1_table() const {
        return M_codev1_table;
    }

private:
    void x_clear_kinetics_info() {
        M_query_id = -1;
        M_query_strand = F_R;
        M_query_size = 0;
        M_fwd_seq.clear();
        M_rev_seq.clear();
        M_fwd_qual.clear();
        M_rev_qual.clear();
        M_fwd_ipd.clear();
        M_rev_ipd.clear();
        M_fwd_pw.clear();
        M_rev_pw.clear();
        M_fwd_pass = 0;
        M_rev_pass = 0;
    }

    void x_clear_map_info() {
        M_is_mapped = false;
        M_sid = -1;
        M_sname.clear();
        M_mapQ = 0;
        M_identity = 0.0;
        M_qas.clear();
        M_sas.clear();
        M_as_size = 0;
        M_qas_pos_list.clear();
        M_sas_pos_list.clear();
    }

private:
    kstring_t                                   M_bam_aux;
    std::vector<std::pair<const char*, int>>    M_aux_cols;
    std::vector<std::pair<char, int>>           M_cigar_op_list;

    std::string                                 M_query_name;
    int                                         M_query_id;
    int                                         M_query_strand;
    int                                         M_query_size;
    std::string                                 M_fwd_seq;
    std::string                                 M_rev_seq;
    std::vector<fq_qual_type>                   M_fwd_qual;
    std::vector<fq_qual_type>                   M_rev_qual;
    std::vector<kinetics_signal_value_type>     M_fwd_ipd;
    std::vector<kinetics_signal_value_type>     M_rev_ipd;
    std::vector<kinetics_signal_value_type>     M_fwd_pw;
    std::vector<kinetics_signal_value_type>     M_rev_pw;
    int                                         M_fwd_pass;
    int                                         M_rev_pass;
    double                                      M_effective_coverage;

    /// mapping info
    bool                                        M_is_mapped;
    int                                         M_qb;
    int                                         M_qe;
    int                                         M_sid;
    std::string                                 M_sname;
    int                                         M_sb;
    int                                         M_se;
    int                                         M_ss;
    int                                         M_mapQ;
    double                                      M_identity;
    std::string                                 M_qas;
    std::string                                 M_sas;
    int                                         M_as_size;
    std::vector<int>                            M_qas_pos_list;
    std::vector<int>                            M_sas_pos_list;

    static const int                            kCodeV1TableSize = 256;
    int                                         M_codev1_table[kCodeV1TableSize];
};

#endif // __SAM_PARSER_HPP
