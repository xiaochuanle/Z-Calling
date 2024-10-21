#include "sam_parser.hpp"

#include "split_string_by_char.hpp"

using namespace std;

static const int kMinRawSignalValue = 0;
static const int kMaxRawSignalValue = 952;
static const int kMinEncodedSignalValue = 0;
static const int kMaxEncodedSignalValue = 256;

static void
validate_signal_value(bool signal_is_encoded, int& s)
{
    if (signal_is_encoded) {
        if (s < kMinEncodedSignalValue || s >= kMaxEncodedSignalValue) {
            fprintf(stderr, "ERROR: In encoded signal setting, signal value %d is out of plausible range [%d, %d)\n", s, kMinEncodedSignalValue, kMaxEncodedSignalValue);
            abort();
        }
    } else {
	    if (s < kMinRawSignalValue) {
            fprintf(stderr, "ERROR: In raw signal setting, signal value must be >= 0: %d\n", s);
            abort();
	    }
    }
}

/*
def fill_codev1_to_frame_table():

    CodeV1ToFrameTableSize = 256
    codev1_to_frame_table = [0] * CodeV1ToFrameTableSize
    idx = 0

    for i in range(64):
        codev1_to_frame_table[idx] = i 
        idx += 1

    for i in range(64, 128):
        codev1_to_frame_table[idx] = (i - 64) * 2 + 64
        idx += 1

    for i in range(128, 192):
        codev1_to_frame_table[idx] = (i - 128) * 4 + 192
        idx += 1

    for i in range(192, 256):
        codev1_to_frame_table[idx] = (i - 192) * 8 + 448
        idx += 1

    return codev1_to_frame_table
*/

static int s_encode_signal_value(int s)
{
    s = min(952, s);
    int t = 0;
    /// s = (t - 192) * 8 + 448
    /// t = (s - 448) / 8 + 192
    if (s >= 448) {
        t = (s - 448) / 8 + 192;
    }
    /// s = (t - 128) * 4 + 192
    /// t = (s - 192) / 4 + 128
    else if (s >= 192) {
        t = (s - 192) / 4 + 128;
    }
    /// s = (t - 64) * 2 + 64
    /// t = (s - 64) / 2 + 64
    else if (s >= 64) {
        t = (s - 64) / 2 + 64;
    }
    else {
        t = s;
    }
    return t;
}

#if 0
CIGAR operation

OP      Description                                                 consume-query           consume-reference

M       alignment match (can be match or mismatch)                  yes                     yes
I       insertion to the reference                                  yes                     no 
D       deletion from the reference                                 no                      yes 
N       skipped region from the reference                           no                      yes
S       soft clipping (clipped sequences represent in SEQ)          yes                     no 
H       hard clipping (clipped sequences NOT represent in SEQ)      no                      no 
P       padding (silent deletion from padded reference)             no                      no 
=       sequence match                                              yes                     yes
X       sequence mismatch                                           yes                     yes
#endif 

static void cigar_to_alignment(const char* query,
    const int query_size,
    const u8* subject,
    const int subject_size,
    vector<pair<char, int>>& op_list,
    int& qb,
    int& qe,
    int& sb,
    int& se,
    std::string& qas,
    std::string& sas,
    std::vector<int>& q_pos_list,
    std::vector<int>& s_pos_list)
{
    pair<char, int>* opa = op_list.data();
    int opc = op_list.size();
    int opi = 0;
    qb = 0;
    sb = 0;
    if (opa[0].first == 'S') {
        qb = opa[0].second;
        opi = 1;
    } else if (opa[0].first == 'H') {
        opi = 1;
    }
    int qi = qb - 1;
    int si = sb - 1;

    qas.clear();
    sas.clear();
    q_pos_list.clear();
    s_pos_list.clear();
    for (; opi < opc; ++opi) {
        char op = opa[opi].first;
        int num = opa[opi].second;

        if (op == 'M') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += DECODE_RESIDUE(subject[si]);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'I') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                qas += query[qi];
                sas += GAP_CHAR;
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'D') {
            for (int i = 0; i < num; ++i) {
                ++si;
                qas += GAP_CHAR;
                sas += DECODE_RESIDUE(subject[si]);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'N') {
            for (int i = 0; i < num; ++i) {
                ++si;
                qas += GAP_CHAR;
                sas += DECODE_RESIDUE(subject[si]);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }            
        } else if (op == 'S') {
            continue;
        } else if (op == 'H') {
            continue;
        } else if (op == 'P') {
            continue;
        } else if (op == '=') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += DECODE_RESIDUE(subject[si]);
		        int qc = query[qi];
		        int sc = DECODE_RESIDUE(subject[si]);
		        if (qc != sc) {
			        fprintf(stderr, "qb = %d, sb = %d, qi = %d, si = %d, qc = %c, sc = %c\n", qb, sb, qi, si, qc, sc);
			        //cerr << NStr::CTempString(cigar, cigar_size) << '\n';
		        }
		        //hbn_assert(qc == sc);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }            
        } else if (op == 'X') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += DECODE_RESIDUE(subject[si]);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else {
            fprintf(stderr, "ERROR: Unrecognised CIGAR operation '%c' in\n", op);
            //cerr << NStr::CTempString(cigar, cigar_size) << '\n';
            exit(1);
        }
    }
    qe = qi;
    se = si;
}

double
calc_effective_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size)
{
	const int E = 20;
	int eff_len = 0;
	int eff_mat = 0;
	int i = 0;
	while (i < align_size) {
		int qc = query_mapped_string[i];
		int tc = target_mapped_string[i];

		if (qc != GAP_CHAR && tc != GAP_CHAR) {
			if (qc == tc) ++eff_mat;
			++eff_len;
			++i;
			continue;
		}

		if (qc == GAP_CHAR && tc == GAP_CHAR) {
			++i;
			continue;
		}

		if (qc == GAP_CHAR) {
			int j = i + 1;
			while (j < align_size) {
				if (query_mapped_string[j] == GAP_CHAR && target_mapped_string[j] == GAP_CHAR) {
					++j;
					continue;
				}
				if (query_mapped_string[j] != GAP_CHAR) break;
				++j;
			}
			if (j - i < E) {
				for (int k = i; k < j; ++k) {
					int qc = query_mapped_string[k];
					int tc = target_mapped_string[k];
					if (qc == GAP_CHAR && tc == GAP_CHAR) continue;
					if (qc == tc) ++eff_mat;
					++eff_len;
				}
			}
			i = j;
			continue;
		}

		hbn_assert(tc == GAP_CHAR);
		if (tc == GAP_CHAR) {
			int j = i + 1;
			while (j < align_size) {
				if (query_mapped_string[j] == GAP_CHAR && target_mapped_string[j] == GAP_CHAR) {
					++j;
					continue;
				}
				if (target_mapped_string[j] != GAP_CHAR) break;
				++j;
			}
			if (j - i < E) {
				for (int k = i; k < j; ++k) {
					qc = query_mapped_string[k];
					tc = target_mapped_string[k];
					if (qc == GAP_CHAR && tc == GAP_CHAR) continue;
					if (qc == tc) ++eff_mat;
					++eff_len;
				}
			}
			i = j;
			continue;
		}
	}
	if (eff_len == 0) return 0.0;
	return 100.0 * eff_mat / eff_len;
}

#define dump_align_string(qaln, saln, aln_size, stream) do { \
	hbn_fwrite(qaln, 1, aln_size, stream); \
	fprintf(stream, "\n"); \
	for (int _i = 0; _i < aln_size; ++_i) { \
		if (qaln[_i] == saln[_i]) fprintf(stream, "|"); \
		else fprintf(stream, "*"); \
	} \
	fprintf(stream, "\n"); \
	hbn_fwrite(saln, 1, aln_size, stream); \
	fprintf(stream, "\n"); \
} while (0)

void
validate_aligned_string(const char* source_file,
						const char* source_func,
						const int source_line,
						int qid,
						const char* query,
						const int qoff,
						const int qend,
						const char* query_mapped_string,
						int tid,
						const u8* target,
						const int toff,
						const int tend,
						const char* target_mapped_string,
						const size_t align_size,
					    const BOOL right_extend)
{
	//return;
	int x = qoff, y = toff;
	for (size_t i = 0; i != align_size; ++i) {
		const char qc = query_mapped_string[i];
		if (qc != GAP_CHAR) {
            const char qc1 = right_extend ? query[x] : query[-x];
            if (qc != qc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, qc = %c, qc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d, align_size = %lu\n",
					  source_file,
					  source_func,
					  source_line,
					  qid,
					  tid,
					  right_extend,
					  i,
					  x,
					  y,
					  qc,
					  qc1,
					  qoff,
					  qend,
					  toff,
					  tend,
					  align_size);
                abort();
            }		  
			++x;
		}
		const char tc = target_mapped_string[i];
		if (tc != GAP_CHAR) {
            int c = right_extend ? target[y] : target[-y];
			const char tc1 = DECODE_RESIDUE(c);
            if (tc != tc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, tc = %c, tc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d\n",
						  source_func,
						  source_func,
						  source_line,
						  qid,
						  tid,
						  right_extend,
						  i,
						  x,
						  y,
						  tc,
						  tc1,
						  qoff,
						  qend,
						  toff,
						  tend);
                dump_align_string(query_mapped_string, target_mapped_string, align_size, stderr);
                abort();
            }
			++y;
		}
	}
}

char s_decode_bam_query_base(int c)
{
    // 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
    switch (c)
    {
    case 1:
        return 'A';
        break;
    case 2:
        return 'C';
        break;
    case 4:
        return 'G';
        break;
    case 8:
        return 'T';
        break;
    case 15:
        return 'N';
        break;
    default:
        fprintf(stderr, "ERROR: Illegal encoded base value %d\n", c);
        abort();
        break;
    }
}

char completement_residue(const char c)
{
    switch (c)
    {
    case 'A':
    case 'a':
        return 'T';
        break;
    case 'C':
    case 'c':
        return 'G';
        break;
    case 'G':
    case 'g':
        return 'C';
        break;
    case 'T':
    case 't':
        return 'A';
    case 'N':
    case 'n':
        return 'N';
    
    default:
        return 'A';
        break;
    }    
}

SAM_Parser::SAM_Parser()
{
    int p = 0;
    for (int i = 0; i < 64; ++i) M_codev1_table[p++] = i;
    for (int i = 64; i < 128; ++i) M_codev1_table[p++] = (i - 64) * 2 + 64;
    for (int i = 128; i < 192; ++i) M_codev1_table[p++] = (i - 128) * 4 + 192;
    for (int i = 192; i < 256; ++i) M_codev1_table[p++] = (i - 192) * 8 + 448;
    hbn_assert(p == kCodeV1TableSize);

    ks_initialize(&M_bam_aux);
}

static bool
s_extract_kinetic_values(bam1_t* bam, char tagname[], vector<SAM_Parser::kinetics_signal_value_type>& s_values)
{
    s_values.clear();
    uint8_t* tag = bam_aux_get(bam, tagname);
    if (!tag) {
        if (errno == ENOENT) {
            fprintf(stderr, "ERROR: BAM record is corruped\n");
            abort();
        }
        return false;
    }
    char type = tag[0];
    char subtype = tag[1];
    hbn_assert(type == 'B');
    hbn_assert(subtype == 'C' || subtype == 'S');
    const bool signal_is_encoded = (subtype == 'C');
    const int N = bam_auxB_len(tag);
    for (int i = 0; i < N; ++i) {
        int s = bam_auxB2i(tag, i);
        validate_signal_value(signal_is_encoded, s);
        if (!signal_is_encoded) s = s_encode_signal_value(s);
        s_values.push_back(s);
    }

    return true;
}

bool SAM_Parser::parse(sam_hdr_t* hdr, bam1_t* bam, const int read_id, const int min_query_size, HbnUnpackedDatabase* ref, SeqName2IdMap* ref_name2id)
{
    x_clear_kinetics_info();
    x_clear_map_info();

    char kFwdIpdTag[2] = { 'f', 'i' };
    char kRevIpdTag[2] = { 'r', 'i' };
    char kFwdPwTag[2] = { 'f', 'p' };
    char kRevPwTag[2] = { 'r', 'p' };
    bool check_kinetic_tags = bam_aux_get(bam, kFwdIpdTag) && bam_aux_get(bam, kRevIpdTag) && bam_aux_get(bam, kFwdPwTag) && bam_aux_get(bam, kRevPwTag);
    if (!check_kinetic_tags) return false;

    M_query_id = read_id;
    M_query_name = bam_get_qname(bam);
    M_query_strand = bam_is_rev(bam) ? REV : FWD;
    M_query_size = bam->core.l_qseq;
    if (M_query_size < min_query_size) return false;

    if (M_query_strand == FWD) {
        uint8_t* seq = bam_get_seq(bam);
        uint8_t* qual = bam_get_qual(bam);
        for (int i = 0; i < M_query_size; ++i) {
            int c = bam_seqi(seq, i);
            c = s_decode_bam_query_base(c);
            int rc = completement_residue(c);
            M_fwd_seq += c;
            M_rev_seq += rc;
            M_fwd_qual.push_back(qual[i] + 33);
            M_rev_qual.push_back(qual[i] + 33);
        }
        reverse(M_rev_seq.begin(), M_rev_seq.end());
        reverse(M_rev_qual.begin(), M_rev_qual.end());
    } else {
        uint8_t* seq = bam_get_seq(bam);
        uint8_t* qual = bam_get_qual(bam);
        for (int i = 0; i < M_query_size; ++i) {
            int rc = bam_seqi(seq, i);
            rc = s_decode_bam_query_base(rc);
            int c = completement_residue(rc);
            M_fwd_seq += c;
            M_rev_seq += rc;
            M_fwd_qual.push_back(qual[i] + 33);
            M_rev_qual.push_back(qual[i] + 33);
        }        
        reverse(M_fwd_seq.begin(), M_fwd_seq.end());
        reverse(M_fwd_qual.begin(), M_fwd_qual.end());
    }

    if (!s_extract_kinetic_values(bam, kFwdIpdTag, M_fwd_ipd)) return false;
    if (!s_extract_kinetic_values(bam, kRevIpdTag, M_rev_ipd)) return false;
    if (!s_extract_kinetic_values(bam, kFwdPwTag, M_fwd_pw)) return false;
    if (!s_extract_kinetic_values(bam, kRevPwTag, M_rev_pw)) return false;

    char kFwdPassTag[2] = { 'f', 'n' };
    char kRevPassTag[2] = { 'r', 'n' };
    uint8_t* aux = bam_aux_get(bam, kFwdPassTag);
    if (aux) M_fwd_pass = bam_aux2i(aux);
    aux = bam_aux_get(bam, kRevPassTag);
    if (aux) M_rev_pass = bam_aux2i(aux);

    M_effective_coverage = 0.0;
    char kECTag[2] = { 'e', 'c' };
    aux = bam_aux_get(bam, kECTag);
    if (aux) M_effective_coverage = bam_aux2f(aux);

    ////////////////// sanity check

    const bool dump_fail_info = false;
    if (M_fwd_seq.size() != M_fwd_qual.size()) {
        if (dump_fail_info) fprintf(stderr, "ERROR: at line %d (0-based, SAM/BAM header not counted):", read_id);
        if (dump_fail_info) fprintf(stderr, "[%s, %d] sequence length (%zu) and quality length (%zu) do not match\n", M_query_name.c_str(), M_query_id, M_fwd_seq.size(), M_fwd_qual.size());
        return false;
    }

    int n_ipd = M_fwd_ipd.size();
    int m_seq_size = M_fwd_seq.size();
    if (n_ipd != m_seq_size) {
        if (dump_fail_info) fprintf(stderr, "ERROR: at line %d (0-based, SAM/BAM header not counted):", read_id);
        if (dump_fail_info) fprintf(stderr, "[%s, %d] forward ipd values (%d) and sequence length (%d) do not match\n", M_query_name.c_str(), M_query_id, n_ipd, m_seq_size);
        return false;
    }

    n_ipd = M_rev_ipd.size();
    if (n_ipd != m_seq_size) {
        if (dump_fail_info) fprintf(stderr, "ERROR: at line %d (0-based, SAM/BAM header not counted):", read_id);
        if (dump_fail_info) fprintf(stderr, "[%s, %d] reverse ipd values (%d) and sequence length (%d) do not match\n", M_query_name.c_str(), M_query_id, n_ipd, m_seq_size);
        return false;
    }

    int n_pw = M_fwd_pw.size();
    if (n_pw != m_seq_size) {
        if (dump_fail_info) fprintf(stderr, "ERROR: at line %d (0-based, SAM/BAM header not counted):", read_id);
        if (dump_fail_info) fprintf(stderr, "[%s, %d] forward pw values (%d) and sequence length (%d) do not match\n", M_query_name.c_str(), M_query_id, n_pw, m_seq_size);
        return false;        
    }

    n_pw = M_rev_pw.size();
    if (n_pw != m_seq_size) {
        if (dump_fail_info) fprintf(stderr, "ERROR: at line %d (0-based, SAM/BAM header not counted):", read_id);
        if (dump_fail_info) fprintf(stderr, "[%s, %d] reverse pw values (%d) and sequence length (%d) do not match\n", M_query_name.c_str(), M_query_id, n_pw, m_seq_size);
        return false;        
    }

    //fprintf(stderr, "n_cigar = %d\n", bam->core.n_cigar);
    //// mapping info
    if (ref == nullptr || ref_name2id == nullptr || bam->core.n_cigar == 0) return true;
    M_is_mapped = true;
    uint32_t* cigar = bam_get_cigar(bam);
    M_cigar_op_list.clear();
    for (int i = 0; i < bam->core.n_cigar; ++i) {
        char op = bam_cigar_opchr(cigar[i]);
        int num = bam_cigar_oplen(cigar[i]);
        //fprintf(stderr, "%d\t%c\t%d\n", i, op, num);
        M_cigar_op_list.emplace_back(op, num);
    }
    //fprintf(stderr, "sid = %d\n", bam->core.tid);
    M_sname = hdr->target_name[bam->core.tid];
    M_sid = ref_name2id->GetIdFromNameSafe(M_sname);
    const u8* chr_seq = ref->GetSequence(M_sid, FWD);
    M_ss = ref->SeqSize(M_sid);
    const char* query_seq = (M_query_strand == FWD) ? fwd_seq() : rev_seq();
    cigar_to_alignment(query_seq, M_query_size, chr_seq + bam->core.pos, M_ss, M_cigar_op_list, M_qb, M_qe, M_sb, M_se, M_qas, M_sas, M_qas_pos_list, M_sas_pos_list);
    M_as_size = M_qas.size();
    M_sb += bam->core.pos;
    M_se += bam->core.pos;
    for (auto& p : M_sas_pos_list) p += bam->core.pos;
    ++M_qe;
    ++M_se;

    M_mapQ = bam->core.qual;
    M_identity = calc_effective_ident_perc(M_qas.c_str(), M_sas.c_str(), M_as_size);

    int qi = M_qb - 1, si = M_sb - 1;
    for (int i = 0; i < M_as_size; ++i) {
        if (M_qas[i] != GAP_CHAR) {
            int qc = M_qas[i];
            ++qi;
            hbn_assert(qi == M_qas_pos_list[i]);
            int qc1 = query_seq[qi];
            hbn_assert(qc == qc1);
        }
        if (M_sas[i] != GAP_CHAR) {
            int sc = M_sas[i];
            ++si;
            hbn_assert(si == M_sas_pos_list[i]);
            int sc1 = chr_seq[si];
            sc1 = DECODE_RESIDUE(sc1);
            hbn_assert(sc == sc1);
        }
    }

    return true;
}
