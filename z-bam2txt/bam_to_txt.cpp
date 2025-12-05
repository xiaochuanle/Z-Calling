#include "../corelib/sam_batch.hpp"
#include "../corelib/sam_parser.hpp"
#include "../corelib/seq_name2id_map.hpp"
#include "../corelib/unpacked_seqdb.hpp"
#include "../parameter/parameters.hpp"
#include "necat_info.hpp"
#include "sam_map_info.hpp"

#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct Options
{
    int min_read_size { MIN_READ_SIZE };
    int min_mapq { MIN_MAPQ };
    double min_identity { MIN_IDENTITY };
    //int min_strand_pass { MIN_STRAND_PASS };
    double min_ec { MIN_EC };
    bool skip_unmapped_motifs { false };
    int num_threads { NUM_THREADS };
    bool use_zp_prob_tag { false };

    const char* bam_path {nullptr} ;
    const char* output {nullptr};
    const char* ref_path { nullptr };

    void dump_usage(int argc, char* argv[]) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s [OPTIONS] BAM OUTPUT [Reference]\n", argv[0]);

        fprintf(stderr, "\n\n");
        fprintf(stderr, "OPTIONAL ARGUMENTS:\n");
        fprintf(stderr, "  -l <integer>\n");
        fprintf(stderr, "    Minumum read length\n");
        fprintf(stderr, "    Default = '%d'\n", MIN_READ_SIZE);
        fprintf(stderr, "  -q <integer>\n");
        fprintf(stderr, "    Minimum mapping quality\n");
        fprintf(stderr, "    Default = '%d'\n", MIN_MAPQ);
        fprintf(stderr, "  -i <real>\n");
        fprintf(stderr, "    Minimum mapping identity\n");
        fprintf(stderr, "    Default = '%g'\n", MIN_IDENTITY);
        //fprintf(stderr, "  -p <integer>\n");
        //fprintf(stderr, "    Skip reads with number of passes of forward/reversed strands less than INT\n");
        //fprintf(stderr, "    Default = '%d'\n", MIN_STRAND_PASS);
        fprintf(stderr, "  -ec <real>\n");
        fprintf(stderr, "    Minimum effective coverage\n");
        fprintf(stderr, "    Default = '%g'\n", MIN_EC);
        fprintf(stderr, "  -zp\n");
        fprintf(stderr, "    Use the ZP BAM tag to parse the modification probability\n");
        fprintf(stderr, "  -a\n");
        fprintf(stderr, "    Skip unmapped motifs in reads\n");
        fprintf(stderr, "  -t <integer>\n");
        fprintf(stderr, "    Number of CPU threads\n");
        fprintf(stderr, "    Default = '%d'\n", NUM_THREADS);
    }

    void dump() {
        fprintf(stderr, "\n");
        fprintf(stderr, "====> PARAMETERS:\n");
        fprintf(stderr, "Min-read-length: %d\n", min_read_size);
        fprintf(stderr, "Min-mapQ: %d\n", min_mapq);
        fprintf(stderr, "Min-identity: %g\n", min_identity);
        //fprintf(stderr, "Min-strand-pass: %d\n", min_strand_pass);
        fprintf(stderr, "Min-ec: %g\n", min_ec);
        fprintf(stderr, "Skip-unmapped-motifs: %d\n", skip_unmapped_motifs);
        fprintf(stderr, "CPU-threads: %d\n", num_threads);
        fprintf(stderr, "BAM: %s\n", bam_path);
        fprintf(stderr, "Output: %s\n", output);
        fprintf(stderr, "Reference: %s\n", ref_path);
    }

    bool parse(int argc, char* argv[])
    {
        int i = 1;
        while (i < argc) {
            bool r = (argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
            if (r) break;

            if (strcmp(argv[i], "-l") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                min_read_size = atoi(argv[i+1]);
                i += 2;
                continue;
            }

            if (strcmp(argv[i], "-q") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                min_mapq = atoi(argv[i+1]);
                i += 2;
                continue;
            }

            if (strcmp(argv[i], "-i") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                min_identity = atof(argv[i+1]);
                i += 2;
                continue;
            }

            if (strcmp(argv[i], "-ec") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                min_ec = atof(argv[i+1]);
                i += 2;
                continue;
            }

#if 0
            if (strcmp(argv[i], "-p") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                min_strand_pass = atoi(argv[i+1]);
                i += 2;
                continue;
            }
#endif

            if (strcmp(argv[i], "-t") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                num_threads = atoi(argv[i+1]);
                i += 2;
                continue;
            }

            if (strcmp(argv[i], "-a") == 0) {
                skip_unmapped_motifs = true;
                i += 1;
                continue;
            }

            if (strcmp(argv[i], "-zp") == 0) {
                use_zp_prob_tag = true;
                i += 1;
                continue;
            }

            fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
            return false;
        }

        if (i == argc) return false;
        bam_path = argv[i];
        ++i;

        if (i == argc) return false;
        output = argv[i];
        ++i;

        if (i == argc) return true;
        ref_path = argv[i];

        return true;
    }
};

class ThreadWorkData
{
public:
    ThreadWorkData(Options* options, int argc, char* argv[])
    {
        M_options = options;
        M_sam = new SAM_Batch(options->bam_path, 10000);//READ_BATCH_SIZE); 
        hbn_fopen(M_out, options->output, "w");

        M_ref = nullptr;
	M_ref_name2id = nullptr;
        if (options->ref_path) {
            M_ref = new HbnUnpackedDatabase(options->ref_path, true);
            M_ref->load_next_batch();
            M_ref_name2id = new SeqName2IdMap();
            int num_chr = M_ref->NumSeqs();
            for (int i = 0; i < num_chr; ++i) M_ref_name2id->add_one_name(M_ref->SeqName(i));
            M_ref_name2id->build_name2id_map();
        }

        M_batch_reads = 0;
        M_batch_bases = 0;
        M_batch_samples = 0;

	M_thread_id = 0;
    }

    ~ThreadWorkData()
    {
        delete M_sam;
        hbn_fclose(M_out);
        delete M_ref;
        delete M_ref_name2id;
    }

    size_t batch_samples()
    {
        return M_batch_samples;
    }

    int batch_reads()
    {
        return M_batch_reads;
    }

    size_t batch_bases()
    {
        return M_batch_bases;
    }

    int processed_reads()
    {
        return M_sam->loaded_reads();
    }

    bool get_next_sam(int& read_id, bam1_t* bam)
    {
        return M_sam->get_next_sam(read_id, bam);
    }

    bool reset_batch_read_idx()
    {
	    M_thread_id = 0;
        M_batch_reads = 0;
        M_batch_bases = 0;
        M_batch_samples = 0;
        return M_sam->reset_batch_read_idx();
    }

    Options* options()
    {
        return M_options;
    }

    sam_hdr_t* sam_hdr()
    {
        return M_sam->sam_hdr();
    }

    HbnUnpackedDatabase* reference()
    {
        return M_ref;
    }

    SeqName2IdMap* ref_name2id()
    {
        return M_ref_name2id;
    }

    void add_txt_calls(const char* s, size_t sl, int num_reads, size_t num_samples, size_t num_bases)
    {
        lock_guard<mutex> lg(M_out_mutex);
        hbn_fwrite(s, 1, sl, M_out);
        M_batch_samples += num_samples;
        M_batch_reads += num_reads;
        M_batch_bases += num_bases;
    }

    int thread_id()
    {
	    lock_guard<mutex> lg(M_out_mutex);
	    return M_thread_id++;
    }

private:
    Options*        M_options;
    SAM_Batch*      M_sam;
    int             M_thread_id;
    HbnUnpackedDatabase*    M_ref;
    SeqName2IdMap*          M_ref_name2id;

    size_t          M_batch_samples;
    int             M_batch_reads;
    size_t          M_batch_bases;

    FILE*           M_out;
    mutex           M_out_mutex;
};

struct BaseModInfo
{
    int qoff;
    int soff;
    double prob;
    bool is_valid;
    char base;
};

static bool
s_extract_mod_base_counts(bam1_t* bam, vector<int>& fwd_base_counts, vector<int>& rev_base_counts)
{
    fwd_base_counts.clear();
    rev_base_counts.clear();
    char mmtagname[2] = { 'M', 'M' };
    uint8_t* mmtag = bam_aux_get(bam, mmtagname);
    if (!mmtag) return false;

    const char* mms = bam_aux2Z(mmtag);
    const char* p = mms;
    while (p[0]) {
        if (p[0] == MOD_BASE_Z) {
            p += strlen("A+m");
            while (p[0] != ';') {
                ++p;
                int cnt = 0;
                while (isdigit(p[0])) {
                    cnt = cnt * 10 + p[0] - '0';
                    ++p;
                }
                fwd_base_counts.push_back(cnt);
            }
        } else if (p[0] == REV_MOD_BASE_Z) {
            p += strlen("T-m");
            while (p[0] != ';') {
                ++p;
                int cnt = 0;
                while (isdigit(p[0])) {
                    cnt = cnt * 10 + p[0] - '0';
                    ++p;
                }
                rev_base_counts.push_back(cnt);
            }
        } else {
            fprintf(stderr, "ERROR: Ill-formated MM BAM tag '%s'\n", mms);
            abort();
        }
        hbn_assert(p[0] == ';');
        ++p;
    }

    return true;
}

static bool
s_extract_mod_probs_from_bam(bam1_t* bam, const bool use_zp_prob_tag, vector<int>& fwd_base_counts, vector<int>& rev_base_counts, vector<double>& probs)
{
    if (!s_extract_mod_base_counts(bam, fwd_base_counts, rev_base_counts)) return false;

    probs.clear();
    if (use_zp_prob_tag) {
        char zptagname[2] = { 'Z', 'P' };
        uint8_t* zptag = bam_aux_get(bam, zptagname);
        if (!zptag) return false;
        int n = bam_auxB_len(zptag);
        for (int i = 0; i < n; ++i) {
            double prob = bam_auxB2f(zptag, i);
            probs.push_back(prob);
        }
    } else {
        char mltagname[2] = { 'M', 'L' };
        uint8_t* mltag = bam_aux_get(bam, mltagname);
        if (!mltagname) return false;
        int n = bam_auxB_len(mltag);
        for (int i = 0; i < n; ++i) {
            double prob = bam_auxB2i(mltag, i);
            prob /= 255.0;
            probs.push_back(prob);
        }
    }

    hbn_assert(fwd_base_counts.size() + rev_base_counts.size() == probs.size());
    return true;
}

static void*
work_thread(void* params)
{
    ThreadWorkData* data = (ThreadWorkData*)(params);
    HbnUnpackedDatabase* ref = data->reference();
    SeqName2IdMap* ref_name2id = data->ref_name2id();
    Options* options = data->options();
    int thread_id = data->thread_id();
    sam_hdr_t* sam_hdr = data->sam_hdr();
    int read_id;
    bam1_t* bam = bam_init1();
    SAM_MapInfo readinfo;
    vector<BaseModInfo> read_bases;
    vector<int> mm_fwd_base_counts, mm_rev_base_counts;
    vector<double> ml_probs;
    ostringstream os;
    string oss;

    while (data->get_next_sam(read_id, bam)) {
        if ((read_id % 1000) == 0) HBN_LOG("%d\t%10d reads processed", thread_id, read_id);
        if (!readinfo.parse(sam_hdr, bam, read_id, options->min_read_size, ref, ref_name2id)) continue;
        if (!s_extract_mod_probs_from_bam(bam, options->use_zp_prob_tag, mm_fwd_base_counts, mm_rev_base_counts, ml_probs)) continue;
        bool is_mapped = readinfo.is_mapped();
        if (is_mapped) is_mapped = (readinfo.mapQ() >= options->min_mapq) && (readinfo.identity() >= options->min_identity);
        if (options->skip_unmapped_motifs && is_mapped == false) continue;
        //if (readinfo.fwd_pass() < options->min_strand_pass || readinfo.rev_pass() < options->min_strand_pass) continue;
        if (readinfo.effective_coverage() < options->min_ec) continue;

        const int qd = readinfo.query_strand();
        const char* seq = (qd == FWD) ? readinfo.fwd_seq() : readinfo.rev_seq();
        const int seq_size = readinfo.query_size();

        read_bases.resize(seq_size);
        for (auto& rb : read_bases) {
            rb.qoff = -1;
            rb.soff = -1;
            rb.base = 'N';
            rb.is_valid = false;
            rb.prob = 0.0;
        }

        for (int i = 0; i <= seq_size - MOTIF_Z_LEN; ++i) {
            read_bases[i].qoff = -1;
            read_bases[i].soff = -1;
            read_bases[i].base = 'N';
            read_bases[i].is_valid = false;
            read_bases[i].prob = -1.0;

            if (seq[i] == MOD_BASE_Z) {
                read_bases[i].qoff = i;
                read_bases[i].base = MOD_BASE_Z;
            } else if (seq[i] == REV_MOD_BASE_Z) {
                read_bases[i].qoff = i;
                read_bases[i].base = REV_MOD_BASE_Z;
            }
        }

        if (is_mapped) {
            const char* qas = readinfo.qas();
            const char* sas = readinfo.sas();
            const int as_size = readinfo.as_size();
            const int* qas_pos_list = readinfo.qas_pos_list();
            const int* sas_pos_list = readinfo.sas_pos_list();
            for (int i = 0; i <= as_size - MOTIF_Z_LEN; ++i) {
                int qc = qas[i];
                int sc = sas[i];
                int qoff = qas_pos_list[i];
                int soff = sas_pos_list[i];

                if (qc == MOD_BASE_Z && sc == MOD_BASE_Z) {
                    hbn_assert(read_bases[qoff].qoff == qoff);
                    hbn_assert(read_bases[qoff].soff == -1);
                    read_bases[qoff].soff = soff;
                } else if (qc == REV_MOD_BASE_Z && sc == REV_MOD_BASE_Z) {
                    hbn_assert(read_bases[qoff].qoff == qoff);
                    hbn_assert(read_bases[qoff].soff == -1);
                    read_bases[qoff].soff = soff;
                }
            }
        }

        const char* fwd_seq = readinfo.fwd_seq();
        int qoff = 0;
        int ml_idx = 0;
        int fwd_counts = mm_fwd_base_counts.size();
        for (int i = 0; i < fwd_counts; ++i) {
            int cnt = 0;
            while (cnt < mm_fwd_base_counts[i]) {
                if (fwd_seq[qoff] == MOD_BASE_Z) ++cnt;
                ++qoff;
            }
            hbn_assert(cnt == mm_fwd_base_counts[i]);
            while (fwd_seq[qoff] != MOD_BASE_Z) ++qoff;
            hbn_assert(fwd_seq[qoff] == MOD_BASE_Z);
            if (qd == FWD) {
                read_bases[qoff].prob = ml_probs[ml_idx++];
                read_bases[qoff].is_valid = true;
            } else {
                int xqoff = reverse_strand_motif_pos(qoff, seq_size, MOTIF_Z_LEN);
                hbn_assert(read_bases[xqoff].base == REV_MOD_BASE_Z);
                read_bases[xqoff].prob = ml_probs[ml_idx++];
                read_bases[xqoff].is_valid = true;
            }
            ++qoff;
        }

        qoff = 0;
        int rev_counts = mm_rev_base_counts.size();
        for (int i = 0; i < rev_counts; ++i) {
            int cnt =0;
            while (cnt < mm_rev_base_counts[i]) {
                if (fwd_seq[qoff] == REV_MOD_BASE_Z) ++cnt;
                ++qoff;
            }
            hbn_assert(cnt == mm_rev_base_counts[i]);
            while (fwd_seq[qoff] != REV_MOD_BASE_Z) ++qoff;
            hbn_assert(fwd_seq[qoff] == REV_MOD_BASE_Z);
            if (qd == FWD) {
                read_bases[qoff].prob = ml_probs[ml_idx++];
                read_bases[qoff].is_valid = true;
            } else {
                int xqoff = reverse_strand_motif_pos(qoff, seq_size, MOTIF_Z_LEN);
                hbn_assert(read_bases[xqoff].base == MOD_BASE_Z);
                read_bases[xqoff].prob = ml_probs[ml_idx++];
                read_bases[xqoff].is_valid = true;
            }
            ++qoff;
        }
        hbn_assert(ml_idx == ml_probs.size());

        const char* qname = readinfo.query_name();
        const char* sname = readinfo.sname();
        const char strand = (qd == FWD) ? '+' : '-';
        os.str("");
        int num_samples = 0;
        for (int i = 0; i < seq_size; ++i) {
            if (read_bases[i].qoff == -1) continue;
            if (options->skip_unmapped_motifs && read_bases[i].soff == -1) continue;
            if (!read_bases[i].is_valid) continue;
	    //if (!num_samples) os << qname; else os << '*';
	        os << qname;
            os << '\t' << read_id << '\t' << strand << '\t' << i << '\t' << read_bases[i].base;
            if (read_bases[i].soff >= 0) os << '\t' << sname << '\t' << read_bases[i].soff;
            os << '\t' << read_bases[i].prob << '\n';
            ++num_samples;
        }
        oss = os.str();
        data->add_txt_calls(oss.c_str(), oss.size(), 1, num_samples, seq_size);
    }
  
    bam_destroy1(bam);
    return nullptr;
}

int main(int argc, char* argv[])
{
    Options options;
    if (!options.parse(argc, argv)) {
        options.dump_usage(argc, argv);
        return EXIT_FAILURE;
    }
    options.dump();
    HbnRunningInfo hbnrun;

    ThreadWorkData data(&options, argc, argv);
    const int num_threads = options.num_threads;
    pthread_t jobs[num_threads];

    int total_reads = 0;
    size_t total_samples = 0;
    size_t total_bases = 0;
    string base_size;
    while (data.reset_batch_read_idx()) {

        for (int i = 0; i < num_threads; ++i) {
            pthread_create(jobs + i, NULL, work_thread, &data);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(jobs[i], NULL);
        }

        int batch_reads = data.batch_reads();
        size_t batch_samples = data.batch_samples();
        size_t batch_bases = data.batch_bases();
        int processed_reads = data.processed_reads();

        base_size = NStr::UInt8ToString_DataSize(batch_bases);
        HBN_LOG("Extract %zu samples in %d reads (%s)", batch_samples, batch_reads, base_size.c_str());

        total_reads += batch_reads;
        total_samples += batch_samples;
        total_bases += batch_bases;
    }

    base_size = NStr::UInt8ToString_DataSize(total_bases);

    return 0;
}
