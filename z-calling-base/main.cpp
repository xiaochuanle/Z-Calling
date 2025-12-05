#include "../corelib/sam_batch.hpp"
#include "../corelib/sam_parser.hpp"
#include "../parameter/parameters.hpp"
#include "build_mod_bam.hpp"
#include "make_read_features.hpp"
#include "necat_info.hpp"

#include <ATen/Parallel.h>
#include <torch/csrc/api/include/torch/nn/functional/activation.h>
#include <torch/script.h>

#include <iostream>
#include <mutex>
#include <set>
#include <vector>

using namespace std;

struct Options
{
    int min_read_size { MIN_READ_SIZE };
    int kmer_size { KMER_SIZE };
    int sample_batch_size { 1024 };
    bool keep_hifi_kinetics { false };
    int num_threads { NUM_THREADS };

    const char* bam_path {nullptr};
    const char* model_path {nullptr};
    const char* mod_bam_path {nullptr};

    void dump_usage(int argc, char* argv[]) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s [OPTIONS] BAM MODEL MOD\n", argv[0]);

        fprintf(stderr, "\n\n");
        fprintf(stderr, "OPTIONAL ARGUMENTS:\n");
        fprintf(stderr, "  -l <integer>\n");
        fprintf(stderr, "    Minumum read length\n");
        fprintf(stderr, "    Default = '%d'\n", MIN_READ_SIZE);
        fprintf(stderr, "  -k <integer>\n");
        fprintf(stderr, "    Kmer\n");
        fprintf(stderr, "    Default = '%d'\n", kmer_size);
        fprintf(stderr, "  -t <integer>\n");
        fprintf(stderr, "    Number of CPU threads\n");
        fprintf(stderr, "    Default = '%d'\n", NUM_THREADS);
        fprintf(stderr, "  -b <integer>\n");
        fprintf(stderr, "    Number of samples in one batch\n");
        fprintf(stderr, "    Default = '%d'\n", EVAL_SAMPLE_BATCH_SIZE);
        fprintf(stderr, "  -e\n");
        fprintf(stderr, "    Save Hifi kinecitcs in mod BAM\n");
    }

    void dump() {
        fprintf(stderr, "\n");
        fprintf(stderr, "====> PARAMETERS:\n");
        fprintf(stderr, "Min-read-length: %d\n", min_read_size);
        fprintf(stderr, "KMER: %d\n", kmer_size);
        fprintf(stderr, "Save-hifi-kinetics: %d\n", keep_hifi_kinetics);
        fprintf(stderr, "Batch-size: %d\n", sample_batch_size);
        fprintf(stderr, "CPU-threads: %d\n", num_threads);
        fprintf(stderr, "BAM: %s\n", bam_path);
        fprintf(stderr, "Model: %s\n", model_path);
        fprintf(stderr, "Output: %s\n", mod_bam_path);
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

            if (strcmp(argv[i], "-k") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                kmer_size = atoi(argv[i+1]);
                i += 2;
                continue;
            }

            if (strcmp(argv[i], "-t") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                num_threads = atoi(argv[i+1]);
                i += 2;
                continue;
            }

            if (strcmp(argv[i], "-e") == 0) {
                keep_hifi_kinetics = false;
                i += 1;
                continue;
            }

            if (strcmp(argv[i], "-b") == 0) {
                if (i + 1 == argc) {
                    fprintf(stderr, "ERROR: argument to mandatory option '%s' is missing\n", argv[i]);
                    return false;
                }
                sample_batch_size = atoi(argv[i+1]);
                i += 2;
                continue;
            }

            fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
            return false;
        }

        if (i == argc) return false;
        bam_path = argv[i];
        ++i;

        if (i == argc) return false;
        model_path = argv[i];
        ++i;

        if (i == argc) return false;
        mod_bam_path = argv[i];

        return true;
    }
};

static void
add_cmd_to_sam_hdr(int argc, char* argv[], sam_hdr_t* hdr)
{
    ostringstream os;
    os << "@PG"
       << "\tID:DNA-mods"
       << "\tPN:DNA-mods"
       << "\tPP:DNA-mods"
       << "\tVN:1.0.0"
       << "\tCL:";

    os << argv[0];
    for (int i = 1; i < argc; ++i) os << ' ' << argv[i];
    os << '\n';

    string cmd = os.str();
    sam_hdr_add_lines(hdr, cmd.c_str(), cmd.size());
}

class ThreadWorkData
{
public:
    ThreadWorkData(Options* options, int argc, char* argv[])
    {
        M_options = options;
        M_sam = new SAM_Batch(options->bam_path, 100000);//READ_BATCH_SIZE);

        sam_hdr_t* hdr_in = M_sam->sam_hdr();
        sam_hdr_t* hdr_out = sam_hdr_dup(hdr_in);
        M_out = sam_open(options->mod_bam_path, "wb");
	hts_set_threads(M_out, 8);
	hts_set_cache_size(M_out, 100000000);
        add_cmd_to_sam_hdr(argc, argv, hdr_out);
        sam_hdr_write(M_out, hdr_out);
        sam_hdr_destroy(hdr_out);

        M_batch_reads = 0;
        M_batch_bases = 0;
        M_batch_samples = 0;
    }

    ~ThreadWorkData()
    {
        delete M_sam;
        sam_close(M_out);
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

    void add_mod_bams(vector<pair<int, bam1_t*>>& mod_bam_list, int num_reads, size_t num_samples, size_t num_bases)
    {
        lock_guard<mutex> lg(M_mod_list_mutex);
        M_mod_bam_list.insert(M_mod_bam_list.end(), mod_bam_list.begin(), mod_bam_list.end());
        M_batch_reads += num_reads;
        M_batch_samples += num_samples;
        M_batch_bases += num_bases;
    }

    void save_mod_bam()
    {
        HBN_LOG("Save %zu mod bams", M_mod_bam_list.size());
        sort(M_mod_bam_list.begin(), M_mod_bam_list.end(), [](const pair<int, bam1_t*>& x, const pair<int, bam1_t*>& y) { return x.first < y.first; });
        sam_hdr_t* hdr = M_sam->sam_hdr();
        for (auto& rdmod : M_mod_bam_list) {
            sam_write1(M_out, hdr, rdmod.second);
            bam_destroy1(rdmod.second);
        }
        string size = NStr::UInt8ToString_DataSize(M_batch_bases);
        HBN_LOG("Extract %zu samples in %d reads (%s)", M_batch_samples, M_batch_reads, size.c_str());
        M_batch_reads = 0;
        M_batch_bases = 0;
        M_batch_samples = 0;
        M_mod_bam_list.clear();
    }

private:
    Options*        M_options;
    SAM_Batch*      M_sam;

    size_t          M_batch_samples;
    int             M_batch_reads;
    size_t          M_batch_bases;

    samFile*                        M_out;
    vector<pair<int, bam1_t*>>      M_mod_bam_list;
    mutex                           M_mod_list_mutex;
};

static void
fill_skipped_tags(set<int>& skipped_tags)
{
    int t0, t1, tx;

    t0 = 'f';
    t1 = 'p';
    tx = (t0 << 8) | t1;
    skipped_tags.insert(tx);

    t0 = 'r';
    t1 = 'p';
    tx = (t0 << 8) | t1;
    skipped_tags.insert(tx);

    t0 = 'f';
    t1 = 'i';
    tx = (t0 << 8) | t1;
    skipped_tags.insert(tx);

    t0 = 'r';
    t1 = 'i';
    tx = (t0 << 8) | t1;
    skipped_tags.insert(tx);
}

static void*
work_thread(void* params)
{
    ThreadWorkData* data = (ThreadWorkData*)(params);
    Options* options = data->options();
    sam_hdr_t* sam_hdr = data->sam_hdr();
    hbn_assert(options->kmer_size & 1);
    const int flanking_bases = options->kmer_size / 2;
    namespace F = torch::nn::functional;
    SAM_Parser readinfo;
    int read_id;
    bam1_t* bam = bam_init1();
    double phred33_to_prob_table[256];
    build_phred33_to_prob_table(phred33_to_prob_table);
    auto tensor_options = torch::TensorOptions().device(torch::kCPU).dtype(torch::kFloat32).requires_grad(false);
    vector<float> batch_features;
    set<int> skipped_bam_tags; if (!options->keep_hifi_kinetics) fill_skipped_tags(skipped_bam_tags);

    vector<pair<int, bam1_t*>> mod_bam_list;
    vector<MolModCallz> all_mod_calls, batch_mod_calls;
    int num_batch_samples = 0;
    int num_reads = 0;
    size_t num_bases = 0;
    std::vector<double> fwd_qual_probs;
    std::vector<double> rev_qual_probs;
    std::vector<double> fwd_ipd_features;
    std::vector<double> rev_ipd_features;
    std::vector<double> fwd_pw_features;
    std::vector<double> rev_pw_features;
    vector<double> kmer_features;

    torch::jit::script::Module model = torch::jit::load(options->model_path);
    model.eval();

    while (data->get_next_sam(read_id, bam)) {
        if ((read_id % 1000) == 0) HBN_LOG("%10d reads processed", read_id);
        mod_bam_list.emplace_back(read_id, bam);
        if (!readinfo.parse(sam_hdr, bam, read_id, options->min_read_size)) {
            bam = bam_init1();
            continue;
        }

        build_read_base_features(readinfo, 
            phred33_to_prob_table,
            fwd_qual_probs,
            rev_qual_probs,
            fwd_ipd_features,
            rev_ipd_features,
            fwd_pw_features,
            rev_pw_features);

        const char* fwd_seq = readinfo.fwd_seq();
        const int seq_size = readinfo.query_size();
        //HBN_LOG("process %d:%s:%d", read_id, readinfo.query_name(), seq_size);
        for (int i = 0; i <= seq_size - MOTIF_Z_LEN; ++i) {
            if (!motif_match_z(fwd_seq, i)) continue;
            MolModCallz mod;
            mod.id = read_id;
            mod.strand = FWD;
            mod.offset = i;
            mod.prob = 0.0;
            batch_mod_calls.push_back(mod);

            int strand = 1 - mod.strand;
            int offset = seq_size - 1 - mod.offset;
            build_kmer_features(readinfo.fwd_seq(),
                readinfo.rev_seq(),
                fwd_qual_probs,
                rev_qual_probs,
                fwd_ipd_features,
                rev_ipd_features,
                fwd_pw_features,
                rev_pw_features,
                readinfo.fwd_pass(),
                readinfo.rev_pass(),
                strand,
                offset,
                seq_size,
                flanking_bases,
                batch_features);

            ++num_batch_samples;
            if (num_batch_samples == options->sample_batch_size) {
                auto batch_of_tensors = torch::from_blob(batch_features.data(), {num_batch_samples, options->kmer_size * EVAL_KMER_BASE_FEATURES}, tensor_options);
	        auto out = model.forward({batch_of_tensors}).toTensor();
                at::Tensor prob = F::softmax(out, F::SoftmaxFuncOptions(1));
                auto prob_a = prob.accessor<float, 2>();
                for (int k = 0; k < num_batch_samples; ++k) batch_mod_calls[k].prob = prob_a[k][1];
                all_mod_calls.insert(all_mod_calls.end(), batch_mod_calls.begin(), batch_mod_calls.end());

                batch_features.clear();
                batch_mod_calls.clear();
                num_batch_samples = 0;
            }
        }

        const char* rev_seq = readinfo.rev_seq();
        for (int i = 0; i <= seq_size - MOTIF_Z_LEN; ++i) {
            if (!motif_match_z(rev_seq, i)) continue;
            MolModCallz mod;
            mod.id = read_id;
            mod.strand = REV;
            mod.offset = i;
            mod.prob = 0.0;
            batch_mod_calls.push_back(mod);

            int strand = 1 - mod.strand;
            int offset = seq_size - 1 - mod.offset;
            build_kmer_features(readinfo.fwd_seq(),
                readinfo.rev_seq(),
                fwd_qual_probs,
                rev_qual_probs,
                fwd_ipd_features,
                rev_ipd_features,
                fwd_pw_features,
                rev_pw_features,
                readinfo.fwd_pass(),
                readinfo.rev_pass(),
                strand,
                offset,
                seq_size,
                flanking_bases,
                batch_features);

            ++num_batch_samples;
            if (num_batch_samples == options->sample_batch_size) {
                auto batch_of_tensors = torch::from_blob(batch_features.data(), {num_batch_samples, options->kmer_size * EVAL_KMER_BASE_FEATURES}, tensor_options);
	            auto out = model.forward({batch_of_tensors}).toTensor();
                at::Tensor prob = F::softmax(out, F::SoftmaxFuncOptions(1));
                auto prob_a = prob.accessor<float, 2>();
                for (int k = 0; k < num_batch_samples; ++k) batch_mod_calls[k].prob = prob_a[k][1];
                all_mod_calls.insert(all_mod_calls.end(), batch_mod_calls.begin(), batch_mod_calls.end());

                batch_features.clear();
                batch_mod_calls.clear();
                num_batch_samples = 0;
            }
        }

        ++num_reads;
        num_bases += seq_size;
        bam = bam_init1();
    }

    if (num_batch_samples) {
        auto batch_of_tensors = torch::from_blob(batch_features.data(), {num_batch_samples, options->kmer_size * EVAL_KMER_BASE_FEATURES}, tensor_options);
	    auto out = model.forward({batch_of_tensors}).toTensor();
        at::Tensor prob = F::softmax(out, F::SoftmaxFuncOptions(1));
        auto prob_a = prob.accessor<float, 2>();
        for (int k = 0; k < num_batch_samples; ++k) batch_mod_calls[k].prob = prob_a[k][1];
        all_mod_calls.insert(all_mod_calls.end(), batch_mod_calls.begin(), batch_mod_calls.end());

        batch_features.clear();
        batch_mod_calls.clear();
        num_batch_samples = 0;
    }

    map<int, pair<MolModCallz*, int>> read_mod_maps;
    MolModCallz* ma = all_mod_calls.data();
    size_t mc = all_mod_calls.size();
    size_t i = 0;
    while (i < mc) {
        size_t j = i + 1;
        while (j < mc && ma[i].id == ma[j].id) ++j;
        read_mod_maps[ma[i].id] = pair<MolModCallz*, int>(ma + i, j - i);
        i = j;
    }

    for (auto& rdmod : mod_bam_list) {
        MolModCallz* fwd_ma = nullptr;
        int fwd_mc = 0;
        MolModCallz* rev_ma = nullptr;
        int rev_mc = 0;
        auto pos = read_mod_maps.find(rdmod.first);
        if (pos != read_mod_maps.end()) {
            ma = pos->second.first;
            mc = pos->second.second;
            fwd_ma = ma;
            for (int k = 0; k < mc; ++k) {
                if (ma[k].strand != FWD) break;
                ++fwd_mc;
            }
            rev_ma = ma + fwd_mc;
            rev_mc = mc - fwd_mc;
            int seq_size = rdmod.second->core.l_qseq;
            for (int k = 0; k < rev_mc; ++k) {
                rev_ma[k].offset = seq_size - 1 - rev_ma[k].offset;
            }
	    if (rev_mc) reverse(rev_ma, rev_ma + rev_mc);
        }
        build_one_mod_bam(rdmod.second, skipped_bam_tags, fwd_ma, fwd_mc, rev_ma, rev_mc);
    }
    data->add_mod_bams(mod_bam_list, num_reads, all_mod_calls.size(), num_bases);    

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

    at::set_num_threads(1);
    at::set_num_interop_threads(1);

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

        data.save_mod_bam();
        HBN_LOG("%10d reads processed", processed_reads);

        total_reads += batch_reads;
        total_samples += batch_samples;
        total_bases += batch_bases;
    }

    base_size = NStr::UInt8ToString_DataSize(total_bases);

    return 0;
}
