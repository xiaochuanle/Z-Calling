#ifndef __SAM_BATCH_HPP
#define __SAM_BATCH_HPP

#include "../3rdparty/htslib/include/sam.h"

#include <mutex>
#include <string>

class SAM_Batch
{
public:
    SAM_Batch(const char* bam_path, const int batch_size) 
    {
        M_batch_size = batch_size;
        M_batch_read_idx = 0;
        M_read_id = 0;

        M_sam = sam_open(bam_path, "rb");
	    hts_set_threads(M_sam, 8);
	    hts_set_cache_size(M_sam, 500000000);
        M_sam_hdr = sam_hdr_read(M_sam);
        M_sam_eof = false;

#if 0
        samFile* sam = sam_open(bam_path, "r");
        M_sam_hdr = sam_hdr_read(sam);
        sam_close(sam);
        M_bam_fp = hts_open(bam_path, "rb");
	    hts_set_threads(M_bam_fp, 4);
	    hts_set_cache_size(M_bam_fp, 100000000);
        M_bam_idx = sam_index_load(M_bam_fp, bam_path);
        M_bam_itr = sam_itr_queryi(M_bam_idx, HTS_IDX_START, 0, 0);
#endif
    }

    ~SAM_Batch()
    {
        sam_close(M_sam);
        sam_hdr_destroy(M_sam_hdr);

        //hts_itr_destroy(M_bam_itr);
        //hts_idx_destroy(M_bam_idx);
        //hts_close(M_bam_fp);
    }

    bool reset_batch_read_idx()
    {
        if (M_sam_eof) return false;
        M_batch_read_idx = 0;
        return true;
    }

    bool get_next_sam(int& read_id, bam1_t*& bam)
    {
        std::lock_guard<std::mutex> lg(M_batch_read_idx_mutex);

        if (M_sam_eof) return false;
        if (M_batch_read_idx >= M_batch_size) return false;
#if 1
        if (sam_read1(M_sam, M_sam_hdr, bam) < 0) {
            M_sam_eof = true;
            return false;
        }
#else
        if (sam_itr_next(M_bam_fp, M_bam_itr, bam) < 0) {
            M_sam_eof = true;
            return false;
        }
#endif

        read_id = M_read_id;
        ++M_read_id;
        ++M_batch_read_idx;
        return true;
    }

    int loaded_reads() {
        return M_read_id;
    }

    sam_hdr_t* sam_hdr() {
        return M_sam_hdr;
    }

private:
    int             M_batch_size;
    int             M_batch_read_idx;
    std::mutex      M_batch_read_idx_mutex;

    int             M_read_id;
    samFile*        M_sam;
    sam_hdr_t*      M_sam_hdr;
    bool            M_sam_eof;

    //htsFile*        M_bam_fp;
    //hts_idx_t*      M_bam_idx;
    //hts_itr_t*     M_bam_itr;
};

#endif // __SAM_BATCH_HPP
