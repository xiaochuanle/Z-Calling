#include "../parameter/parameters.hpp"
#include "../3rdparty/htslib/include/hts_endian.h"
#include "build_mod_bam.hpp"

#include <sstream>
#include <string>

using namespace std;

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end)
{
    int size;
    uint32_t n;
    if (s >= end) return end;
    size = aux_type2size(*s); ++s; // skip type
    switch (size) {
    case 'Z':
    case 'H':
        while (s < end && *s) ++s;
        return s < end ? s + 1 : end;
    case 'B':
        if (end - s < 5) return NULL;
        size = aux_type2size(*s); ++s;
        n = le_to_u32(s);
        s += 4;
        if (size == 0 || end - s < size * n) return NULL;
        return s + size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size) return NULL;
        return s + size;
    }
}

static void
iterate_bam_tags(bam1_t* bam, const set<int>& skipped_tags)
{
    uint8_t* end = bam->data + bam->l_data;
    uint8_t* s_from = bam_get_aux(bam);
    uint8_t* s_to = s_from;
    kstring_t tag = KS_INITIALIZE;
    char tagname[2];

    while (s_from < end) {
        int x = (int)s_from[0]<<8 | s_from[1];
        int t0 = s_from[0];
        int t1 = s_from[1];
        int tx = (t0 << 8) | t1;
        if (skipped_tags.find(tx) == skipped_tags.end()) {
            tagname[0] = t0;
            tagname[1] = t1;
            ks_clear(&tag);
            bam_aux_get_str(bam, tagname, &tag);
            fprintf(stderr, "%s\n", ks_c_str(&tag));
        }
        uint8_t* s = skip_aux(s_from + 2, end);
        s_from = s;
    }

    ks_free(&tag);
}

static void
s_remove_skipped_tags(bam1_t* bam, const std::set<int>& skipped_tags)
{
    char tagname[2];
    for (auto s : skipped_tags) {
        int t1 = s & 255;
        int t0 = s >> 8;
        tagname[0] = t0;
        tagname[1] = t1;
        uint8_t* aux = bam_aux_get(bam, tagname);
        if (aux) bam_aux_del(bam, aux);
    }

    tagname[0] = 'M';
    tagname[1] = 'L';
    uint8_t* aux = bam_aux_get(bam, tagname);
    if (aux) bam_aux_del(bam, aux);

    tagname[0] = 'M';
    tagname[1] = 'M';
    aux = bam_aux_get(bam, tagname);
    if (aux) bam_aux_del(bam, aux);

    tagname[0] = 'Z';
    tagname[1] = 'P';
    aux = bam_aux_get(bam, tagname);
    if (aux) bam_aux_del(bam, aux);
}

static void
s_add_one_mm_cnt(int cnt, vector<uint8_t>& mmtags)
{
    char buf[64];
    int p = 0;
    do {
        buf[p] = (cnt % 10) + '0';
        ++p;
        cnt /= 10;
    } while (cnt);
    reverse(buf, buf + p);
    for (int i = 0; i < p; ++i) mmtags.push_back(buf[i]);
}

extern
char s_decode_bam_query_base(int c);

extern
char completement_residue(const char c);

static inline char
s_get_bam_fwd_base(uint8_t* bam_seq, int l_seq, int strand, int offset)
{
        offset = (strand == FWD) ? offset : (l_seq - 1 - offset);
        char c = s_decode_bam_query_base(bam_seqi(bam_seq, offset));
        return (strand == FWD) ? c : completement_residue(c);
}

void build_one_mod_bam(bam1_t* bam, const std::set<int>& skipped_tags,
    MolModCallz* fwd_ma, const int fwd_mc,
    MolModCallz* rev_ma, const int rev_mc)
{
    s_remove_skipped_tags(bam, skipped_tags);
    if (fwd_mc == 0 && rev_mc == 0) return;

    uint8_t* bam_seq = bam_get_seq(bam);
    int strand = bam_is_rev(bam);
    vector<uint8_t> mmtags;

    if (fwd_mc) {
        mmtags.push_back(MOD_BASE_Z);
        mmtags.push_back('+');
        mmtags.push_back('m');
        int last_qoff = 0;
        for (int i = 0; i < fwd_mc; ++i) {
	    if (i < fwd_mc - 1) hbn_assert(fwd_ma[i].offset < fwd_ma[i+1].offset);
            int c = s_get_bam_fwd_base(bam_seq, bam->core.l_qseq, strand, fwd_ma[i].offset);
            hbn_assert(c == MOD_BASE_Z);
            int cnt = 0;
            for (int k = last_qoff; k < fwd_ma[i].offset; ++k) {
                int c = s_get_bam_fwd_base(bam_seq, bam->core.l_qseq, strand, k);
                if (c == MOD_BASE_Z) ++cnt;
            }
            mmtags.push_back(',');
            s_add_one_mm_cnt(cnt, mmtags);
            last_qoff = fwd_ma[i].offset + 1;
        }
        mmtags.push_back(';');
    }

    if (rev_mc) {
        mmtags.push_back(REV_MOD_BASE_Z);
        mmtags.push_back('-');
        mmtags.push_back('m');
        int last_qoff = 0;
        for (int i = 0; i < rev_mc; ++i) {
	    if (i < rev_mc - 1) hbn_assert(rev_ma[i].offset < rev_ma[i+1].offset);
            int c = s_get_bam_fwd_base(bam_seq, bam->core.l_qseq, strand, rev_ma[i].offset);
            hbn_assert(c == REV_MOD_BASE_Z);
            int cnt = 0;
            for (int k = last_qoff; k < rev_ma[i].offset; ++k) {
                int c = s_get_bam_fwd_base(bam_seq, bam->core.l_qseq, strand, k);
                if (c == REV_MOD_BASE_Z) ++cnt;
            }
            mmtags.push_back(',');
            s_add_one_mm_cnt(cnt, mmtags);
            last_qoff = rev_ma[i].offset + 1;
        }
        mmtags.push_back(';');
    }

    mmtags.push_back('\0');
    char tagname[2];
    tagname[0] = 'M';
    tagname[1] = 'M';
    bam_aux_append(bam, tagname, 'Z', mmtags.size(), mmtags.data());

    vector<uint8_t>& mltags = mmtags;
    mltags.clear();
    mltags.push_back('C');
    vector<float> zptags;
    uint8_t buf[4];
    u32_to_le(fwd_mc + rev_mc, buf);
    mltags.insert(mltags.end(), buf, buf + 4);
    for (int i = 0; i < fwd_mc; ++i) {
        uint32_t n = fwd_ma[i].prob * 256;
        uint8_t c = min<int>(255, n);
        mltags.push_back(c);
        zptags.push_back(fwd_ma[i].prob);
    }
    for (int i = 0; i < rev_mc; ++i) {
        uint32_t n = rev_ma[i].prob * 256;
        uint8_t c = min<int>(255, n);
        mltags.push_back(c);
        zptags.push_back(rev_ma[i].prob);
    }
    tagname[0] = 'M';
    tagname[1] = 'L';
    bam_aux_append(bam, tagname, 'B', mltags.size(), mltags.data());

    tagname[0] = 'Z';
    tagname[1] = 'P';
    bam_aux_update_array(bam, tagname, 'f', zptags.size(), zptags.data());
}