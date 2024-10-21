#ifndef __UNPACKED_SEQDB_HPP
#define __UNPACKED_SEQDB_HPP

#include "fasta.hpp"
#include "hbn_aux.h"
#include "../str_util/ncbistr.hpp"

#include <string>
#include <vector>

class CUnpackedDatabase
{
public:
    struct CUnpackedSeqInfo
    {
        int seq_id;
        size_t seq_offset;
        int seq_size;
        size_t name_offset;
    };

public:
    CUnpackedDatabase(const char* fasta_path,
        bool load_fwd_seq = true,
        bool load_rev_seq = true,
        size_t batch_size = std::numeric_limits<size_t>::max()) {
        m_seqid_offset = 0;
        m_seqinfo.seq_id = 0;
        m_seqinfo.seq_offset = 0;
        m_seqinfo.seq_size = 0;
        m_seqinfo.name_offset = 0;
        m_load_fwd_seq = load_fwd_seq;
        m_load_rev_seq = load_rev_seq;
        m_batch_size = batch_size;
        m_input_list.push_back(fasta_path);
        m_input_idx = 0;
        m_in = NULL;
        open_next_stream();
    }

    CUnpackedDatabase(std::vector<std::string>& input_list,
        bool load_fwd_seq = true,
        bool load_rev_seq = true,
        size_t batch_size = std::numeric_limits<size_t>::max()) {
        m_seqid_offset = 0;
        m_seqinfo.seq_id = 0;
        m_seqinfo.seq_offset = 0;
        m_seqinfo.seq_size = 0;
        m_seqinfo.name_offset = 0;
        m_base_offset = 0;
        m_load_fwd_seq = load_fwd_seq;
        m_load_rev_seq = load_rev_seq;
        m_batch_size = batch_size;
        m_input_list = input_list;
        m_input_idx = 0;
        m_in = NULL;
        open_next_stream();
    }

    void clear() {
        m_seqid_offset = m_seqinfo.seq_id;
        m_base_offset = 0;
        m_name_list.clear();
        m_seqinfo_list.clear();
        m_fwd_seq_list.clear();
        m_rev_seq_list.clear();
        m_seqinfo.name_offset = 0;
        m_seqinfo.seq_offset = 0;
        m_seqinfo.seq_size = 0;
    }

    bool load_next_batch() {
        clear();
        size_t bases_loaded = 0;
        int seqs_loaded = 0;
        while (load_next_sequence()) {
            m_seqinfo.name_offset = m_name_list.size();
            m_seqinfo.seq_offset = m_base_offset;
            m_seqinfo.seq_size = m_in->sequence().size();
            m_seqinfo_list.push_back(m_seqinfo);
            ++m_seqinfo.seq_id;

            m_name_list.insert(m_name_list.end(), m_in->name().begin(), m_in->name().end());
            m_name_list.push_back('\0');

            bool update_base_offset = false;
            if (m_load_fwd_seq) {
                update_base_offset = true;
                for (auto c : m_in->sequence()) {
                    int ec = c;
                    ec = nst_nt4_table[ec];
                    if (ec > 3) ec = 0;
                    m_fwd_seq_list.push_back(ec);
                }
            }
            if (m_load_rev_seq) {
                update_base_offset = true;
                for (auto iter = m_in->sequence().rbegin(); iter != m_in->sequence().rend(); ++iter) {
                    int ec = *iter;
                    ec = nst_nt4_table[ec];
                    if (ec > 3) ec = 0;
                    ec = 3 - ec;
                    m_rev_seq_list.push_back(ec);
                }
            }
            if (update_base_offset) m_base_offset += m_in->sequence().size();

            ++seqs_loaded;
            bases_loaded += m_in->sequence().size();
            if (bases_loaded >= m_batch_size) break;
        }
        if (!seqs_loaded) return false;

        std::string loaded_size = NStr::UInt8ToString_DataSize(bases_loaded);
        HBN_LOG("Load %d sequences (%s)", seqs_loaded, loaded_size.c_str());
        return true;
    }

    int NumSeqs() const {
        return m_seqinfo_list.size();
    }

    const char* SeqName(const int seq_id) const {
        return m_name_list.data() + m_seqinfo_list[seq_id].name_offset;
    }

    const u8* GetSequence(const int seq_id, const int strand) const {
        const u8* s = (strand == FWD)
                      ? 
                      m_fwd_seq_list.data()
                      :
                      m_rev_seq_list.data();
        if (strand == REV) hbn_assert(!m_rev_seq_list.empty());
        return s + m_seqinfo_list[seq_id].seq_offset;
    }

    int SeqSize(const int seq_id) const {
        return m_seqinfo_list[seq_id].seq_size;
    }

    int SeqIdOffset() const {
        return m_seqid_offset;
    }

private:
    bool load_next_sequence() {
        if (!m_in) return false;
        while (m_in->AtEOF()) {
            if (!open_next_stream()) return false;
        }
        if (m_in == NULL || m_in->AtEOF()) return false;
        m_in->ReadOneSeq();
        return true;
    }

    bool open_next_stream() {
        if (m_in) {
            delete m_in;
            m_in = NULL;
        }
        if (m_input_idx >= m_input_list.size()) return false;
        m_in = new HbnFastaReader(m_input_list[m_input_idx++].c_str());
        return true;
    }

private:
    int                 m_seqid_offset;
    CUnpackedSeqInfo    m_seqinfo;
    std::vector<char>   m_name_list;
    std::vector<CUnpackedSeqInfo> m_seqinfo_list;
    size_t          m_base_offset;
    std::vector<u8> m_fwd_seq_list;
    std::vector<u8> m_rev_seq_list;
    bool            m_load_fwd_seq;
    bool            m_load_rev_seq;
    size_t          m_batch_size;
    
    std::vector<std::string> m_input_list;
    size_t                   m_input_idx;
    HbnFastaReader*          m_in;
};

typedef CUnpackedDatabase HbnUnpackedDatabase;

#endif // __UNPACKED_SEQDB_HPP