#ifndef __BUILD_MOD_BAM_HPP
#define __BUILD_MOD_BAM_HPP

#include "../corelib/sam_batch.hpp"
#include "../corelib/sam_parser.hpp"
#include "make_read_features.hpp"

#include <set>

void build_one_mod_bam(bam1_t* bam, const std::set<int>& skipped_tags,
    MolModCallz* fwd_ma, const int fwd_mc,
    MolModCallz* rev_ma, const int rev_mc);

#endif // __BUILD_MOD_BAM_HPP