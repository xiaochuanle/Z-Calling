//
// Created by admin on 2024/09/18.
//
#include <iostream>
#include "../3rdparty/htslib/include/sam.h"
#include <filesystem>
#include "../3rdparty/libsvm/svm.h"
#include "../3rdparty/argparse/argparse.h"
#include <vector>
#include <algorithm>
#include <random>

int train(size_t num_threshold,
          int feature_num,
          std::filesystem::path bam_file_pos,
          std::filesystem::path bam_file_neg,
          std::filesystem::path model_path,
          float C,
          float gamma);

int valid(size_t num_threshold,
          int feature_num,
          std::filesystem::path bam_file_pos,
          std::filesystem::path bam_file_neg,
          std::filesystem::path model_path);

int predict(size_t num_threshold,
            int feature_num,
            std::filesystem::path input_bam_file,
            std::filesystem::path output_file,
            std::filesystem::path model_path);
