//
// Created by admin on 2024/09/18.
//

#include <fstream>
#include "func.h"
#include "BamRead.h"

int train(size_t num_threshold,
          int feature_num,
          std::filesystem::path bam_file_pos,
          std::filesystem::path bam_file_neg,
          std::filesystem::path model_path,
          float C,
          float gamma) {
    std::cout << "Loading positive reads bam" << std::endl;
    samFile *bam_in = sam_open(bam_file_pos.c_str(), "r");
    bam_hdr_t *bam_header = sam_hdr_read(bam_in);
    bam1_t *aln = bam_init1();
    std::vector<std::shared_ptr<BamRead>> read_vec;
    while (sam_read1(bam_in, bam_header, aln) >= 0) {
        std::shared_ptr<BamRead> bam_ptr = std::make_shared<BamRead>(aln, 1, feature_num);
        if (bam_ptr->p_num < num_threshold || !(bam_ptr->is_valid)) {
            continue;
        }
        read_vec.push_back(bam_ptr);
    }

    std::cout << "Loading negative reads bam" << std::endl;
    bam_in = sam_open(bam_file_neg.c_str(), "r");
    bam_header = sam_hdr_read(bam_in);
    aln = bam_init1();
    while (sam_read1(bam_in, bam_header, aln) >= 0) {
        std::shared_ptr<BamRead> bam_ptr = std::make_shared<BamRead>(aln, -1, feature_num);
        if (bam_ptr->p_num < num_threshold || !(bam_ptr->is_valid)) {
            continue;
        }
        read_vec.push_back(bam_ptr);
    }

    std::cout << "Shuffling..." << std::endl;
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(read_vec.begin(), read_vec.end(), g);

    std::cout << "Setting feature..." << std::endl;
    svm_problem prob;
    prob.l = read_vec.size();
    prob.y = new double[read_vec.size()];
    prob.x = new svm_node *[read_vec.size()];

    for (int i = 0; i < read_vec.size(); ++i) {
        prob.x[i] = new svm_node[feature_num + 1];
        std::vector<double> feature = read_vec[i]->probability_proportions;
        for (int j = 0; j < feature_num; ++j) {
            prob.x[i][j].index = j + 1;
            prob.x[i][j].value = feature[j];
        }
        prob.x[i][feature_num].index = -1;
        prob.y[i] = read_vec[i]->label;
    }


    // 设置 SVM 参数
    svm_parameter param;
    param.svm_type = C_SVC;
    param.kernel_type = RBF;
    param.gamma = gamma;
    param.C = C;
    param.cache_size = 10000;
    param.eps = 1e-3;

    // 检查参数是否正确
    const char *err = svm_check_parameter(&prob, &param);
    if (err) {
        std::cerr << "Parameter error: " << err << std::endl;
        return -1;
    }

    std::cout << "Training..." << std::endl;
    // 训练 SVM 模型
    svm_model *model = svm_train(&prob, &param);
    svm_save_model(model_path.c_str(), model);

    std::cout << "Valid..." << std::endl;
    int32_t accurate_num = 0;
    double epsilon = 1e-6;
    for (int i = 0; i < 50000 && i < read_vec.size(); ++i) {
        svm_node test_node[feature_num + 1];
        std::vector<double> feature = read_vec[i]->probability_proportions;
        for (int j = 0; j < feature_num; ++j) {
            test_node[j].index = j + 1;
            test_node[j].value = feature[j];
        }
        test_node[feature_num].index = -1;
        double result = svm_predict(model, test_node);
        if (abs(result - read_vec[i]->label) < epsilon) {
            accurate_num++;
        }
    }
    float accuracy = (float) accurate_num / 50000;
    std::cout << "Prediction accuracy: " << accuracy << std::endl;

    svm_free_and_destroy_model(&model);
    delete[] prob.y;
    for (int i = 0; i < read_vec.size(); ++i) {
        delete[] prob.x[i];
    }
    delete[] prob.x;
    return 0;
}

int valid(size_t num_threshold,
          int feature_num,
          std::filesystem::path bam_file_pos,
          std::filesystem::path bam_file_neg,
          std::filesystem::path model_path) {
    std::cout << "Loading positive reads bam" << std::endl;
    samFile *bam_in = sam_open(bam_file_pos.c_str(), "r");
    bam_hdr_t *bam_header = sam_hdr_read(bam_in);
    bam1_t *aln = bam_init1();
    std::vector<std::shared_ptr<BamRead>> read_vec;
    while (sam_read1(bam_in, bam_header, aln) >= 0) {
        std::shared_ptr<BamRead> bam_ptr = std::make_shared<BamRead>(aln, 1, feature_num);
        if (bam_ptr->p_num < num_threshold || !(bam_ptr->is_valid)) {
            continue;
        }
        read_vec.push_back(bam_ptr);
    }

    std::cout << "Loading negative reads bam" << std::endl;
    bam_in = sam_open(bam_file_neg.c_str(), "r");
    bam_header = sam_hdr_read(bam_in);
    aln = bam_init1();
    while (sam_read1(bam_in, bam_header, aln) >= 0) {
        std::shared_ptr<BamRead> bam_ptr = std::make_shared<BamRead>(aln, -1, feature_num);
        if (bam_ptr->p_num < num_threshold || !(bam_ptr->is_valid)) {
            continue;
        }
        read_vec.push_back(bam_ptr);
    }

    std::cout << "Setting feature..." << std::endl;
    svm_problem prob;
    prob.l = read_vec.size();
    prob.y = new double[read_vec.size()];
    prob.x = new svm_node *[read_vec.size()];
    for (int i = 0; i < read_vec.size(); ++i) {
        prob.x[i] = new svm_node[feature_num + 1];
        std::vector<double> feature = read_vec[i]->probability_proportions;
        for (int j = 0; j < feature_num; ++j) {
            prob.x[i][j].index = j + 1;
            prob.x[i][j].value = feature[j];
        }
        prob.x[i][feature_num].index = -1;
        prob.y[i] = read_vec[i]->label;
    }

    svm_model *model = svm_load_model(model_path.c_str());
    int32_t accurate_num = 0;
    int32_t pos_total_num = 0;
    int32_t neg_total_num = 0;
    int32_t true_pos_num = 0;
    int32_t true_neg_num = 0;
    double epsilon = 1e-6;
    std::cout << "Valid..." << std::endl;
    for (int i = 0; i < read_vec.size(); ++i) {
        svm_node test_node[feature_num + 1];
        std::vector<double> feature = read_vec[i]->probability_proportions;
        for (int j = 0; j < feature_num; ++j) {
            test_node[j].index = j + 1;
            test_node[j].value = feature[j];
        }
        test_node[feature_num].index = -1;
        double result = svm_predict(model, test_node);
        if (read_vec[i]->label == -1) {
            neg_total_num++;
            if (abs(result - read_vec[i]->label) < epsilon) {
                accurate_num++;
                true_neg_num++;
            }
        } else {
            pos_total_num++;
            if (abs(result - read_vec[i]->label) < epsilon) {
                accurate_num++;
                true_pos_num++;
            }
        }
    }
    float accuracy = (float) accurate_num / read_vec.size();
    std::cout << "Prediction accuracy: " << accuracy << std::endl;
    std::cout << "true pos count: " << true_pos_num << ", total: " << pos_total_num << std::endl;
    std::cout << "true neg count: " << true_neg_num << ", total: " << neg_total_num << std::endl;

    svm_free_and_destroy_model(&model);
    delete[] prob.y;
    for (int i = 0; i < read_vec.size(); ++i) {
        delete[] prob.x[i];
    }
    delete[] prob.x;
    return 0;
}

int predict(size_t num_threshold,
            int feature_num,
            std::filesystem::path input_bam_file,
            std::filesystem::path output_file,
            std::filesystem::path model_path) {
    svm_model *model = svm_load_model(model_path.c_str());
    samFile *bam_in = sam_open(input_bam_file.c_str(), "r");
    bam_hdr_t *bam_header = sam_hdr_read(bam_in);
    bam1_t *aln = bam_init1();
    std::ofstream out;
    out.open(output_file);
    double epsilon = 1e-6;
    while (sam_read1(bam_in, bam_header, aln) >= 0) {
        std::shared_ptr<BamRead> bam_ptr = std::make_shared<BamRead>(aln, 0, feature_num);
        if (bam_ptr->p_num < num_threshold || !(bam_ptr->is_valid)) {
            out << bam_ptr->query_name << "\t" << "*" << std::endl;
            continue;
        }
        svm_node node[feature_num + 1];
        std::vector<double> feature = bam_ptr->probability_proportions;
        for (int j = 0; j < feature_num; ++j) {
            node[j].index = j + 1;
            node[j].value = feature[j];
        }
        node[feature_num].index = -1;
        double result = svm_predict(model, node);
        if (abs(result - (-1)) < epsilon) {
            out << bam_ptr->query_name << "\t" << "-" << std::endl;
        } else {
            out << bam_ptr->query_name << "\t" << "+" << std::endl;
        }
    }
    out.close();
    return 0;

}