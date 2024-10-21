#include <iostream>
#include "BamRead.h"
#include "../3rdparty/argparse/argparse.h"
#include "func.h"


int main(int argc, char **argv) {
    argparse::ArgumentParser program("z-calling-read", "0.1.0");

    argparse::ArgumentParser train_args("train");
    train_args.add_argument("pos_bam_path");
    train_args.add_argument("neg_bam_path");
    train_args.add_argument("model_path");
    train_args.add_argument("length_threshold").scan<'i', int>();
    train_args.add_argument("section_num").scan<'i', int>();
    train_args.add_argument("C").default_value(1.0).scan<'f', float>();
    train_args.add_argument("gamma").default_value(0.5).scan<'f', float>();
    program.add_subparser(train_args);

    argparse::ArgumentParser valid_args("valid");
    valid_args.add_argument("pos_bam_path");
    valid_args.add_argument("neg_bam_path");
    valid_args.add_argument("model_path");
    valid_args.add_argument("length_threshold").scan<'i', int>();
    valid_args.add_argument("section_num").scan<'i', int>();
    program.add_subparser(valid_args);

    argparse::ArgumentParser predict_args("predict");
    predict_args.add_argument("input_bam_path");
    predict_args.add_argument("output_file_path");
    predict_args.add_argument("model_path");
    predict_args.add_argument("length_threshold").scan<'i', int>();
    predict_args.add_argument("section_num").scan<'i', int>();
    program.add_subparser(predict_args);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    if (program.is_subcommand_used("train")) {
        std::filesystem::path pos_bam_path = train_args.get<std::string>("pos_bam_path");
        std::filesystem::path neg_bam_path = train_args.get<std::string>("neg_bam_path");
        std::filesystem::path model_path = train_args.get<std::string>("model_path");
        size_t length_threshold = train_args.get<int32_t>("length_threshold");
        int feature_num = train_args.get<int>("section_num");
        float C = train_args.get<float>("C");
        float gamma = train_args.get<float>("gamma");
        train(length_threshold, feature_num, pos_bam_path, neg_bam_path, model_path, C, gamma);
    }

    if (program.is_subcommand_used("valid")) {
        std::filesystem::path pos_bam_path = valid_args.get<std::string>("pos_bam_path");
        std::filesystem::path neg_bam_path = valid_args.get<std::string>("neg_bam_path");
        std::filesystem::path model_path = valid_args.get<std::string>("model_path");
        size_t length_threshold = valid_args.get<int32_t>("length_threshold");
        int feature_num = valid_args.get<int>("section_num");
        valid(length_threshold, feature_num, pos_bam_path, neg_bam_path, model_path);
    }

    if (program.is_subcommand_used("predict")) {
        std::filesystem::path input_bam_path = predict_args.get<std::string>("input_bam_path");
        std::filesystem::path output_file_path = predict_args.get<std::string>("output_file_path");
        std::filesystem::path model_path = predict_args.get<std::string>("model_path");
        size_t length_threshold = predict_args.get<int32_t>("length_threshold");
        int feature_num = predict_args.get<int>("section_num");
        predict(length_threshold, feature_num, input_bam_path, output_file_path, model_path);
    }

    return 0;
}

