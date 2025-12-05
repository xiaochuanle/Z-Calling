#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <algorithm>

struct Record {
    int count_above_threshold = 0;
    int count_below_threshold = 0;
};

struct ChrCoor {
    std::string chr;
    int coor;

    bool operator<(const ChrCoor& other) const {
        if (chr == other.chr) {
            return coor < other.coor;
        }
        return chr < other.chr;
    }
};

struct ChrCoorRecord {
    ChrCoor chr_coor;
    Record record;
};

int main(int argc, char* argv[]) {
    std::string input_filename;
    std::string output_filename;
    double threshold = 0.5;

    // Parse command line arguments for input, output, and threshold
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            input_filename = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            output_filename = argv[++i];
        } else if (arg == "--threshold" && i + 1 < argc) {
            threshold = std::stod(argv[++i]);
        }
    }

    // Check if input and output filenames are provided
    if (input_filename.empty() || output_filename.empty()) {
        std::cerr << "Error: Input and output file names must be provided using -i and -o options." << std::endl;
        return 1;
    }

    std::ifstream input_file(input_filename);
    if (!input_file.is_open()) {
        std::cerr << "Error: Unable to open input file: " << input_filename << std::endl;
        return 1;
    }

    std::unordered_map<std::string, Record> chr_coor_map;
    std::string line;

    // Read input file line by line
    while (std::getline(input_file, line)) {
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        // Split line into columns
        while (std::getline(ss, token, '\t')) {
            columns.push_back(token);
        }

        // Check if there are enough columns
        // Mapped reads have 8 columns, Unmapped have 6. We only process mapped.
        if (columns.size() < 8) {
            continue;
        }

        // Extract relevant columns based on z-bam2txt output format:
        // Col 5: Reference Name (chr)
        // Col 6: Reference Coordinate (coor)
        // Col 7: Probability (prob)
        std::string chr = columns[5];
        int coor = std::stoi(columns[6]);
        double prob = std::stod(columns[7]);

        std::string chr_coor_key = chr + "_" + std::to_string(coor);
        
        // Update record counts based on probability threshold
        if (prob >= threshold) {
            chr_coor_map[chr_coor_key].count_above_threshold++;
        } else {
            chr_coor_map[chr_coor_key].count_below_threshold++;
        }
    }

    input_file.close();

    // Convert unordered_map to a vector for sorting
    std::vector<ChrCoorRecord> sorted_records;
    for (const auto& entry : chr_coor_map) {
        const std::string& chr_coor_key = entry.first;
        const Record& record = entry.second;

        // FIX: Use rfind (reverse find) to find the LAST underscore.
        // This ensures correct parsing even if chromosome names contain underscores 
        // (e.g., "chr1_gl000191_random_100").
        size_t delimiter_pos = chr_coor_key.rfind('_');
        
        if (delimiter_pos == std::string::npos) {
             std::cerr << "Warning: Malformed key found: " << chr_coor_key << std::endl;
             continue;
        }

        std::string chr = chr_coor_key.substr(0, delimiter_pos);
        int coor = std::stoi(chr_coor_key.substr(delimiter_pos + 1));

        sorted_records.push_back({{chr, coor}, record});
    }

    // Sort records by chr and coor
    std::sort(sorted_records.begin(), sorted_records.end(), [](const ChrCoorRecord& a, const ChrCoorRecord& b) {
        return a.chr_coor < b.chr_coor;
    });

    // Write output file
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Error: Unable to open output file: " << output_filename << std::endl;
        return 1;
    }

    output_file << "chr\tcoor\tnum_records_above_threshold\tnum_records_below_threshold\tratio_above_threshold\n";
    for (const auto& entry : sorted_records) {
        const ChrCoor& chr_coor = entry.chr_coor;
        const Record& record = entry.record;

        int total_records = record.count_above_threshold + record.count_below_threshold;
        double ratio = (total_records > 0) ? static_cast<double>(record.count_above_threshold) / total_records : 0.0;

        output_file << chr_coor.chr << "\t" << chr_coor.coor << "\t" << record.count_above_threshold << "\t"
                    << record.count_below_threshold << "\t" << std::fixed << std::setprecision(2) << ratio << "\n";
    }

    output_file.close();
    return 0;
}
