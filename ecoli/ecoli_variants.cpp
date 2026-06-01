#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <cctype>
#include <algorithm>
#include <chrono>
#include <set>
#include <iomanip>

// reverse complement of DNA sequence
std::string reverse_complement(const std::string &seq) {
    std::string rc;
    rc.reserve(seq.size());

    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (*it) {
            case 'A': rc += 'T'; break;
            case 'T': rc += 'A'; break;
            case 'C': rc += 'G'; break;
            case 'G': rc += 'C'; break;
            default:  rc += *it;
        }
    }

    return rc;
}

// load variants from CSV file into a set
std::set<std::tuple<char, int, char>>
load_variants(const std::string &filename) {

    std::set<std::tuple<char, int, char>> vars;
    std::ifstream file(filename);

    if (!file.is_open()) return vars;

    std::string line;

    // CSV format: type,pos,base
    while (std::getline(file, line)) {

        if (line.empty()) continue;

        std::stringstream ss(line);

        std::string t, p, b;
        std::getline(ss, t, ',');
        std::getline(ss, p, ',');
        std::getline(ss, b, ',');

        if (t.empty() || p.empty() || b.empty()) continue;

        vars.insert({t[0], std::stoi(p), b[0]});
    }

    return vars;
}

// compute TP/FP/FN metrics between predicted and truth sets
double compute_metrics(
    const std::set<std::tuple<char,int,char>> &predicted,
    const std::set<std::tuple<char,int,char>> &truth,
    int &TP, int &FP, int &FN
) {
    TP = FP = FN = 0;

    for (auto &v : predicted)
        (truth.count(v) ? TP : FP)++;

    for (auto &v : truth)
        if (!predicted.count(v)) FN++;

    return (TP + FP + FN)
        ? (double)TP / (TP + FP + FN)
        : 0.0;
}

int main() {

    // start timer for variant calling
    auto start = std::chrono::steady_clock::now();

    std::ifstream sam("data/aln.sam");

    if (!sam.is_open()) {
        std::cerr << "Ne mogu otvoriti aln.sam\n";
        return 1;
    }

    std::map<std::pair<int,char>, std::map<char,int>> variant_counts;
    std::map<int,int> coverage;

    std::string line;
    
    // compute coverage and collect variant evidence
    while (std::getline(sam, line)) {

        // skip SAM header and empty lines
        if (line.empty() || line[0] == '@') continue;

        std::istringstream ss(line);

        std::string qname, flag_str, rname, pos_str, mapq;
        std::string cigar, rnext, pnext, tlen;
        std::string seq, qual;

        ss >> qname >> flag_str >> rname >> pos_str >> mapq
           >> cigar >> rnext >> pnext >> tlen >> seq >> qual;

        if (cigar == "*" || rname == "*") continue;

        int flag = std::stoi(flag_str);

        // skip unmapped and secondary/supplementary alignments
        if (flag & 4 || flag & 256 || flag & 2048) continue;

        int pos = std::stoi(pos_str) - 1;
        int seq_idx = 0;

        bool rev = flag & 16;

        // reverse complement sequence if read is on reverse strand
        if (rev) seq = reverse_complement(seq);

        int tpos = pos;

        // first pass: compute coverage per reference position
        for (size_t i = 0; i < cigar.size();) {

            int num = 0;

            while (i < cigar.size() && isdigit(cigar[i]))
                num = num * 10 + (cigar[i++] - '0');

            if (i >= cigar.size()) break;

            char op = cigar[i++];

            if (op=='=' || op=='M' || op=='X') {
                for (int j=0;j<num;j++) coverage[tpos++]++;
                seq_idx += num;
            }
            else if (op=='D') {
                for (int j=0;j<num;j++) coverage[tpos++]++;
            }
            else if (op=='I' || op=='S') {
                seq_idx += num;
            }
        }

        seq_idx = 0;

        // second pass: collect variant evidence
        for (size_t i = 0; i < cigar.size();) {

            int num = 0;

            while (i < cigar.size() && isdigit(cigar[i]))
                num = num * 10 + (cigar[i++] - '0');

            if (i >= cigar.size()) break;

            char op = cigar[i++];

            // handle matches
            if (op=='=' || op=='M') {
                pos += num;
                seq_idx += num;
            }

            // handle mismatches, insertions, deletions
            else if (op=='X') {
                for (int j=0;j<num;j++) {
                    char base = std::toupper(seq[seq_idx]);
                    variant_counts[{pos,'X'}][base]++;
                    pos++; seq_idx++;
                }
            }

            else if (op=='I') {
                for (int j=0;j<num;j++) {
                    char base = std::toupper(seq[seq_idx]);
                    variant_counts[{pos - 1,'I'}][base]++;
                    seq_idx++;
                }
            }

            else if (op=='D') {
                for (int j=0;j<num;j++) {
                    variant_counts[{pos,'D'}]['-']++;
                    pos++;
                }
            }

            else if (op=='S') seq_idx += num;
        }
    }

    sam.close();

    std::vector<std::tuple<int,char,char>> variants;

    // convert raw counts into final variant calls
    for (auto &kv : variant_counts) {

        int pos = kv.first.first;
        char type = kv.first.second;

        int total = 0, best = 0;
        char best_base = 0;

        // find most common base and total count for this variant
        for (auto &bc : kv.second) {
            total += bc.second;
            if (bc.second > best) {
                best = bc.second;
                best_base = bc.first;
            }
        }

        int cov = coverage[pos];
        if (!cov) cov = 1;

        double freq = (double)total / cov;

        // call variant if frequency >= 70% and supported by at least 3 reads
        if (freq >= 0.7 && total >= 3)
            variants.push_back({pos,type,best_base});
    }

    std::set<std::tuple<char,int,char>> predicted;

    for (auto &v : variants)
        predicted.insert({std::get<1>(v), std::get<0>(v), std::get<2>(v)});

    // end timer for variant calling
    auto end = std::chrono::steady_clock::now();

    double time_sec =
        std::chrono::duration<double>(end - start).count();
        
    // load truth sets and compute metrics
    auto truth1 = load_variants("stmpileup_var.csv");

    int TP1, FP1, FN1;
    double acc1 = compute_metrics(predicted, truth1, TP1, FP1, FN1);

    double precision1 = TP1 + FP1 ? (double)TP1 / (TP1 + FP1) : 0;
    double recall1    = TP1 + FN1 ? (double)TP1 / (TP1 + FN1) : 0;
    double f11        = (precision1 + recall1) ? 2 * precision1 * recall1 / (precision1 + recall1) : 0;

    auto truth2 = load_variants("data/ecoli_mutated.csv");

    int TP2, FP2, FN2;
    double acc2 = compute_metrics(predicted, truth2, TP2, FP2, FN2);

    double precision2 = TP2 + FP2 ? (double)TP2 / (TP2 + FP2) : 0;
    double recall2    = TP2 + FN2 ? (double)TP2 / (TP2 + FN2) : 0;
    double f12        = (precision2 + recall2) ? 2 * precision2 * recall2 / (precision2 + recall2) : 0;

    std::cout << std::fixed << std::setprecision(2);

    std::cout << "=== TIME (variant calling) ===\n";
    std::cout << "Time: " << time_sec << " s\n\n";

    std::cout << "=== Check (stmpileup_var) ===\n";
    std::cout << "TP: " << TP1 << " FP: " << FP1 << " FN: " << FN1 << "\n";
    std::cout << "Accuracy: " << acc1 * 100 << "%\n";
    std::cout << "Precision: " << precision1 * 100 << "%\n";
    std::cout << "Recall: " << recall1 * 100 << "%\n";
    std::cout << "F1: " << f11 * 100 << "%\n\n";

    std::cout << "=== Ground Truth (ecoli_mutated) ===\n";
    std::cout << "TP: " << TP2 << " FP: " << FP2 << " FN: " << FN2 << "\n";
    std::cout << "Accuracy: " << acc2 * 100 << "%\n";
    std::cout << "Precision: " << precision2 * 100 << "%\n";
    std::cout << "Recall: " << recall2 * 100 << "%\n";
    std::cout << "F1: " << f12 * 100 << "%\n";

    return 0;
}
