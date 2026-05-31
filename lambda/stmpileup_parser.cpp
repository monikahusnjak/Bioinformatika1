#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <cstring>

using namespace std;

// Represents a detected genomic variant.
struct Variant {

    // Variant type:
    // X = SNP/substitution, I = insertion, D = deletion
    char type;

    // 0-based genomic position
    int pos;

    // Alternative base/sequence or "-" for deletions
    string alt;
};

int main() {

    // Open mpileup output file for reading
    ifstream infile("pileup.txt");

    if (!infile.is_open()) {
        cerr << "Cannot open pileup.txt\n";
        return 1;
    }

    // Stores all detected variants
    vector<Variant> variants;

    string line;

    // Process mpileup file line by line
    while (getline(infile, line)) {

        // Skip empty lines
        if (line.empty())
            continue;

        stringstream ss(line);

        string chrom;
        int pos1;           // 1-based position from mpileup
        char ref;
        int depth;
        string bases;
        string quals;

        // Parse standard mpileup columns
        ss >> chrom >> pos1 >> ref >> depth >> bases >> quals;

        // Convert genomic position to 0-based indexing
        int pos0 = pos1 - 1;

        // Skip positions without read coverage
        if (depth == 0)
            continue;

        // Stores support counts for each variant type
        map<char, int> snpCount;
        map<string, int> insertionCount;
        map<string, int> deletionCount;

        // Parse the mpileup read bases field
        // Mpileup format specification is used here
        for (size_t i = 0; i < bases.size(); i++) {

            char c = bases[i];

            // '^' marks the start of a read segment
            // The next character stores mapping quality
            if (c == '^') {
                i++;
                continue;
            }

            // '$' marks the end of a read segment
            if (c == '$') {
                continue;
            }

            // '.' and ',' represent reference matches
            if (c == '.' || c == ',') {
                continue;
            }

            // '*' represents a deleted base continuation
            if (c == '*') {
                continue;
            }

            // Process SNP candidates
            if (strchr("ACGTacgt", c)) {

                char base = toupper(c);

                // Count only bases different from reference
                if (base != ref) {
                    snpCount[base]++;
                }

                continue;
            }

            // Process insertion or deletion events
            if (c == '+' || c == '-') {

                bool isInsertion = (c == '+');

                i++;

                // Parse indel length
                string number;

                while (i < bases.size() && isdigit(bases[i])) {
                    number += bases[i];
                    i++;
                }

                if (number.empty())
                    continue;

                int len = stoi(number);

                // Parse inserted/deleted sequence
                string seq;

                for (int j = 0; j < len && i < bases.size(); j++, i++) {
                    seq += toupper(bases[i]);
                }

                // Step back because the outer loop increments i
                i--;

                // Count insertion or deletion support
                if (isInsertion) {
                    insertionCount[seq]++;
                } else {
                    deletionCount[seq]++;
                }
            }
        }

        // Variant is accepted if at least 70% of reads support it

        // SNP calling
        for (auto &p : snpCount) {

            int support = p.second;

            if ((double)support / depth >= 0.7) {

                Variant v;
                v.type = 'X';
                v.pos = pos0;
                v.alt = string(1, p.first);

                variants.push_back(v);
            }
        }

        // Insertion calling
        for (auto &p : insertionCount) {

            int support = p.second;

            if ((double)support / depth >= 0.7) {

                Variant v;
                v.type = 'I';
                v.pos = pos0;
                v.alt = p.first;

                variants.push_back(v);
            }
        }

        // Deletion calling
        for (auto &p : deletionCount) {

            int support = p.second;

            if ((double)support / depth >= 0.7) {

                // Output one deletion entry per deleted base
                // Each deleted base gets its own genomic position
                for (size_t k = 0; k < p.first.size(); k++) {

                    Variant v;
                    v.type = 'D';
                    v.pos = pos0 + 1 + (int)k;
                    v.alt = "-";

                    variants.push_back(v);
                }
            }
        }
    }

    infile.close();

    // Create CSV output file
    ofstream outfile("stmpileup_var.csv");

    if (!outfile.is_open()) {
        cerr << "Cannot create stmpileup_var.csv\n";
        return 1;
    }

    // Write all detected variants to CSV
    for (const auto &v : variants) {

        outfile << v.type
                << ","
                << v.pos
                << ","
                << v.alt
                << "\n";
    }

    outfile.close();

    cout << "Output written to stmpileup_var.csv\n";

    return 0;
}