#include <iostream>
#include <fstream>
#include <string>

#include "fasta.h"

using namespace std;

Fasta readFasta(Fasta fasta, string filename) {
    std::ifstream input(filename);
    if (!input.good()) {
        cerr << "Error: can't open " << filename << endl;
        return fasta;
    }

    int c = 0;
    string line;
    while (getline(input, line).good()) {
        if (line[0] == '>') {
            c++;
        }
    }
    if (c == 0) {
        cerr << "Error: no headers detected in " << filename << endl;
        return fasta;
    }
    fasta.count = c;

    input.clear();
    input.seekg(0);
    fasta.headers = new string[fasta.count];
    fasta.sequences = new string[fasta.count];
    int sid = 0;
    string name, content;
    while (getline(input, line).good()) {
        if (line.empty() || line[0] == '>') {
            if (!name.empty()) {
                fasta.headers[sid] = name;
                fasta.sequences[sid] = content;
                sid++;
                name.clear();
            }
            if (!line.empty()) {
                name = line.substr(1);
            }
            content.clear();
        }
        else if (!name.empty()) {
            if (line.find(' ') != string::npos) {
                name.clear();
                content.clear();
            }
            else {
                content += line;
            }
        }
    }
    if (!name.empty()) {
        fasta.headers[sid] = name;
        fasta.sequences[sid] = content;
    }

    for (int i = 0; i < fasta.count; ++i) {
        cout << ">" << fasta.headers[i] << endl;
        cout << fasta.sequences[i] << endl;
    }

    return fasta;
}
