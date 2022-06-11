#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <time.h>

#include "fasta.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

void to_index(int& i, int count, int* dim, int rdim, int* coords) {
    i = 0;
    rdim = 1;
    for (int s = 0; s < count; ++s) {
        i += coords[s]*rdim;
        rdim *= dim[s];
    }
}

void to_index(int& i, int count, int* dim, int rdim, int* coords, int** bin, int a) {
    i = 0;
    rdim = 1;
    for (int s = 0; s < count; ++s) {
        i += (coords[s] - bin[a][s]) * rdim;
        rdim *= dim[s];
    }
}

void get_coords(int count, int* dim, int index, int denom, int* coords) {
    denom = 1;
    for (int i = 0; i < count; ++i) {
        coords[i] = (index/denom)%dim[i];
        denom = denom * dim[i];
    }
}

int are_equal(int count, int score, string sequences[], int* coords) {
    char c = sequences[0][coords[0]];
    for (int i = 1; i < count; ++i) {
        if (sequences[i][coords[i]] != c) { return -score; }
    }
    return score;
}

Fasta readFasta(Fasta fasta, string filename);

int main(int argc, char** argv) {
    if (argc <= 1) {
        cerr << "Usage: " << argv[0] << " [infile] gap_penalty(1) similarity_score(5)" << endl;
        return -1;
    }

    int gap = 1;
    int sco = 5;
    if (argc > 2) {
        gap = atoi(argv[2]);
        if (argc > 3) {
            sco = atoi(argv[3]);
        }
    }

    Fasta fasta{};
    fasta = readFasta(fasta, argv[1]);

    int total = 1;
    int* dim = new int[fasta.count];
    cout << "dimensions: ";
    for (int i = 0; i < fasta.count; ++i) {
        dim[i] = fasta.sequences[i].length();
        total *= dim[i];
        cout << dim[i] << "*";
    }
    cout << " = " << total << endl;

    int* swm = new int[total];
    for (int i = 0; i < total; ++i) {
        swm[total] = 0;
    }

    int score_max = 0;
    int max_index = 0;

    int d = 0;
    int array_size = pow(2, fasta.count);
    int* max_array = new int[array_size];
    max_array[0] = 0;

    int** bin = new int* [array_size];
    for (int i = 0; i < array_size; ++i) {
        bin[i] = new int[fasta.count];
        for (int b = 0; b < fasta.count; ++b) {
            bin[i][b] = ((i & 1 << b) > 0);
        }
    }

    int ret = 0;
    int rdim = 1;
    int denom = 1;
    int* coords = new int[fasta.count];
    int* prev_coords = new int[fasta.count];
    bool border = false;
    int index = 0;
    for (int c = 0; c < total; ++c) {
        get_coords(fasta.count, dim, c, denom, coords);
        
        for (int s = 0; s < fasta.count; ++s) {
            prev_coords[s] = coords[s] - 1;
            if ( coords[s] == 0 ) {
                border = true;
            }
        }
        
        if (border) {border = false; continue; }
        
        for (int a = 1; a < array_size - 1; ++a) {
            to_index(index, fasta.count, dim, rdim, coords, bin, a);
            max_array[a] = swm[index] - gap;
        }
        to_index(index, fasta.count, dim, rdim, coords, bin, array_size - 1);
        max_array[array_size-1] = swm[index] + are_equal(fasta.count, sco, fasta.sequences, prev_coords);
        
        swm[c] = *max_element(max_array, max_array + array_size);
        if (score_max < swm[c]) {
            score_max = swm[c];
            max_index = c;
        }
    }

    int c = max_index;
    int* cpos = new int[fasta.count];
    get_coords(fasta.count, dim, c, denom, cpos);

    int strl = *max_element(cpos, cpos + fasta.count) + 1;

    char** nsequences = new char*[fasta.count];
    for (int i = 0; i < fasta.count; ++i) {
        nsequences[i] = (char*)malloc(strl);
        if (!nsequences[i]) {
            cerr << "Error: Could not allocate memory" << endl;
            return 1;
        }
        for (int n = 0; n < strl; ++n) {
            nsequences[i][n] = '-';
        }
        nsequences[i][strl] = '\0';
    }

    int l = strl-1;
    int score = score_max;
    int m = 0;
    int nm = 0;
    int ma = 0;
    while (l > -1) {
        if (score == 0) {
            for (int s = 0; s < fasta.count; ++s) {
                if (bin[ma][s] == 1) {
                    cpos[s] = cpos[s] - 1;
                    nsequences[s][l] = fasta.sequences[s][cpos[s]];
                }
                else {
                    nsequences[s][l] = '-';
                }
            }
            break;
        }
        
        for (int s = 0; s < fasta.count; ++s) {
            cout << cpos[s] << " ";
        }
        for (int s = 0; s < fasta.count; ++s) {
            cout << fasta.sequences[s][cpos[s]] << " ";
        }
        cout << score << " " << ma << " " << endl;

        m = 0;
        nm = 0;
        for (int a = array_size-1; a > 0; --a) {
            to_index(index, fasta.count, dim, rdim, cpos, bin, a);
            nm = max(m, swm[index]);
            if (m < nm) {
                m = nm;
                ma = a;
            }
        }

        for (int s = 0; s < fasta.count; ++s) {
            if (bin[ma][s] == 1) {
                cpos[s] = cpos[s] - 1;
                nsequences[s][l] = fasta.sequences[s][cpos[s]];
            }
            else {
                nsequences[s][l] = '-';
            }
        }
        to_index(c, fasta.count, dim, rdim, cpos);
        score = swm[c];
        --l; 
    }

    for (int i = 0; i < fasta.count; ++i) {
        cout << ">" << fasta.headers[i] << "-align" << endl;
        cout << nsequences[i] << endl;
    }
    return 0;
}
