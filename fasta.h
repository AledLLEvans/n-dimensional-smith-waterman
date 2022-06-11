#pragma once

#include <string>

using namespace std;

struct Fasta {
    string* headers;
    string* sequences;
    int count;
};