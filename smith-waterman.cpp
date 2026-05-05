/**
 * Julian Mendoza & Noah Villanueva -- serial code for Smith Waterman algorithm
 * 
 * This program is meant to serve as the baseline for the Smith Waterman algorithm
 * where it will be compared with the parallelized version of the algorithm to compare
 * performance outcome.
 * 
 * $> g++ smith-waterman.cpp -o xxx
 * $> ./xxx n
 *  Where n is the size of the random strings that will be comapared (e.g. 20000)
 *  (e.g. ./xxx 20000)
 * Thank you!
 */

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <random>

using namespace std;

// Scoring constants
const int MATCH    = 2;
const int MISMATCH = -1;
const int GAP      = -1;

// Ported random string generator
std::string generateSequence(int length, int id) {
    string bases = "ATGC";
    string seq = "";
    mt19937 rng(12345 + id); // Fixed seed for reproducibility
    uniform_int_distribution<int> dist(0, 3);
    for (int i = 0; i < length; ++i) {
        seq += bases[dist(rng)];
    }
    return seq;
}

void smithWaterman(const string& seq1, const string& seq2) {
    int m = seq1.length();
    int n = seq2.length();
    
    // Create (m+1) x (n+1) matrix initialized to 0
    vector<vector<int>> score(m + 1, vector<int>(n + 1, 0));
    int maxScore = 0;

    // Start High-Resolution Timer
    auto start = chrono::high_resolution_clock::now();

    // 1. Fill the matrix
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int scoreSub = (seq1[i - 1] == seq2[j - 1]) ? MATCH : MISMATCH;
            
            score[i][j] = max({
                0,                                   // Floor at zero
                score[i - 1][j - 1] + scoreSub,      // Diagonal
                score[i - 1][j] + GAP,               // Up
                score[i][j - 1] + GAP                // Left
            });

            // Track highest score
            if (score[i][j] > maxScore) {
                maxScore = score[i][j];
            }
        }
    }

    // End High-Resolution Timer
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    cout << "Computation complete." << endl;
    cout << "Alignment Score: " << maxScore << endl;
    cout << "Execution time: " << duration.count() << " (sec)" << endl;

    // 2. Traceback (Removed to prevent I/O bottlenecks during benchmarking)
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Error: ./filename n" <<  endl;
        return -1;
    }

    int n = atoi(argv[1]);

    string s1 = generateSequence(n, 1);
    string s2 = generateSequence(n, 2);
    
    smithWaterman(s1, s2);
    
    return 0;
}