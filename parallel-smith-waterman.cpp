#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>
#include <unistd.h>
#include <random>
#include <new>

using namespace std;

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

// Adjust CHUNK_SIZE based on your cache line (e.g., 64 bytes) to manage cache-line alignment
const int CHUNK_SIZE = (int)(sysconf(_SC_LEVEL1_DCACHE_LINESIZE) / 2);

// Scoring parameters
const int MATCH = 2;
const int MISMATCH = -1;
const int GAP = -1;

void smithWatermanParallel(const string& seqA, const string& seqB) {
    int n = seqA.length();
    int m = seqB.length();

    // Changed short to int to prevent overflow with large strings and high match scores
    vector<vector<int>> H(n + 1, vector<int>(m + 1, 0));

    // Calculate number of chunks in each dimension
    int chunks_i = (n + CHUNK_SIZE - 1) / CHUNK_SIZE;
    int chunks_j = (m + CHUNK_SIZE - 1) / CHUNK_SIZE;

    // Total number of anti-diagonals in the tiled matrix
    int total_diagonals = chunks_i + chunks_j - 1;
    
    int globalMaxScore = 0;

    // Wavefront: Process diagonal by diagonal
    double start = omp_get_wtime();
    for (int d = 0; d < total_diagonals; ++d) {
        
        // Find which blocks belong to the current anti-diagonal 'd'
        int start_bi = max(0, d - chunks_j + 1);
        int end_bi = min(d, chunks_i - 1);

        #pragma omp parallel for schedule(dynamic)
        for (int bi = start_bi; bi <= end_bi; ++bi) {
            int bj = d - bi;

            // Calculate actual pixel boundaries for this chunk
            int row_start = bi * CHUNK_SIZE + 1;
            int row_end = min(row_start + CHUNK_SIZE, n + 1);
            
            int col_start = bj * CHUNK_SIZE + 1;
            int col_end = min(col_start + CHUNK_SIZE, m + 1);

            int localMax = 0;

            // Standard Smith-Waterman within the chunk
            for (int i = row_start; i < row_end; ++i) {
                for (int j = col_start; j < col_end; ++j) {
                    int score = (seqA[i - 1] == seqB[j - 1]) ? MATCH : MISMATCH;
                    
                    H[i][j] = max({
                        0,
                        H[i - 1][j - 1] + score, // Diagonal
                        H[i - 1][j] + GAP,       // Up
                        H[i][j - 1] + GAP        // Left
                    });
                    
                    if (H[i][j] > localMax) {
                        localMax = H[i][j];
                    }
                }
            }
            
            // Safely evaluate the highest score across threads
            if (localMax > globalMaxScore) {
                #pragma omp critical
                {
                    if (localMax > globalMaxScore) {
                        globalMaxScore = localMax;
                    }
                }
            }
        }
    }
    double end = omp_get_wtime();

    cout << "Computation complete." << endl;
    cout << "Alignment Score: " << globalMaxScore << endl;
    cout << "Execution time: " << end - start << " (sec)" << endl;
}

int main() {
    cout << "Chunk size: " << CHUNK_SIZE << endl;

    int str_size = 20000;
    string s1 = generateSequence(str_size, 1);
    string s2 = generateSequence(str_size, 2);

    smithWatermanParallel(s1, s2);
    
    return 0;
}