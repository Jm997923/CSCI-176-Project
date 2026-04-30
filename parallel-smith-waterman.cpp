#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>
#include <unistd.h>
#include <random>
#include <new>

using namespace std;

std::string generateSequence(int length) {
    string bases = "ATGC";
    string seq = "";
    mt19937 rng(12345); // Fixed seed for reproducibility
    uniform_int_distribution<int> dist(0, 3);
    for (int i = 0; i < length; ++i) {
        seq += bases[dist(rng)];
    }
    return seq;
}

// Adjust CHUNK_SIZE based on your cache line (e.g., 64 bytes)
// If using ints (4 bytes), 16 would fit exactly in a 64-byte line.
// However, for Smith-Waterman, slightly larger blocks (32-64) 
// often perform better due to reduced scheduling overhead.
const int CHUNK_SIZE = (int)(sysconf(_SC_LEVEL1_DCACHE_LINESIZE) / 2);

// Scoring parameters
const int MATCH = 2;
const int MISMATCH = -1;
const int GAP = -1;

void smithWatermanParallel(const string& seqA, const string& seqB) {
    int n = seqA.length();
    int m = seqB.length();

    // Score matrix (n+1 x m+1)
    vector<vector<short>> H(n + 1, vector<short>(m + 1, 0));

    // Calculate number of chunks in each dimension
    int chunks_i = (n + CHUNK_SIZE - 1) / CHUNK_SIZE;
    int chunks_j = (m + CHUNK_SIZE - 1) / CHUNK_SIZE;

    // Total number of anti-diagonals in the tiled matrix
    int total_diagonals = chunks_i + chunks_j - 1;

    // Wavefront: Process diagonal by diagonal
    double start = omp_get_wtime();
    for (int d = 0; d < total_diagonals; ++d) {
        
        // Find which blocks belong to the current anti-diagonal 'd'
        // Block coordinates (bi, bj) must satisfy: bi + bj = d
        int start_bi = max(0, d - chunks_j + 1);
        int end_bi = min(d, chunks_i - 1);

        // This loop handles the "ramp up/down" of threads automatically.
        // On d=0, only 1 block is available. In the middle, many are.
        #pragma omp parallel for schedule(dynamic)
        for (int bi = start_bi; bi <= end_bi; ++bi) {
            int bj = d - bi;

            // Calculate actual pixel boundaries for this chunk
            int row_start = bi * CHUNK_SIZE + 1;
            int row_end = min(row_start + CHUNK_SIZE, n + 1);
            
            int col_start = bj * CHUNK_SIZE + 1;
            int col_end = min(col_start + CHUNK_SIZE, m + 1);

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
                }
            }
        }
    }
    double end = omp_get_wtime();

    // Results would be extracted here (Traceback or Max Score)
    cout << "Computation complete." << endl;
    cout << "Execution time: " << end - start << " (sec)" << endl;
}

int main() {
    cout << "Chunk size: " << CHUNK_SIZE << endl;

    int str_size = 20000;
    string s1 = generateSequence(str_size);
    string s2 = generateSequence(str_size);

    smithWatermanParallel(s1, s2);
    
    return 0;
}