#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// Scoring constants
const int MATCH    = 2;
const int MISMATCH = -1;
const int GAP      = -1;

void smithWaterman(string seq1, string seq2) {
    int m = seq1.length();
    int n = seq2.length();
    
    // Create (m+1) x (n+1) matrix initialized to 0
    vector<vector<int>> score(m + 1, vector<int>(n + 1, 0));
    
    int maxScore = 0;
    int maxI = 0, maxJ = 0;

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

            // Track highest score for traceback start
            if (score[i][j] > maxScore) {
                maxScore = score[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    cout << "Alignment Score: " << maxScore << endl;

    // 2. Traceback (Simplified)
    string align1 = "", align2 = "";
    int i = maxI, j = maxJ;

    while (i > 0 && j > 0 && score[i][j] > 0) {
        int current = score[i][j];
        int diag = score[i - 1][j - 1];
        int up = score[i - 1][j];
        int left = score[i][j - 1];

        if (current == diag + ((seq1[i - 1] == seq2[j - 1]) ? MATCH : MISMATCH)) {
            align1 = seq1[i - 1] + align1;
            align2 = seq2[j - 1] + align2;
            i--; j--;
        } else if (current == up + GAP) {
            align1 = seq1[i - 1] + align1;
            align2 = "-" + align2;
            i--;
        } else {
            align1 = "-" + align1;
            align2 = seq2[j - 1] + align2;
            j--;
        }
    }

    cout << "Local Alignment:" << endl;
    cout << align1 << endl;
    cout << align2 << endl;
}

int main() {
    string s1 = "GGTTAGACTTTGATGAAAAAAAAATGATAAAAAA";
    string s2 = "GGTTUACGTAGGAAAGGAGAAAA";
    
    smithWaterman(s1, s2);
    
    return 0;
}