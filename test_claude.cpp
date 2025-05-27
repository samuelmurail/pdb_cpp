#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <iomanip>
#include <climits>

class NeedlemanWunsch {
private:
    // BLOSUM62 substitution matrix
    std::unordered_map<char, std::unordered_map<char, int>> blosum62 = {
        {'A', {{'A', 4}, {'R', -1}, {'N', -2}, {'D', -2}, {'C', 0}, {'Q', -1}, {'E', -1}, {'G', 0}, {'H', -2}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1}, {'T', 0}, {'W', -3}, {'Y', -2}, {'V', 0}}},
        {'R', {{'A', -1}, {'R', 5}, {'N', 0}, {'D', -2}, {'C', -3}, {'Q', 1}, {'E', 0}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -2}, {'K', 2}, {'M', -1}, {'F', -3}, {'P', -2}, {'S', -1}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -3}}},
        {'N', {{'A', -2}, {'R', 0}, {'N', 6}, {'D', 1}, {'C', -3}, {'Q', 0}, {'E', 0}, {'G', 0}, {'H', 1}, {'I', -3}, {'L', -3}, {'K', 0}, {'M', -2}, {'F', -3}, {'P', -2}, {'S', 1}, {'T', 0}, {'W', -4}, {'Y', -2}, {'V', -3}}},
        {'D', {{'A', -2}, {'R', -2}, {'N', 1}, {'D', 6}, {'C', -3}, {'Q', 0}, {'E', 2}, {'G', -1}, {'H', -1}, {'I', -3}, {'L', -4}, {'K', -1}, {'M', -3}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -3}}},
        {'C', {{'A', 0}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', 9}, {'Q', -3}, {'E', -4}, {'G', -3}, {'H', -3}, {'I', -1}, {'L', -1}, {'K', -3}, {'M', -1}, {'F', -2}, {'P', -3}, {'S', -1}, {'T', -1}, {'W', -2}, {'Y', -2}, {'V', -1}}},
        {'Q', {{'A', -1}, {'R', 1}, {'N', 0}, {'D', 0}, {'C', -3}, {'Q', 5}, {'E', 2}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -2}, {'K', 1}, {'M', 0}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -2}, {'Y', -1}, {'V', -2}}},
        {'E', {{'A', -1}, {'R', 0}, {'N', 0}, {'D', 2}, {'C', -4}, {'Q', 2}, {'E', 5}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -3}, {'K', 1}, {'M', -2}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}}},
        {'G', {{'A', 0}, {'R', -2}, {'N', 0}, {'D', -1}, {'C', -3}, {'Q', -2}, {'E', -2}, {'G', 6}, {'H', -2}, {'I', -4}, {'L', -4}, {'K', -2}, {'M', -3}, {'F', -3}, {'P', -2}, {'S', 0}, {'T', -2}, {'W', -2}, {'Y', -3}, {'V', -3}}},
        {'H', {{'A', -2}, {'R', 0}, {'N', 1}, {'D', -1}, {'C', -3}, {'Q', 0}, {'E', 0}, {'G', -2}, {'H', 8}, {'I', -3}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -1}, {'P', -2}, {'S', -1}, {'T', -2}, {'W', -2}, {'Y', 2}, {'V', -3}}},
        {'I', {{'A', -1}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -3}, {'E', -3}, {'G', -4}, {'H', -3}, {'I', 4}, {'L', 2}, {'K', -3}, {'M', 1}, {'F', 0}, {'P', -3}, {'S', -2}, {'T', -1}, {'W', -3}, {'Y', -1}, {'V', 3}}},
        {'L', {{'A', -1}, {'R', -2}, {'N', -3}, {'D', -4}, {'C', -1}, {'Q', -2}, {'E', -3}, {'G', -4}, {'H', -3}, {'I', 2}, {'L', 4}, {'K', -2}, {'M', 2}, {'F', 0}, {'P', -3}, {'S', -2}, {'T', -1}, {'W', -2}, {'Y', -1}, {'V', 1}}},
        {'K', {{'A', -1}, {'R', 2}, {'N', 0}, {'D', -1}, {'C', -3}, {'Q', 1}, {'E', 1}, {'G', -2}, {'H', -1}, {'I', -3}, {'L', -2}, {'K', 5}, {'M', -1}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}}},
        {'M', {{'A', -1}, {'R', -1}, {'N', -2}, {'D', -3}, {'C', -1}, {'Q', 0}, {'E', -2}, {'G', -3}, {'H', -2}, {'I', 1}, {'L', 2}, {'K', -1}, {'M', 5}, {'F', 0}, {'P', -2}, {'S', -1}, {'T', -1}, {'W', -1}, {'Y', -1}, {'V', 1}}},
        {'F', {{'A', -2}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -2}, {'Q', -3}, {'E', -3}, {'G', -3}, {'H', -1}, {'I', 0}, {'L', 0}, {'K', -3}, {'M', 0}, {'F', 6}, {'P', -4}, {'S', -2}, {'T', -2}, {'W', 1}, {'Y', 3}, {'V', -1}}},
        {'P', {{'A', -1}, {'R', -2}, {'N', -2}, {'D', -1}, {'C', -3}, {'Q', -1}, {'E', -1}, {'G', -2}, {'H', -2}, {'I', -3}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -4}, {'P', 7}, {'S', -1}, {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -2}}},
        {'S', {{'A', 1}, {'R', -1}, {'N', 1}, {'D', 0}, {'C', -1}, {'Q', 0}, {'E', 0}, {'G', 0}, {'H', -1}, {'I', -2}, {'L', -2}, {'K', 0}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 4}, {'T', 1}, {'W', -3}, {'Y', -2}, {'V', -2}}},
        {'T', {{'A', 0}, {'R', -1}, {'N', 0}, {'D', -1}, {'C', -1}, {'Q', -1}, {'E', -1}, {'G', -2}, {'H', -2}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1}, {'T', 5}, {'W', -2}, {'Y', -2}, {'V', 0}}},
        {'W', {{'A', -3}, {'R', -3}, {'N', -4}, {'D', -4}, {'C', -2}, {'Q', -2}, {'E', -3}, {'G', -2}, {'H', -2}, {'I', -3}, {'L', -2}, {'K', -3}, {'M', -1}, {'F', 1}, {'P', -4}, {'S', -3}, {'T', -2}, {'W', 11}, {'Y', 2}, {'V', -3}}},
        {'Y', {{'A', -2}, {'R', -2}, {'N', -2}, {'D', -3}, {'C', -2}, {'Q', -1}, {'E', -2}, {'G', -3}, {'H', 2}, {'I', -1}, {'L', -1}, {'K', -2}, {'M', -1}, {'F', 3}, {'P', -3}, {'S', -2}, {'T', -2}, {'W', 2}, {'Y', 7}, {'V', -1}}},
        {'V', {{'A', 0}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -2}, {'E', -2}, {'G', -3}, {'H', -3}, {'I', 3}, {'L', 1}, {'K', -2}, {'M', 1}, {'F', -1}, {'P', -2}, {'S', -2}, {'T', 0}, {'W', -3}, {'Y', -1}, {'V', 4}}}
    };
    
    int gapOpenPenalty = -10;    // Penalty for opening a new gap
    int gapExtensionPenalty = -1; // Penalty for extending an existing gap
    
    int getScore(char a, char b) {
        if (blosum62.find(a) != blosum62.end() && blosum62[a].find(b) != blosum62[a].end()) {
            return blosum62[a][b];
        }
        return -4; // Default penalty for unknown characters
    }

public:
    struct AlignmentResult {
        std::string seq1Aligned;
        std::string seq2Aligned;
        int score;
    };
    
    AlignmentResult align(const std::string& seq1, const std::string& seq2) {
        int m = seq1.length();
        int n = seq2.length();
        
        // Create three matrices for affine gap penalties
        // M[i][j] = best score ending with match/mismatch
        // Ix[i][j] = best score ending with gap in seq1 (insertion in seq2)
        // Iy[i][j] = best score ending with gap in seq2 (deletion from seq1)
        std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, INT_MIN));
        std::vector<std::vector<int>> Ix(m + 1, std::vector<int>(n + 1, INT_MIN));
        std::vector<std::vector<int>> Iy(m + 1, std::vector<int>(n + 1, INT_MIN));
        
        // Initialize base cases
        M[0][0] = 0;
        
        // Initialize first row (gaps in seq1)
        for (int j = 1; j <= n; j++) {
            Ix[0][j] = gapOpenPenalty + (j - 1) * gapExtensionPenalty;
        }
        
        // Initialize first column (gaps in seq2)
        for (int i = 1; i <= m; i++) {
            Iy[i][0] = gapOpenPenalty + (i - 1) * gapExtensionPenalty;
        }
        
        // Fill the matrices
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                // Match/mismatch matrix
                M[i][j] = getScore(seq1[i-1], seq2[j-1]) + std::max({M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]});
                
                // Gap in seq1 (insertion in seq2)
                Ix[i][j] = std::max(M[i][j-1] + gapOpenPenalty, Ix[i][j-1] + gapExtensionPenalty);
                
                // Gap in seq2 (deletion from seq1)
                Iy[i][j] = std::max(M[i-1][j] + gapOpenPenalty, Iy[i-1][j] + gapExtensionPenalty);
            }
        }
        
        // Find the best final score
        int finalScore = std::max({M[m][n], Ix[m][n], Iy[m][n]});
        
        // Traceback to find the alignment
        std::string aligned1 = "";
        std::string aligned2 = "";
        
        int i = m, j = n;
        
        // Determine which matrix has the optimal score
        char currentMatrix;
        if (M[m][n] == finalScore) currentMatrix = 'M';
        else if (Ix[m][n] == finalScore) currentMatrix = 'X';
        else currentMatrix = 'Y';
        
        while (i > 0 || j > 0) {
            if (currentMatrix == 'M') {
                if (i > 0 && j > 0) {
                    aligned1 = seq1[i-1] + aligned1;
                    aligned2 = seq2[j-1] + aligned2;
                    
                    // Determine previous matrix
                    int prevScore = M[i][j] - getScore(seq1[i-1], seq2[j-1]);
                    if (M[i-1][j-1] == prevScore) currentMatrix = 'M';
                    else if (Ix[i-1][j-1] == prevScore) currentMatrix = 'X';
                    else currentMatrix = 'Y';
                    
                    i--; j--;
                } else break;
            }
            else if (currentMatrix == 'X') {
                // Gap in seq1
                aligned1 = '-' + aligned1;
                aligned2 = seq2[j-1] + aligned2;
                
                // Determine if gap continues or starts new
                if (j > 1 && Ix[i][j-1] + gapExtensionPenalty == Ix[i][j]) {
                    currentMatrix = 'X'; // Continue gap
                } else {
                    currentMatrix = 'M'; // Gap started from match
                }
                j--;
            }
            else { // currentMatrix == 'Y'
                // Gap in seq2
                aligned1 = seq1[i-1] + aligned1;
                aligned2 = '-' + aligned2;
                
                // Determine if gap continues or starts new
                if (i > 1 && Iy[i-1][j] + gapExtensionPenalty == Iy[i][j]) {
                    currentMatrix = 'Y'; // Continue gap
                } else {
                    currentMatrix = 'M'; // Gap started from match
                }
                i--;
            }
        }
        
        return {aligned1, aligned2, finalScore};
    }
    
    void printAlignment(const AlignmentResult& result) {
        std::cout << "Alignment Score: " << result.score << std::endl;
        std::cout << "Sequence 1: " << result.seq1Aligned << std::endl;
        std::cout << "Sequence 2: " << result.seq2Aligned << std::endl;
        std::cout << "Matches:    ";
        
        for (size_t i = 0; i < result.seq1Aligned.length(); i++) {
            if (result.seq1Aligned[i] == result.seq2Aligned[i]) {
                std::cout << '|';
            } else if (result.seq1Aligned[i] == '-' || result.seq2Aligned[i] == '-') {
                std::cout << ' ';
            } else {
                std::cout << ':';
            }
        }
        std::cout << std::endl;
    }
    
    // Setter methods to customize gap penalties
    void setGapPenalties(int openPenalty, int extensionPenalty) {
        gapOpenPenalty = openPenalty;
        gapExtensionPenalty = extensionPenalty;
    }
};

int main() {
    NeedlemanWunsch nw;
    
    // You can customize gap penalties
    nw.setGapPenalties(-10, -1); // Gap open: -10, Gap extension: -1
    
    // Example usage
    std::string seq1 = "AQDMVSPPXPIADEPLTVXSLSWKDRRL";
    std::string seq2 = "AQDMVSPPPPIADEPLTVNTGIYLIECYSLDDKAETFKVNAFLSLSWKDRRLAFDPV";
    
    

    std::cout << "Gap Open Penalty: -10, Gap Extension Penalty: -1" << std::endl;
    std::cout << "Input sequences:" << std::endl;
    std::cout << "Sequence 1: " << seq1 << std::endl;
    std::cout << "Sequence 2: " << seq2 << std::endl;
    std::cout << std::endl;
    
    auto result = nw.align(seq1, seq2);
    nw.printAlignment(result);
    
    std::cout << std::endl << "Another example with different gap penalties:" << std::endl;
    
    // Try with different gap penalties
    nw.setGapPenalties(-5, -2);
    std::cout << "Gap Open Penalty: -5, Gap Extension Penalty: -2" << std::endl;
    
    seq1 = "ACGTACGT";
    seq2 = "ACGACGT";
    
    std::cout << "Input sequences:" << std::endl;
    std::cout << "Sequence 1: " << seq1 << std::endl;
    std::cout << "Sequence 2: " << seq2 << std::endl;
    std::cout << std::endl;
    
    result = nw.align(seq1, seq2);
    nw.printAlignment(result);
    
    return 0;
}