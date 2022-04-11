#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <random>
#include <iomanip>

#include "cryptocontextgen.h"
#include "palisade.h"

using namespace std;


int64_t bin2Int(vector<int64_t> binVect);
bool areEqual(vector<int64_t> v, vector<int64_t> u);
vector<int> vectIntPermutation(int seed, int begin, int end);
map<int, vector<int>> generatePermutations(int seed, size_t nPerm, size_t len);
const vector<string> explode(const string& s, const char& c); 
map<int, vector<vector<int64_t>>> readIrisCodeAsBlocks(string filename, size_t nB, size_t lenB);
map<int, vector<int64_t>> readIrisCodeAsLines(string filename);
map<int, vector<int64_t>> genKeys(size_t hb, size_t nBlocks);
void writeProcTime(string filename, double processingTime);
vector<vector<int64_t> > transpose2DVect(vector<vector<int64_t> > v);
map<int, vector<vector<int64_t>>> splitRegions( map<int, vector<int64_t>> sampleLines);
map<int, vector<vector<int64_t>>> permuteRegions(map<int, vector<vector<int64_t>>> regions, map<int, vector<int>> permutations);
map<int, vector<int64_t>> joinRegions( map<int, vector<vector<int64_t>>> regions);
map<int, vector<vector<int64_t>>> getPermutedBlocks(map<int, vector<int64_t>> regionsJoined);


vector<int64_t> genBloomFilter(vector<vector<int64_t>> transfBlock);
vector<vector<int64_t>> xorWithColKey(vector<vector<int64_t>> block, vector<int64_t> key);
map<int, vector<int64_t>> gen1stBFbasedTemplate(map<int, vector<vector<int64_t>>> blocks, map<int, vector<int64_t>> keys);
map<int, vector<int64_t>> gen2ndBFbasedTemplate(map<int, vector<int64_t>> sampleLines, map<int, vector<int>> permutations);
map<int, vector<int64_t>> genBFbasedTemplate(map<int, vector<vector<int64_t>>> transfBlocks);






double wHDBFunderHE(CryptoContext<DCRTPoly> cryptoContext,
                        LPKeyPair<DCRTPoly> keyPair,
                        vector<Ciphertext<DCRTPoly>> cipherBF1, 
                        vector<Ciphertext<DCRTPoly>> cipherBF2,
                        Ciphertext<DCRTPoly> ctTwo);

double wHDBFunderHE(CryptoContext<DCRTPoly> cryptoContext,
                        LPKeyPair<DCRTPoly> keyPair,
                        vector<Ciphertext<DCRTPoly>> cipherBF1, 
                        vector<Ciphertext<DCRTPoly>> cipherBF2);

double wHDPackedBFsunderHE(CryptoContext<DCRTPoly> cryptoContext,
                        LPKeyPair<DCRTPoly> keyPair,
                        vector<Ciphertext<DCRTPoly>> cipherBF1, 
                        vector<Ciphertext<DCRTPoly>> cipherBF2,
                        vector<int32_t> indexList);

double shiftedHDunderHE(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair, 
                        vector<Ciphertext<DCRTPoly>> cipher1, vector<Ciphertext<DCRTPoly>> cipher2,
                        Ciphertext<DCRTPoly> ctTwo, int nShift);

double shiftedHDunderHE(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair, vector<Ciphertext<DCRTPoly>> cipher1, vector<Ciphertext<DCRTPoly>> cipher2, int nShift);

double shiftedPackedHDunderHE(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair, vector<Ciphertext<DCRTPoly>> cipher1, vector<Ciphertext<DCRTPoly>> cipher2, int nShift);

double wHDBFplain(map<int, vector<int64_t>> cBF1, map<int, vector<int64_t>> cBF2);


