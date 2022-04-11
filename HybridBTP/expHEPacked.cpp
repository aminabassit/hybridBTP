#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>


#include "palisade.h"

#include "include/bloom.hpp"

using namespace std;
using namespace lbcrypto;

const string DATAFOLDER = "../data/IITD_IrisCodes_left_lines/";
const size_t NROWS = 20;
const size_t LENROW = 512;
const int SUBJECTS = 224;

int main(int argc, char *argv[]) {

    TimeVar t, runT;
    double processingTime(0.0);

    int64_t nBits(16);
    uint64_t m(32768);
    auto plaintextMod = FirstPrime<NativeInteger>(nBits, m);
    usint plaintextModulus = (usint) plaintextMod.ConvertToInt();


    std::cout << "plaintextModulus = " << plaintextModulus << std::endl;
    EncodingParams encodingParams( std::make_shared<EncodingParamsImpl>(plaintextModulus));

    double sigma = 3.2;

    


    string sBits = argv[1];
    usint mD;
    cout << "HEStd_"+sBits+"_classic" << endl;
    SecurityLevel securityLevel;

    if (sBits == "128")
    {
        securityLevel = HEStd_128_classic;
        mD = 1;
    }
    if (sBits == "192")
    {
        securityLevel = HEStd_192_classic;
        mD = 2;
    }
    if (sBits == "256")
    {
        securityLevel = HEStd_256_classic;
        mD = 2;
    } 

    string RUNTIMEFILE = "../experiments/runtimes/MSExpHE"+sBits+"Packed.csv";

    CryptoContext<DCRTPoly> cryptoContext =
    CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(
            encodingParams, securityLevel, sigma, 0, mD, 0, OPTIMIZED, 2, 0, 30, 0);



    cryptoContext->Enable(ENCRYPTION);
    cryptoContext->Enable(SHE);

    LPKeyPair<DCRTPoly> keyPair;
    keyPair = cryptoContext->KeyGen();

    cryptoContext->EvalMultKeysGen(keyPair.secretKey);
    cryptoContext->EvalSumKeyGen(keyPair.secretKey);

    // Generate the rotation evaluation keys
    vector<int> shiftIndexes = {-8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8};
    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, shiftIndexes);

    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    size_t nPackedRowsPerCipher = (size_t) ringDim/LENROW;
    size_t numCiphers = (NROWS%nPackedRowsPerCipher == 0)? (size_t) NROWS/nPackedRowsPerCipher : (size_t) NROWS/nPackedRowsPerCipher + 1;

    cout << "ringDim = " << ringDim << endl;

    
    
    
    
    
    /////////////////////////////////////


    vector<Ciphertext<DCRTPoly>> ctSamp1Suj1, ctSamp2Suj1, ctSamp3Suj2;
    


    TIC(t);

    for (int subject1 = 1; subject1 < SUBJECTS+1; subject1++) {
        string sample1Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_1.txt";
        string sample2Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_2.txt";

        cout << "sample1Subj1F = " << sample1Subj1F << endl;
        cout << "sample2Subj1F = " << sample2Subj1F << endl;

        auto sample1Subj1 = readIrisCodeAsLines(sample1Subj1F);
        auto sample2Subj1 = readIrisCodeAsLines(sample2Subj1F);  

        ctSamp1Suj1.resize(numCiphers);
        ctSamp2Suj1.resize(numCiphers);

        for (size_t i = 0; i < numCiphers; i++) {
            vector<int64_t> packedRows;
            for (size_t j = i * nPackedRowsPerCipher; j < (i+1) * nPackedRowsPerCipher; j++)
            {
                packedRows.insert(packedRows.end(), sample1Subj1[j].begin(), sample1Subj1[j].end());
            } 
            Plaintext pTSamp1Suj1 = cryptoContext->MakePackedPlaintext(packedRows);
            ctSamp1Suj1.at(i) = cryptoContext->Encrypt(keyPair.publicKey, pTSamp1Suj1);
        }         

        cout << "ctSamp1Suj1 is encrypted ..." << endl;


        TIC(runT);
        #pragma omp parallel for schedule(static,1)
        for (size_t i = 0; i < numCiphers; i++) {
            vector<int64_t> packedRows;
            for (size_t j = i * nPackedRowsPerCipher; j < (i+1) * nPackedRowsPerCipher; j++)
            {
                packedRows.insert(packedRows.end(), sample2Subj1[j].begin(), sample2Subj1[j].end());
            } 
            Plaintext pPSamp2Suj1 = cryptoContext->MakePackedPlaintext(packedRows);
            ctSamp2Suj1.at(i) = cryptoContext->Encrypt(keyPair.publicKey, pPSamp2Suj1);
        } 
        auto matedshiftdHD = shiftedPackedHDunderHE(cryptoContext, keyPair, ctSamp1Suj1, ctSamp2Suj1, 8);
        writeProcTime(RUNTIMEFILE, TOC_MS(runT));
        cout << "matedshiftdHD = " << matedshiftdHD << endl;

        


        
        for (int subject2 = subject1+1; subject2 < SUBJECTS+1; subject2++) {
            string sample3Subj2F = DATAFOLDER + "lg_"+to_string(subject2)+"_3.txt";
            cout << "sample3Subj2F = " << sample3Subj2F << endl;
            auto sample3Subj2 = readIrisCodeAsLines(sample3Subj2F);
            ctSamp3Suj2.resize(NROWS);

            TIC(runT);
            #pragma omp parallel for schedule(static,1)
            for (size_t i = 0; i < numCiphers; i++)
            {
                vector<int64_t> packedRows;
                for (size_t j = i * nPackedRowsPerCipher; j < (i+1) * nPackedRowsPerCipher; j++)
                {
                    packedRows.insert(packedRows.end(), sample3Subj2[j].begin(), sample3Subj2[j].end());
                }                 
                Plaintext pTSamp3Suj2 = cryptoContext->MakePackedPlaintext(packedRows);
                ctSamp3Suj2.at(i) = cryptoContext->Encrypt(keyPair.publicKey, pTSamp3Suj2);
            }
            // cout << "ctSamp3Suj2 is encrypted ..." << endl;            
            auto nonMatedshiftdHD = shiftedPackedHDunderHE(cryptoContext, keyPair, ctSamp1Suj1, ctSamp3Suj2, 8);
            writeProcTime(RUNTIMEFILE, TOC_MS(runT));
            cout << "nonMatedshiftdHD = " << nonMatedshiftdHD << endl;

        }
        
    }

    processingTime = TOC(t);
    cout << "Time : "   << processingTime << "ms" << endl;
    


    return 0;
}