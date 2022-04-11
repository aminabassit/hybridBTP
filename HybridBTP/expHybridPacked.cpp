#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "cryptocontextgen.h"
#include "palisade.h"

#include "include/bloom.hpp"

using namespace std;
using namespace lbcrypto;

const string DATAFOLDER = "../data/IITD_IrisCodes_left_blocks/";
const size_t LENBF = 1024;
const size_t NBFS = 32;
const size_t NB = 32; 
const size_t LENB = 32;
const int SUBJECTS = 224;

int main(int argc, char *argv[]) {
    


    
    TimeVar t, runT;
    double processingTime(0.0);

    int64_t nBits(16);
    uint64_t m(32768);
    auto plaintextMod = FirstPrime<NativeInteger>(nBits, m);
    usint plaintextModulus = (usint) plaintextMod.ConvertToInt();
    EncodingParams encodingParams( std::make_shared<EncodingParamsImpl>(plaintextModulus));
    std::cout << "plaintextModulus = " << plaintextModulus << std::endl;

    double sigma(3.2);
    string sBits = argv[1];
    cout << "HEStd_"+sBits+"_classic" << endl;
    string RUNTIMEFILE = "../experiments/runtimes/MSExpHybrid"+sBits+"Packed.csv";
    // string TEMPENCRUNTIMEFILE = "../experiments/runtimes/MSExpHybrid-Template"+sBits+"Packed.csv";
    // string PROBENCRUNTIMEFILE = "../experiments/runtimes/MSExpHybrid-Probe"+sBits+"Packed.csv";

    SecurityLevel securityLevel;
    if (sBits == "128")
    {
        securityLevel = HEStd_128_classic;
    }
    if (sBits == "192")
    {
        securityLevel = HEStd_192_classic;
    }
    if (sBits == "256")
    {
        securityLevel = HEStd_256_classic;
    }  


    CryptoContext<DCRTPoly> cryptoContext =
    CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(
            encodingParams, securityLevel, sigma, 0, 1, 0, OPTIMIZED, 2, 0, 30, 0);



    cryptoContext->Enable(ENCRYPTION);
    cryptoContext->Enable(SHE);


    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;

    cout << "p = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus()  << endl;
    cout << "n = " << ringDim << endl;
    cout << "log2 q = " << log2(cryptoContext->GetCryptoParameters()
                            ->GetElementParams()
                            ->GetModulus()
                            .ConvertToDouble()) << endl;

    LPKeyPair<DCRTPoly> keyPair;
    keyPair = cryptoContext->KeyGen();

    cryptoContext->EvalMultKeysGen(keyPair.secretKey);
    cryptoContext->EvalSumKeyGen(keyPair.secretKey);

    vector<int32_t> indexList;
    for (size_t i = LENBF; i < ringDim; i+=LENBF)
    {
        indexList.push_back(i);
    }

    cryptoContext->EvalAtIndexKeyGen(keyPair.secretKey, indexList);

    size_t nPackedBFsPerCipher = (size_t) ringDim/LENBF;
    size_t numCiphers = (size_t) NBFS/nPackedBFsPerCipher;
    


    /////////////////////////////////////






    vector<Ciphertext<DCRTPoly>> ctSamp1Suj1BF, ctSamp2Suj1BF, ctSamp3Suj2BF;
    
    

    TIC(t);

    for (int subject1 = 1; subject1 < SUBJECTS+1; subject1++) {

        string sample1Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_1.txt";
        string sample2Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_2.txt";

        auto sample1Subj1 = readIrisCodeAsBlocks(sample1Subj1F, NB, LENB);
        auto sample2Subj1 = readIrisCodeAsBlocks(sample2Subj1F, NB, LENB); 
        auto keysSub1 = genKeys(10, 32);

        
        
        
        ctSamp1Suj1BF.resize(numCiphers);
        ctSamp2Suj1BF.resize(numCiphers);

        
        
        // TIC(runT);
        auto sample1Subj1BF = gen1stBFbasedTemplate(sample1Subj1, keysSub1);        
        #pragma omp parallel for schedule(static,1)
        for (size_t i = 0; i < numCiphers; i++) {
            vector<int64_t> packedBFs;
            for (size_t j = i * nPackedBFsPerCipher; j < (i+1) * nPackedBFsPerCipher; j++)
            {
                packedBFs.insert(packedBFs.end(), sample1Subj1BF[j].begin(), sample1Subj1BF[j].end());
            }            
            auto pTSamp1Suj1BF = cryptoContext->MakePackedPlaintext(packedBFs);
            ctSamp1Suj1BF.at(i) = cryptoContext->Encrypt(keyPair.publicKey, pTSamp1Suj1BF);
        }         
        // writeProcTime(TEMPENCRUNTIMEFILE, TOC_MS(runT)); 

        


        TIC(runT);
        auto sample2Subj1BF = gen1stBFbasedTemplate(sample2Subj1, keysSub1);
        #pragma omp parallel for schedule(static,1)
        for (size_t i = 0; i < numCiphers; i++) {
            vector<int64_t> packedBFs;
            for (size_t j = i * nPackedBFsPerCipher; j < (i+1) * nPackedBFsPerCipher; j++)
            {
                packedBFs.insert(packedBFs.end(), sample2Subj1BF[j].begin(), sample2Subj1BF[j].end());
            } 
            auto pTSamp2Suj1BF = cryptoContext->MakePackedPlaintext(packedBFs); 
            ctSamp2Suj1BF.at(i) = cryptoContext->Encrypt(keyPair.publicKey, pTSamp2Suj1BF);
        }        
        // writeProcTime(PROBENCRUNTIMEFILE, TOC_MS(runT));  
        auto matedWHD =  wHDPackedBFsunderHE(cryptoContext, keyPair, ctSamp1Suj1BF, ctSamp2Suj1BF, indexList);  
        writeProcTime(RUNTIMEFILE, TOC_MS(runT)); 
        cout << "matedWHD = " << matedWHD << endl;        
        

        
        for (int subject2 = subject1+1; subject2 < SUBJECTS+1; subject2++) {
            

            string sample3Subj2F = DATAFOLDER + "lg_"+to_string(subject2)+"_3.txt";
            auto sample3Subj2 = readIrisCodeAsBlocks(sample3Subj2F, NB, LENB);
            auto keysSub2 = genKeys(10, 32);
                      


            ctSamp3Suj2BF.resize(numCiphers);

            

            TIC(runT);
            auto sample3Subj2BF = gen1stBFbasedTemplate(sample3Subj2, keysSub2);  
            #pragma omp parallel for schedule(static,1)
            for (size_t i = 0; i < numCiphers; i++)
            {
                vector<int64_t> packedBFs;
                for (size_t j = i * nPackedBFsPerCipher; j < (i+1) * nPackedBFsPerCipher; j++)
                {
                    packedBFs.insert(packedBFs.end(), sample3Subj2BF[j].begin(), sample3Subj2BF[j].end());
                } 
                auto pTSamp3Suj2BF = cryptoContext->MakePackedPlaintext(packedBFs); 
                ctSamp3Suj2BF.at(i) = cryptoContext->Encrypt(keyPair.publicKey, pTSamp3Suj2BF);
            }
            // writeProcTime(PROBENCRUNTIMEFILE, TOC_MS(runT)); 
            auto nonMatedWHD = wHDPackedBFsunderHE(cryptoContext, keyPair, ctSamp1Suj1BF, ctSamp3Suj2BF, indexList);  
            writeProcTime(RUNTIMEFILE, TOC_MS(runT));           
            cout << "nonMatedWHD = " << nonMatedWHD << endl;
            

        }
        
    }

    processingTime = TOC(t);
    cout << "Time : "   << processingTime << "ms" << endl;
    


    return 0;
}