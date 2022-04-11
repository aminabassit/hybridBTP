#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "cryptocontextgen.h"
#include "palisade.h"
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "pubkeylp-ser.h"
#include "scheme/bfvrns/bfvrns-ser.h"

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

    
    

    int64_t nBits(16);
    uint64_t m(32768);
    auto plaintextMod = FirstPrime<NativeInteger>(nBits, m);
    usint plaintextModulus = (usint) plaintextMod.ConvertToInt();
    EncodingParams encodingParams( std::make_shared<EncodingParamsImpl>(plaintextModulus));
    
    std::cout << "plaintextModulus = " << plaintextModulus << std::endl;

    double sigma(3.2);

    
    string sBits = argv[1];


    cout << "HEStd_"+sBits+"_classic" << endl;
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
 
    string btpFolder = "../experiments/templates/Hybrid"+sBits+"/";
    string btpFile("TemplateHybrid"+sBits+"-");
    CryptoContext<DCRTPoly> cryptoContext =    CryptoContextFactory<DCRTPoly>::genCryptoContextBFVrns(
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
    





    vector<Ciphertext<DCRTPoly>> ctSamp1Suj1BF;  

    int subject1 = 1;
    string sample1Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_1.txt";

    auto sample1Subj1 = readIrisCodeAsBlocks(sample1Subj1F, NB, LENB);
    auto keysSub1 = genKeys(10, 32);       
        
        
    ctSamp1Suj1BF.resize(numCiphers);
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

    

    for (size_t i = 0; i < ctSamp1Suj1BF.size(); i++)
    {
        if (!Serial::SerializeToFile(btpFolder + btpFile + to_string(i) + ".txt", ctSamp1Suj1BF.at(i), SerType::BINARY)) {
        std::cerr
            << "Error writing serialization of a hybrid template to "<< btpFile <<".txt"
            << std::endl;
        return 1;
        }   
    }
    
    cout << "Hybrid template saved ......." << endl; 

    return 0;
}