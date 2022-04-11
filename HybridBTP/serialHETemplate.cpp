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

const string DATAFOLDER = "../data/IITD_IrisCodes_left_lines/";
const size_t NROWS = 20;
const size_t LENROW = 512;
const int SUBJECTS = 224;


int main(int argc, char *argv[]) {

    
    

    int64_t nBits(16);
    uint64_t m(32768);
    auto plaintextMod = FirstPrime<NativeInteger>(nBits, m);
    usint plaintextModulus = (usint) plaintextMod.ConvertToInt();
    EncodingParams encodingParams( std::make_shared<EncodingParamsImpl>(plaintextModulus));
    std::cout << "plaintextModulus = " << plaintextModulus << std::endl;

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



    string btpFolder = "../experiments/templates/HE"+sBits+"/";
    string btpFile("TemplateHybrid"+sBits+"-");
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

    
    
    
    


    vector<Ciphertext<DCRTPoly>> ctSamp1Suj1;
    int subject1 = 1;
    string sample1Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_1.txt";
    auto sample1Subj1 = readIrisCodeAsLines(sample1Subj1F); 

    ctSamp1Suj1.resize(numCiphers);
    for (size_t i = 0; i < numCiphers; i++) {
        vector<int64_t> packedRows;
        for (size_t j = i * nPackedRowsPerCipher; j < (i+1) * nPackedRowsPerCipher; j++)
        {
            packedRows.insert(packedRows.end(), sample1Subj1[j].begin(), sample1Subj1[j].end());
        } 
        Plaintext pTSamp1Suj1 = cryptoContext->MakePackedPlaintext(packedRows);
        ctSamp1Suj1.at(i) = cryptoContext->Encrypt(keyPair.publicKey, pTSamp1Suj1);
    }


    for (size_t i = 0; i < ctSamp1Suj1.size(); i++)
    {
        if (!Serial::SerializeToFile(btpFolder + btpFile + to_string(i) + ".txt", ctSamp1Suj1.at(i), SerType::BINARY)) {
        std::cerr
            << "Error writing serialization of an HE template to "<< btpFile <<".txt"
            << std::endl;
        return 1;
        }   
    }
        
    cout << "HE template saved ......." << endl; 

    return 0;
}