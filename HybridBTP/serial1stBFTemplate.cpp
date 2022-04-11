#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <cstdlib>


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
const string BTPFOLDER = "../experiments/templates/1stBF/";
const size_t NROWS = 20;
const size_t NB = 32; 
const size_t LENB = 32;
const int SUBJECTS = 224; 



int main(int argc, char *argv[]) {

    
    string btpFile("Template1stBF-");
    int subject1 = 1;
    string sample1Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_1.txt";

    auto sample1Subj1 = readIrisCodeAsBlocks(sample1Subj1F, NB, LENB); 
    auto keysSub1 = genKeys(10, 32);
    auto sample1Subj1BF = gen1stBFbasedTemplate(sample1Subj1, keysSub1);


    for (size_t i = 0; i < sample1Subj1BF.size(); i++)
    {
        if (!Serial::SerializeToFile(BTPFOLDER + btpFile + to_string(i) + ".txt", sample1Subj1BF[i], SerType::BINARY)) {
        std::cerr
            << "Error writing serialization of a 1st BF template to "<< btpFile <<".txt"
            << std::endl;
        return 1;
        }   
    }




    cout << "1st BF template saved ......." << endl; 



    return 0;
}









