#define PROFILE

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
const string RUNTIMEFILE = "../experiments/runtimes/USExp1stBF.csv";
const size_t NROWS = 20;
const size_t NB = 32; 
const size_t LENB = 32;
const int SUBJECTS = 224; 


int main(int argc, char *argv[]) {

    TimeVar t, runT;
    double processingTime(0.0);

    


    TIC(t);

    for (int subject1 = 1; subject1 < SUBJECTS+1; subject1++) {
        string sample1Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_1.txt";
        string sample2Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_2.txt";

        auto sample1Subj1 = readIrisCodeAsBlocks(sample1Subj1F, NB, LENB);
        auto sample2Subj1 = readIrisCodeAsBlocks(sample2Subj1F, NB, LENB);  
        auto keysSub1 = genKeys(10, 32);

        auto sample1Subj1BF = gen1stBFbasedTemplate(sample1Subj1, keysSub1);
        
        

        cout << "Plain WHD on 1st BF ..." << endl;
        TIC(runT);
        auto sample2Subj1BF = gen1stBFbasedTemplate(sample2Subj1, keysSub1);
        auto matedshiftdHD = wHDBFplain(sample1Subj1BF, sample2Subj1BF);
        writeProcTime(RUNTIMEFILE, TOC_US(runT));
        cout << "matedshiftdHD = " << matedshiftdHD << endl;       


        
        for (int subject2 = subject1+1; subject2 < SUBJECTS+1; subject2++) {
            string sample3Subj2F = DATAFOLDER + "lg_"+to_string(subject2)+"_3.txt";
            auto sample3Subj2 = readIrisCodeAsBlocks(sample3Subj2F, NB, LENB);
            auto keysSub2 = genKeys(10, 32);
            

            TIC(runT);
            auto sample3Subj2BF = gen1stBFbasedTemplate(sample3Subj2, keysSub2); 
            auto nonMatedshiftdHD = wHDBFplain(sample1Subj1BF, sample3Subj2BF);
            writeProcTime(RUNTIMEFILE, TOC_US(runT)); 
            cout << "nonMatedshiftdHD = " << nonMatedshiftdHD << endl;

        }
        
    }

    processingTime = TOC(t);
    cout << "Time : "   << processingTime << "ms" << endl;
    


    return 0;
}