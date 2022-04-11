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

const string DATAFOLDER = "../data/IITD_IrisCodes_left_lines/";
const string RUNTIMEFILE = "../experiments/runtimes/USExp2ndBF.csv";
const size_t NROWS = 20;
const int SUBJECTS = 224; 


int main(int argc, char *argv[]) {

    TimeVar t, runT;
    double processingTime(0.0);

    


    TIC(t);

    for (int subject1 = 1; subject1 < SUBJECTS+1; subject1++) {
        int seed = 564169; 
        string sample1Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_1.txt";
        string sample2Subj1F = DATAFOLDER + "lg_"+to_string(subject1)+"_2.txt";

        auto sample1Subj1 = readIrisCodeAsLines(sample1Subj1F);
        auto sample2Subj1 = readIrisCodeAsLines(sample2Subj1F);  
        auto permutationsSub1 = generatePermutations(seed, 4, 80);

        auto sample1Subj1BF = gen2ndBFbasedTemplate(sample1Subj1, permutationsSub1);
        
        

        cout << "Plain WHD on 2nd BF ..." << endl;
        TIC(runT);
        auto sample2Subj1BF = gen2ndBFbasedTemplate(sample2Subj1, permutationsSub1);
        auto matedshiftdHD = wHDBFplain(sample1Subj1BF, sample2Subj1BF);
        writeProcTime(RUNTIMEFILE, TOC_US(runT));
        cout << "matedshiftdHD = " << matedshiftdHD << endl;       


        
        for (int subject2 = subject1+1; subject2 < SUBJECTS+1; subject2++) {
            string sample3Subj2F = DATAFOLDER + "lg_"+to_string(subject2)+"_3.txt";
            auto sample3Subj2 = readIrisCodeAsLines(sample3Subj2F);
            auto permutationsSub2 = generatePermutations(seed, 4, 80);
            

            TIC(runT);
            auto sample3Subj2BF = gen2ndBFbasedTemplate(sample3Subj2, permutationsSub2); 
            auto nonMatedshiftdHD = wHDBFplain(sample1Subj1BF, sample3Subj2BF);
            writeProcTime(RUNTIMEFILE, TOC_US(runT)); 
            cout << "nonMatedshiftdHD = " << nonMatedshiftdHD << endl;

        }
        
    }

    processingTime = TOC(t);
    cout << "Time : "   << processingTime << "ms" << endl;
    


    return 0;
}