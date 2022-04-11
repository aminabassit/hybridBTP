#include "../include/bloom.hpp"
#include <string>
#include <fstream>
#include <iomanip>



using namespace std;



bool areEqual(vector<int64_t> v, vector<int64_t> u){
    auto len = v.size();
    for (size_t i = 0; i < len; i++)
    {
        if (v.at(i) != u.at(i))
        {
            return false;
        }        
    }    
    return true;
}

vector<int> vectIntPermutation(int seed, int begin, int end){
	vector<int> permuted;
	srand(seed + time(0));
	for(auto i = begin; i<end; i++){
		permuted.push_back(i);
	}
	shuffle(permuted.begin(), permuted.end(), default_random_engine(rand()));	
	return permuted;
}

map<int, vector<int>> generatePermutations(int seed, size_t nPerm, size_t len){
	map<int, vector<int>> permutedMap;
	for(size_t i = 0; i<nPerm; i++){
		permutedMap[i] = vectIntPermutation( seed + i, 0, len);
	}
	return permutedMap;
}


map<int, vector<int64_t>> genKeys(size_t hb, size_t nBlocks){

    map<int, vector<int64_t>> keys;
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> distBin(0,1); 
    for (size_t i = 0; i < nBlocks; i++)
    {
        for (size_t j = 0; j < hb; j++)
        {
            keys[i].push_back(distBin(rng));
        }
        
    }
    return keys;
}


int64_t bin2Int(vector<int64_t> binVect){
    int64_t dec = 0;
    int n = binVect.size();
    for(int i = n-1; i>=0 ; i--)
    {           
        dec += binVect.at(i) *(1<<(n-i-1)); 
    }   
    return dec;
}


const vector<string> explode(const string& s, const char& c)
{
	string buff{ "" };
	vector<string> v;

	for (auto n : s)
	{
		if (n != c) buff += n; else
			if (n == c && buff != "") { v.push_back(buff); buff = ""; }
	}
	if (buff != "") v.push_back(buff);

	return v;
}


// vector<vector<int> > trans_vec(b[0].size(), vector<int>());

//     for (int i = 0; i < b.size(); i++)
//     {
//         for (int j = 0; j < b[i].size(); j++)
//         {
//             trans_vec[j].push_back(b[i][j]);
//         }
//     }

vector<vector<int64_t> > transpose2DVect(vector<vector<int64_t> > v)
{
    vector<vector<int64_t> > transVec(v.at(0).size(), vector<int64_t>());

    for (size_t i = 0; i < v.size(); i++)
    {   
        for (size_t j = 0; j < v.at(i).size(); j++)
        {
            transVec.at(j).push_back(v.at(i).at(j));
        }
    }

    return transVec;
}


map<int, vector<vector<int64_t>>> readIrisCodeAsBlocks(string filename, size_t nB, size_t lenB){
    ifstream in_file; 
    string line;
    vector<string> strVect;
    map<int, vector<vector<int64_t>>> irisCode;
    string str;
    in_file.open(filename);
    if (!in_file)
    {
        throw invalid_argument("file open error");
    }
    
    for (size_t b = 0; b < nB; b++){
        irisCode[b].resize(lenB);
        for (size_t i = 0; i < lenB; i++){
            getline(in_file, line);
            strVect = explode(line, ',');
            for(auto s : strVect){
            irisCode[b].at(i).push_back( stoi(s)); 
            }
            strVect.resize(0);
        } 
    }    
    in_file.close();
    return irisCode;
}








map<int, vector<int64_t>> readIrisCodeAsLines(string filename){
    ifstream in_file; 
    string line;
    vector<string> strVect;
    map<int, vector<int64_t>> irisCode;
    string str;
    in_file.open(filename);
    if (!in_file)
    {
        throw invalid_argument("file open error");
    }
    
    for (size_t b = 0; b < 20; b++){
        irisCode[b].resize(0);
        getline(in_file, line);
        strVect = explode(line, ',');
        for(auto s : strVect){
            irisCode[b].push_back(stoi(s)); 
        }
        strVect.resize(0);
    }    
    in_file.close();
    return irisCode;
}



void writeProcTime(string filename, double processingTime){
    ofstream out_file;
    out_file.open(filename, ios::app);
    if (!out_file)
    {
        throw "Error creating file";
    }
    out_file << processingTime << endl;  
    out_file.close();
}








/////////////////////////////////////////////


map<int, vector<vector<int64_t>>> splitRegions( map<int, vector<int64_t>> sampleLines){

    map<int, vector<vector<int64_t>>> regions;
    size_t halfLine = sampleLines[0].size()/2;


    for (size_t i = 0; i < 10; i++)
    {
        
        vector<int64_t> r0(sampleLines.at(i).begin(), sampleLines.at(i).begin()+halfLine);
        vector<int64_t> r1(sampleLines.at(i).begin()+halfLine, sampleLines.at(i).end());
        vector<int64_t> r2(sampleLines.at(i+10).begin(), sampleLines.at(i+10).begin()+halfLine);
        vector<int64_t> r3(sampleLines.at(i+10).begin()+halfLine, sampleLines.at(i+10).end()); 

        for (size_t j = 0; j < halfLine; j+=32)
        {
            vector<int64_t> r00(r0.begin()+i, r0.begin()+i+32);
            vector<int64_t> r11(r1.begin()+i, r1.begin()+i+32);
            vector<int64_t> r22(r2.begin()+i, r2.begin()+i+32);
            vector<int64_t> r33(r3.begin()+i, r3.begin()+i+32);
            regions[0].push_back(r00);
            regions[1].push_back(r11);
            regions[2].push_back(r22);
            regions[3].push_back(r33);
        }


    }
    return regions;
}





map<int, vector<vector<int64_t>>> permuteRegions(map<int, vector<vector<int64_t>>> regions, map<int, vector<int>> permutations){
    map<int, vector<vector<int64_t>>> regionsPerm;
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < regions[i].size(); j++)
        {
            regionsPerm[i].push_back(regions[i].at(permutations[i].at(j)));
        }
        
    }
    return regionsPerm;
}

map<int, vector<int64_t>> joinRegions( map<int, vector<vector<int64_t>>> regions){
    map<int, vector<int64_t>> regionsJoined;
    size_t l = 0;
    for (size_t r = 0; r < regions.size(); r++)
    {
        for (size_t i = 0; i < regions[r].size(); i+=8)
        {   
            vector<int64_t> line;
            for (size_t j = 0; j < 8; j++)
            {
                line.insert(line.end(), regions[r].at(i+j).begin(), regions[r].at(i+j).end());
            }
            regionsJoined[l]= line;
            l++;
        }        
    }    
    return regionsJoined;
}


map<int, vector<vector<int64_t>>> getPermutedBlocks(map<int, vector<int64_t>> regionsJoined){

    map<int, vector<vector<int64_t>>> blocks;
    size_t b = 0;

    
    for (size_t i = 0; i < regionsJoined.size(); i+=10)
    {   
        vector<vector<int64_t>> blockSet;
        for (size_t j = 0; j < 10; j++)
        {
            blockSet.push_back(regionsJoined[i+j]);
        }
        auto blockSetTransp = transpose2DVect(blockSet);

        
        for (size_t j = 0; j < blockSetTransp.size(); j+=32)
        {
            vector<vector<int64_t>> block(blockSetTransp.begin()+j, blockSetTransp.begin()+j+32);
            blocks[b] = block;
            b++;
        }
        
        
    } 

    return blocks;
}




vector<int64_t> genBloomFilter(vector<vector<int64_t>> transfBlock){
    auto nWords = transfBlock.size();
    auto nBits = transfBlock.at(0).size();
    
    vector<int64_t> bf((int64_t)pow(2,nBits), 0);
    vector<int64_t> indeX;
    // #pragma omp parallel for schedule(static,1)
    for (size_t i = 0; i < nWords; i++)
    {
        auto index = bin2Int(transfBlock.at(i));
        bf.at(index) = (int64_t) 1;
        indeX.push_back(index);
    } 
    // cout <<  "indeX = " << indeX << endl;
    return bf;
}



map<int, vector<int64_t>> genBFbasedTemplate(map<int, vector<vector<int64_t>>> transfBlocks){
    map<int, vector<int64_t>> bfs;
    size_t nBlocks = transfBlocks.size();
    #pragma omp parallel for schedule(static,1)
    for (size_t i = 0; i < nBlocks; i++)
    {
        bfs[i]= genBloomFilter(transfBlocks[i]);
    }
    return bfs;
}

vector<vector<int64_t>> xorWithColKey(vector<vector<int64_t>> block, vector<int64_t> key){
    vector<vector<int64_t>> transfBlock;
    size_t nWords = block.size();
    size_t nBits = block.at(0).size();

    for (size_t i = 0; i < nWords; i++)
    {
        vector<int64_t> xoredBlock;
        for (size_t j = 0; j < nBits; j++)
        {
            xoredBlock.push_back(block.at(i).at(j) ^ key.at(j));
        }        
        transfBlock.push_back(xoredBlock);
    }    
    return transfBlock;
}

map<int, vector<int64_t>> gen1stBFbasedTemplate(map<int, vector<vector<int64_t>>> blocks, map<int, vector<int64_t>> keys){
    map<int, vector<int64_t>> bfs;
    size_t nBlocks = blocks.size();
    for (size_t i = 0; i < nBlocks; i++)
    {
        auto transfBlock = xorWithColKey(blocks[i], keys[i]);
        bfs[i]= genBloomFilter(transfBlock);
    }
    return bfs;
}



map<int, vector<int64_t>> gen2ndBFbasedTemplate(map<int, vector<int64_t>> sampleLines, map<int, vector<int>> permutations){
    map<int, vector<int64_t>> bfs;
    auto regions = splitRegions(sampleLines);
    auto regionsPerm = permuteRegions(regions, permutations);
    auto joinedRe = joinRegions(regionsPerm);
    auto blocks = getPermutedBlocks(joinedRe);
    size_t nBlocks = blocks.size();
    for (size_t i = 0; i < nBlocks; i++)
    {        
        bfs[i]= genBloomFilter(blocks[i]);
    }
    return bfs;
}



///////////////////////////////////////
//                                   //
//   Hybrid BTP: 1st BF under HE     //
//                                   //
//////////////////////////////////////






double wHDBFunderHE(CryptoContext<DCRTPoly> cryptoContext,
                        LPKeyPair<DCRTPoly> keyPair,
                        vector<Ciphertext<DCRTPoly>> cipherBF1, 
                        vector<Ciphertext<DCRTPoly>> cipherBF2,
                        Ciphertext<DCRTPoly> ctTwo){
    
    
    size_t nBlocks = 32;
    usint lenBF = 1024;

    vector<Ciphertext<DCRTPoly>> ciphBF1XorBF2, ciphBF1SumBF2;
    ciphBF1XorBF2.resize(nBlocks); 
    ciphBF1SumBF2.resize(nBlocks); 
    


    #pragma omp parallel for schedule(static,1)
    for (size_t i = 0; i < nBlocks; i++)
    {
        auto ciphBF1PlsBF2 = cryptoContext->EvalAdd(cipherBF1.at(i), cipherBF2.at(i));
        auto ciph2multBF1multBF2 = cryptoContext->EvalMultMany({ctTwo, cipherBF1.at(i), cipherBF2.at(i)});
        auto ciphBF1XorBF2ct = cryptoContext->EvalSub(ciphBF1PlsBF2, ciph2multBF1multBF2);
        ciphBF1XorBF2.at(i) = cryptoContext->EvalSum(ciphBF1XorBF2ct, lenBF); // numerator
        ciphBF1SumBF2.at(i) = cryptoContext->EvalSum(ciphBF1PlsBF2, lenBF); // denominator
    }
    vector<double> dhVect(32, 0.0);
    // dhVect.resize(32);

    #pragma omp parallel for schedule(static,1)    
    for (size_t i = 0; i < nBlocks; i++)
    {
        Plaintext ptBF1XorBF2, ptOnesBF1andBF2;

        cryptoContext->Decrypt(keyPair.secretKey, ciphBF1XorBF2.at(i), &ptBF1XorBF2);
        cryptoContext->Decrypt(keyPair.secretKey, ciphBF1SumBF2.at(i), &ptOnesBF1andBF2);
        ptBF1XorBF2->SetLength(1);
        ptOnesBF1andBF2->SetLength(1);
        auto numerator = ptBF1XorBF2->GetPackedValue().at(0);
        auto denominator = ptOnesBF1andBF2->GetPackedValue().at(0);        
        dhVect.at(i) = (double) numerator / denominator;
    }

    
    double wHD = accumulate(dhVect.begin(), dhVect.end(), double(0.0)) / 32;

    return wHD;
}


double wHDBFunderHE(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair,
                        vector<Ciphertext<DCRTPoly>> cipherBF1, vector<Ciphertext<DCRTPoly>> cipherBF2){
    
    
    size_t nBlocks = 32;
    usint lenBF = 1024;

    vector<Ciphertext<DCRTPoly>> ciphBF1XorBF2, ciphBF1SumBF2;
    ciphBF1XorBF2.resize(nBlocks); 
    ciphBF1SumBF2.resize(nBlocks); 
    


    #pragma omp parallel for schedule(static,1)
    for (size_t i = 0; i < nBlocks; i++)
    {
        auto ciphBF1PlsBF2 = cryptoContext->EvalAdd(cipherBF1.at(i), cipherBF2.at(i));
        auto ciphBF1multBF2 = cryptoContext->EvalMult(cipherBF1.at(i), cipherBF2.at(i));
        auto ciph2multBF1multBF2 = cryptoContext->EvalAdd(ciphBF1multBF2, ciphBF1multBF2);
        auto ciphBF1XorBF2ct = cryptoContext->EvalSub(ciphBF1PlsBF2, ciph2multBF1multBF2);
        ciphBF1XorBF2.at(i) = cryptoContext->EvalSum(ciphBF1XorBF2ct, lenBF); // numerator
        ciphBF1SumBF2.at(i) = cryptoContext->EvalSum(ciphBF1PlsBF2, lenBF); // denominator
    }
    vector<double> dhVect(32, 0.0);
    // dhVect.resize(32);

    #pragma omp parallel for schedule(static,1)    
    for (size_t i = 0; i < nBlocks; i++)
    {
        Plaintext ptBF1XorBF2, ptOnesBF1andBF2;

        cryptoContext->Decrypt(keyPair.secretKey, ciphBF1XorBF2.at(i), &ptBF1XorBF2);
        cryptoContext->Decrypt(keyPair.secretKey, ciphBF1SumBF2.at(i), &ptOnesBF1andBF2);
        ptBF1XorBF2->SetLength(1);
        ptOnesBF1andBF2->SetLength(1);
        auto numerator = ptBF1XorBF2->GetPackedValue().at(0);
        auto denominator = ptOnesBF1andBF2->GetPackedValue().at(0);        
        dhVect.at(i) = (double) numerator / denominator;
    }

    
    double wHD = accumulate(dhVect.begin(), dhVect.end(), double(0.0)) / 32;

    return wHD;
}


double wHDPackedBFsunderHE(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair,
                        vector<Ciphertext<DCRTPoly>> cipherBF1, vector<Ciphertext<DCRTPoly>> cipherBF2,
                        vector<int32_t> indexList){
    
    size_t nBlocks = 32;
    usint lenBF = 1024;
    size_t nCiphers = cipherBF1.size();

    

    vector<Ciphertext<DCRTPoly>> ciphBF1XorBF2, ciphBF1SumBF2;
    ciphBF1XorBF2.resize(nBlocks); 
    ciphBF1SumBF2.resize(nBlocks); 
    size_t nPackedBFsPerCipher = (size_t) nBlocks/nCiphers;

    // cout << "nCiphers = " << nCiphers << endl;
    // cout << "nPackedBFsPerCipher = " << nPackedBFsPerCipher << endl;


    #pragma omp parallel for schedule(static,1)
    for (size_t i = 0; i < nCiphers; i++)
    {
        auto ciphBF1PlsBF2 = cryptoContext->EvalAdd(cipherBF1.at(i), cipherBF2.at(i));
        auto ciphBF1multBF2 = cryptoContext->EvalMult(cipherBF1.at(i), cipherBF2.at(i));
        auto ciph2multBF1multBF2 = cryptoContext->EvalAdd(ciphBF1multBF2, ciphBF1multBF2);
        auto ciphBF1XorBF2ct = cryptoContext->EvalSub(ciphBF1PlsBF2, ciph2multBF1multBF2);
        auto ciphBF1XorBF2Packed = cryptoContext->EvalSum(ciphBF1XorBF2ct, lenBF); // numerator
        auto ciphBF1SumBF2Packed = cryptoContext->EvalSum(ciphBF1PlsBF2, lenBF); // denominator

        size_t index = i * nPackedBFsPerCipher;

        ciphBF1XorBF2.at(index) = ciphBF1XorBF2Packed;
        ciphBF1SumBF2.at(index) = ciphBF1SumBF2Packed;        
        
        // exctract packed ciphers 
        for (size_t j = 0 ; j < nPackedBFsPerCipher-1; j++)
        {
            index += 1;
            ciphBF1XorBF2.at(index) = cryptoContext->EvalAtIndex(ciphBF1XorBF2Packed, indexList[j]);
            ciphBF1SumBF2.at(index) = cryptoContext->EvalAtIndex(ciphBF1SumBF2Packed, indexList[j]);
        }
        
    }


    vector<double> dhVect(nBlocks, 0.0);

    #pragma omp parallel for schedule(static,1)    
    for (size_t i = 0; i < nBlocks; i++)
    {
        Plaintext ptBF1XorBF2, ptOnesBF1andBF2;

        cryptoContext->Decrypt(keyPair.secretKey, ciphBF1XorBF2.at(i), &ptBF1XorBF2);
        cryptoContext->Decrypt(keyPair.secretKey, ciphBF1SumBF2.at(i), &ptOnesBF1andBF2);
        ptBF1XorBF2->SetLength(1);
        ptOnesBF1andBF2->SetLength(1);
        auto numerator = ptBF1XorBF2->GetPackedValue().at(0);
        auto denominator = ptOnesBF1andBF2->GetPackedValue().at(0);    
        // cout << " numerator = " << numerator << " denominator = " << denominator << endl; 
        dhVect.at(i) = (double) numerator / denominator;
    }

    
    double wHD = (double) accumulate(dhVect.begin(), dhVect.end(), double(0.0)) / nBlocks;
    // dhVect.resize(nBlocks);

    return wHD;
}


double shiftedHDunderHE(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair, vector<Ciphertext<DCRTPoly>> cipher1, vector<Ciphertext<DCRTPoly>> cipher2, Ciphertext<DCRTPoly> ctTwo, int nShift){

    size_t nRows = 20;
    usint lenRow = 512;

    double minShiftHD(10000.0); 

    
    for (int j = -nShift; j < nShift+1; j++){        

        double shHD(0.0);
        #pragma omp parallel for schedule(static,1)
        for (size_t i = 0; i < nRows; i++){  

            auto shiftedCipher2 = (j == 0)? cipher2.at(i): cryptoContext->EvalAtIndex(cipher2.at(i), j);            
            

            auto ciph1PlsCiph2 = cryptoContext->EvalAdd(cipher1.at(i), shiftedCipher2);
            auto ciph2multCiph1multCiph2 = cryptoContext->EvalMultMany({ctTwo, cipher1.at(i), shiftedCipher2});
            auto ciph1XorCiph2 = cryptoContext->EvalAdd(ciph1PlsCiph2, ciph2multCiph1multCiph2);
            auto hdC1C2 = cryptoContext->EvalSum(ciph1XorCiph2, lenRow); 

            Plaintext plainHD;

            cryptoContext->Decrypt(keyPair.secretKey, hdC1C2, &plainHD);
            plainHD->SetLength(1);

            shHD += plainHD->GetPackedValue().at(0);       
        }

        shHD = shHD/(double)(lenRow*nRows);

        if(shHD <= minShiftHD){
            minShiftHD = shHD;
        }
        
    }

    return minShiftHD;
}



double shiftedHDunderHE(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair, vector<Ciphertext<DCRTPoly>> cipher1, vector<Ciphertext<DCRTPoly>> cipher2, int nShift){

    size_t nRows = 20;
    usint lenRow = 512;

    double minShiftHD(10000.0); 

    
    for (int j = -nShift; j < nShift+1; j++){        

        double shHD(0.0);
        // #pragma omp parallel for schedule(static,1)        
        for (size_t i = 0; i < nRows; i++){  

            auto shiftedCipher2 = (j == 0)? cipher2.at(i): cryptoContext->EvalAtIndex(cipher2.at(i), j);            
            

            auto ciph1PlsCiph2 = cryptoContext->EvalAdd(cipher1.at(i), shiftedCipher2);
            auto ciphCiph1multCiph2 = cryptoContext->EvalMult(cipher1.at(i), shiftedCipher2);
            auto ciph2multCiph1multCiph2 = cryptoContext->EvalAdd(ciphCiph1multCiph2, ciphCiph1multCiph2);
            auto ciph1XorCiph2 = cryptoContext->EvalSub(ciph1PlsCiph2, ciph2multCiph1multCiph2);
            auto hdC1C2 = cryptoContext->EvalSum(ciph1XorCiph2, lenRow); 

            Plaintext plainHD;

            cryptoContext->Decrypt(keyPair.secretKey, hdC1C2, &plainHD);
            plainHD->SetLength(1);

            shHD += plainHD->GetPackedValue().at(0);       
        }

        shHD = shHD/(double)(lenRow*nRows);

        if(shHD <= minShiftHD){
            minShiftHD = shHD;
        }
        
    }

    return minShiftHD;
}

double shiftedPackedHDunderHE(CryptoContext<DCRTPoly> cryptoContext, LPKeyPair<DCRTPoly> keyPair, vector<Ciphertext<DCRTPoly>> cipher1, vector<Ciphertext<DCRTPoly>> cipher2, int nShift){

    size_t nRows = 20;
    usint lenRow = 512;
    size_t nCiphers = cipher1.size();
    
    auto ringDim = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    size_t nLastPackedRows = (size_t) nRows%(ringDim/lenRow);
    size_t lenLastPackedRows = lenRow * nLastPackedRows;

    // cout << "nCiphers = " << nCiphers << endl;
    // cout << "nLastPackedRows = " << nLastPackedRows << endl;
    // cout << "lenLastPackedRows = " << lenLastPackedRows << endl;

    double minShiftHD(1000.0); 

    #pragma omp parallel for schedule(static,1)
    for (int j = -nShift; j < nShift+1; j++){        

        double shHD(0.0);       
        for (size_t i = 0; i < nCiphers; i++){  

            auto shiftedCipher2 = (j == 0)? cipher2.at(i): cryptoContext->EvalAtIndex(cipher2.at(i), j);            
            

            auto ciph1PlsCiph2 = cryptoContext->EvalAdd(cipher1.at(i), shiftedCipher2);
            auto ciphCiph1multCiph2 = cryptoContext->EvalMult(cipher1.at(i), shiftedCipher2);
            auto ciph2multCiph1multCiph2 = cryptoContext->EvalAdd(ciphCiph1multCiph2, ciphCiph1multCiph2);
            auto ciph1XorCiph2 = cryptoContext->EvalSub(ciph1PlsCiph2, ciph2multCiph1multCiph2);
            size_t batchSize = (i == nCiphers-1 && i != 0)? lenLastPackedRows:ringDim; 
            // cout << "batchSize = " << batchSize << endl;
            auto hdC1C2 = cryptoContext->EvalSum(ciph1XorCiph2, batchSize); 

            Plaintext plainHD;

            cryptoContext->Decrypt(keyPair.secretKey, hdC1C2, &plainHD);
            plainHD->SetLength(1);
            // cout << "plainHD = " << plainHD << endl;
            shHD += plainHD->GetPackedValue().at(0);       
        }

        shHD = shHD/(double)(lenRow*nRows);

        if(shHD < minShiftHD){
            minShiftHD = shHD;
        }
        
    }

    return minShiftHD;
}


double wHDBFplain(map<int, vector<int64_t>> cBF1, map<int, vector<int64_t>> cBF2){
    double wHD(0.0);
    size_t nBlocks = 32;
    usint lenBF = 1024;
    
    for (size_t i = 0; i < nBlocks; i++){

        vector<int64_t> xorBF1BF2(lenBF,0);

        transform(cBF1[i].begin(), cBF1[i].end(), cBF2[i].begin(), xorBF1BF2.begin(), bit_xor<int>());
        auto numerator = accumulate(xorBF1BF2.begin(), xorBF1BF2.end(), double(0.0));

        auto sumBF1 = accumulate(cBF1[i].begin(), cBF1[i].end(), double(0.0));
        auto sumBF2 = accumulate(cBF2[i].begin(), cBF2[i].end(), double(0.0));
        auto denominator = sumBF1 + sumBF2;

        wHD += (double) numerator / denominator;
    }

    return wHD / nBlocks;
}







// double wHDBFplain(map<int, vector<int>> cBF1, map<int, vector<int>> cBF2){
//     double wHD(0.0);
//     size_t nBlocks = 32;
//     usint lenBF = 1024;
//     vector<int> vectTow(1024, -2);

    
//     for (size_t i = 0; i < nBlocks; i++){
        
//         auto sumBF1 = accumulate(cBF1[i].begin(), cBF1[i].end(), double(0.0));
//         auto sumBF2 = accumulate(cBF2[i].begin(), cBF2[i].end(), double(0.0));
//         auto denominator = sumBF1 + sumBF2;




//         vector<int> addBF1B2(lenBF,0), multBF1B2(lenBF,0), mult2BF1B2(lenBF,0);
//         transform(cBF1[i].begin(), cBF1[i].end(), cBF2[i].begin(), addBF1B2.begin(), plus<int>());

//         auto denominator = accumulate(addBF1B2.begin(), addBF1B2.end(), double(0.0));

//         transform(cBF1[i].begin(), cBF1[i].end(), cBF2[i].begin(), multBF1B2.begin(), multiplies<int>());
//         transform(multBF1B2.begin(), multBF1B2.end(), vectTow.begin(), mult2BF1B2.begin(), multiplies<int>());



//         auto numerator = accumulate(mult2BF1B2.begin(), mult2BF1B2.end(), double(0.0));
//         wDH += (double) numerator / denominator ;
//     }



//     return wDH / nBlocks;
// }












