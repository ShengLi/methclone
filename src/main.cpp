//
//  main.cpp
//  methclone
//
//  Created by Sheng Li on 7/17/13.
//  Copyright (c) 2013 Weill Medical College of Cornell University. All rights reserved.
//

#include <iostream>
#include "bamByChr.h"

int main(int argc, const char * argv[])
{
    if (argc<5) {
        std::cerr << "Usage:\n./methclone bam1 bam2 out.gz sampleName" << std::endl;
        return 1;
    } else {
        // input files
        std::string bamFile1=std::string(argv[1]);
        std::string bamFile2=std::string(argv[2]);
        const char * outFile=argv[3];
        std::string sample=argv[4];
        //default methdiff: 0 // trying to capture all range of methylation dynamics
        int m=0;
        //default distance cutoff: 72 // trying to match seq length
        int d=72;
        //default read coverage: 60 
        int cov=60;
        
        std::cerr << "Default program paramters for methclone: \nMax loci width: " << d << " bp" << std::endl;
        std::cerr << "Min methylation difference cutoff: " << m << "%" << std::endl;
        std::cerr << "Min read coverage: " << cov << std::endl;
        // input: 2 bam files
        //std::string bamFile1="/Users/shl2018/Awork/project/AML_RRBS/eclone/eclone/test_1.bam";
        //std::string bamFile2="/Users/shl2018/Awork/project/AML_RRBS/eclone/eclone/test_2.bam";
        // ouput: entropy text file
        //const char * outFile="/Users/shl2018/Awork/project/AML_RRBS/eclone/eclone/test.entropy.gz";
        // sample name
        //std::string sample="UA1";
        
        // read-in bam by chromosome for each bam together.
        // 1. read alignment and convert to methylation patterns
        // 2. reformat the methylation format into 16 pattern vector
        // 3. compare first bam with the second bam to calculate the entropy
        // 4. output entropy.
        bamToLociMeth(bamFile1, bamFile2, outFile, sample, d, cov, m);
        return 0;
    }
}

