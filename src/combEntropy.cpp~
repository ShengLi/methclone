//
//  combEntropy.cpp
//  methclone
//
//  Created by Sheng Li on 7/17/13.
//  Copyright (c) 2013 Weill Medical College of Cornell University. All rights reserved.
//

#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <gzstream.h>
// all 16 patterns
int getAllPatterns(std::map <std::string, int> & patterns, std::map <std::string, int> & meth)
{
    for (int i=0; i<2; i++) {
        for (int k=0; k<2; k++) {
            for (int j=0; j<2; j++) {
                for (int n=0; n<2; n++) {
                    std::ostringstream x;
                    x << i << k <<j<<n;
                    patterns[x.str()]=0;
                    meth[x.str()]=i+k+j+n;
                    //std::cout << x.str() << ":" << meth[x.str()]<< std::endl;
                }
            }
        }
    }
    return 0;
}

// take pattern vector and return all 16 patterns' frequency
int methPatterns(std::map <std::string, int> & allp, 
                 std::vector<std::string> p, 
                 std::map <std::string, int> patternMeth, 
                 double & meth,
                 std::string & maxp)
{
    int freq=p.size();
    //count the freqeuncy of different patterns.
    for (int i=0; i<freq; i++) {
        allp[p[i]]++;
    }
    int count=0;
    int f=0;
    int maxf=0;
    for (std::map <std::string, int>::const_iterator it=allp.begin(); it!=allp.end(); it++) {
        f=it->second;
        count += f * (patternMeth[it->first]);
        if (f >= maxf) {
            maxp = it->first;
            maxf = f;
        }
//        std::cout << it -> first << " "  << it->second << "*" << patternMeth[it->first] << " " << freq << std::endl;
    }
    meth=100.0 * double(count)/double(freq)/4.0;
    return 0;
}
int normalize(std::map<std::string, int> allp, std::map<std::string, double> & allpNormed, int freq)
{
    for (std::map <std::string, int>::const_iterator i=allp.begin(); i!=allp.end(); i++) {
        // normalize to 100 first
        allpNormed[i->first]=100.0*double(i->second)/double(freq);
    }
    return 0;
}
// normalization based on library size and uni-coverage across loci.
/*int nestedNormalize(std::map<std::string, int> allp1,int freq1,
                    std::map<std::string, int> allp2,int freq2,
                    double fnorm1k2,
                    std::map<std::string, double> & allpNormed1,
                    std::map<std::string, double> & allpNormed2)
{
    double fnorm2= 200.0/(freq1 + fnorm1k2 * freq2);
    for (std::map <std::string, int>::const_iterator i=allp1.begin(); i!=allp1.end(); i++) {
        allpNormed1[i->first]=fnorm2 * double(i->second);
        allpNormed2[i->first]=fnorm1k2 * fnorm2 * allp2[i->first];
    }
    return 0;
}*/
// group pattern vector by combining two pattern vectors
int groupPatterns(std::map <std::string, double> p1, std::map <std::string, double> p2, std::map <std::string, double> & p)
{
    for (std::map <std::string, double>::const_iterator i=p1.begin(); i!=p1.end(); i++) {
        p[i->first]=0.5* (i->second +p2.find(i->first)->second);
    }
    return 0;
}

// calculate entropy and normalized pattern for one sample
int entropy(std::map<std::string, double> allp,double & entropy, std::string & patterns)
{
    std::map <std::string, double> p;
    entropy=0;
    for (std::map <std::string, double>::const_iterator i=allp.begin(); i!=allp.end(); i++) {
        // entropy function
        double n=i->second;
        entropy=entropy+lgamma(n+1);
        // get normalized pattern
        std::ostringstream oss;
        oss << n << "\t";
        patterns+=oss.str();
    }
    return 0;
}

// calculate combinatorial entropy and normalized patterns
int combEntropy(std::string loci, 
                int freq1, 
                int freq2, 
                std::vector<std::string> p1, 
                std::vector<std::string> p2, 
                std::map <std::string, int> allpatterns,
                std::map <std::string, int> patternMeth,
                int methdiff,
				//double fnorm1k2,
                ogzstream & myfile)
{
    double meth1,meth2;
    std::string maxp1, maxp2;
    std::map <std::string, int> allp1=allpatterns;
    methPatterns(allp1, p1, patternMeth, meth1, maxp1);
    std::map <std::string, int> allp2=allpatterns;
    methPatterns(allp2, p2, patternMeth, meth2, maxp2);
    std::map <std::string, double> allpn, allpn1, allpn2;
    normalize(allp1, allpn1, freq1);
    normalize(allp2, allpn2, freq2);
	//nestedNormalize(allp1, freq1, allp2, freq2, fnorm1k2, allpn1, allpn2);
    groupPatterns(allpn1, allpn2, allpn);
    double e,e1,e2;
    std::string x,x1,x2;
    entropy(allpn, e, x);
    entropy(allpn1, e1, x1);
    entropy(allpn2, e2, x2);
    double combp= 2.0*e - e1 - e2;
    int methd=meth1-meth2;
    if (methd >= methdiff || methd <= - methdiff) {
        myfile << loci << "\t"<< combp << "\t" <<  freq1 << "\t" << freq2 << "\t" << meth1 << "\t" << meth2 << "\t" << maxp1 << "\t" << maxp2 << "\t" << x1 << x2 <<std::endl;
    }
    return 0;
}
// iterate lociMeth1 and lociMeth2, and print loci that shared.
int interLoci(std::map<std::string, std::vector<std::string> > lociMeth1, 
              std::map<std::string, std::vector<std::string> > lociMeth2, 
              std::map <std::string, int> allpatterns,
              std::map <std::string, int> patternMeth,
              int freq, int methdiff, ogzstream & myfile)
{
    std::map<std::string, std::vector<std::string> >::const_iterator it;
    for (it=lociMeth1.begin(); it!=lociMeth1.end(); it++)
    {
        std::string loci=it->first;
        if (lociMeth2.find(loci)!=lociMeth2.end()) {
            std::vector<std::string> p1=it->second;
            int freq1=p1.size();
            std::vector<std::string> p2=lociMeth2.find(loci)->second;
            int freq2=p2.size();
            if(freq1 >= freq & freq2 >= freq)
            {
                combEntropy(loci, freq1, freq2, p1, p2, allpatterns, patternMeth, methdiff, myfile);
            }
        }
    }
    return 0;
}
