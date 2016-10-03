//
//  bamByChr.h
//  methclone
//
//  Created by Sheng Li on 7/17/13.
//  Copyright (c) 2013 Weill Medical College of Cornell University. All rights reserved.
//

#ifndef methclone_bamByChr_h
#define methclone_bamByChr_h
#include "api/BamReader.h"

int bamCheck(std::string bamFile, BamTools::BamReader & reader);
int readerToMeth(BamTools::BamReader & reader, BamTools::RefVector::const_iterator i, int d, std::map<std::string, std::vector<std::string> > & lociMeth, const BamTools::RefVector refs, std::string sample);
int bamToLociMeth(std::string bamFile1, 
                  std::string bamFile2, 
                  const char * outFile, 
                  std::string sample, 
                  int d, int freq, int methdiff);
#endif
