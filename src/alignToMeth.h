//
//  alignToMeth.h
//  methclone
//
//  Created by Sheng Li on 7/17/13.
//  Copyright (c) 2013 Weill Medical College of Cornell University. All rights reserved.
//

#ifndef methclone_alignToMeth_h
#define methclone_alignToMeth_h

int byread(BamTools::BamAlignment al, int d, std::map<std::string, std::vector<std::string> > & lociMeth, const BamTools::RefVector refs, std::string sample);

#endif
