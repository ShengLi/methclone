//
//  alignToMeth.cpp
//  methclone
//
//  Created by Sheng Li on 7/17/13.
//  Copyright (c) 2013 Weill Medical College of Cornell University. All rights reserved.
//

#include <iostream>
#include "api/BamReader.h"
#include <sstream>
#include <map>
#include <algorithm>

int byread(BamTools::BamAlignment al, 
           int d, 
           std::map<std::string, std::vector<std::string> > & lociMeth, 
           const BamTools::RefVector refs, 
           std::string sample)
{
    if (al.IsMapped() && al.HasTag("XM:Z")){
        std::vector<int> loc;
        std::vector<int> m;
        std::string Meth;
        
        loc.clear(); m.clear();// make sure that the vector is clear of content.
        al.GetTag("XM:Z:", Meth);
        std::string strand="+";
        if (al.IsReverseStrand()) {
            std::reverse(Meth.begin(), Meth.end());
            strand="-";
        }
        for (int i = 0; i < al.QueryBases.size(); i++) {
            if (Meth[i] == 'z')
            {
                loc.push_back(al.Position+i+1);
                m.push_back(0);
            } 
            else if (Meth[i] == 'Z')
            {
                loc.push_back(al.Position+i+1);
                m.push_back(1);
            } 
        } 
        if (loc.size()>=4) {
            for (int i=3 ; i <loc.size(); i++) {
                int dist=loc[i]-loc[i-3];
                if (dist<= d) {
                    std::ostringstream loci;
                    loci <<refs.at(al.RefID).RefName<< "\t" << loc[i-3] << "\t" << loc[i] << "\t" << sample << "\t" << loc[i]-loc[i-3]  << "\t" << strand << "\t" << loc[i-3]<<":"<<loc[i-2]<<":"<<loc[i-1]<<":"<<loc[i];
                    std::ostringstream meth;
                    meth << m[i-3]<<m[i-2] <<m[i-1] << m[i];
                    lociMeth[loci.str()].push_back(meth.str());
                    //std::cout << loci.str() << "\t" << meth.str()<< std::endl;
                }
            }
        }
    }
    
    return 0;
}
