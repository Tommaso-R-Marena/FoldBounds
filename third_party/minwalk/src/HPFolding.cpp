//
//  HPFolding.cpp
//  HPFolding
//
//  Created by Stéphane Sercu on 31/10/16.
//  Copyright © 2016 Stéphane Sercu. All rights reserved.
//

#include "HPFolding.hpp"

HPFolding::HPFolding(const string chain):   chain(chain),
                                            len((int)chain.length()),
                                            grid(Grid(0.5*len, chain)),
                                            finalGrid(Grid(0.5*len, chain))

{

    if (len==0)
        return;
    analyse_chain();
}

HPFolding::HPFolding() : HPFolding("") {}

void HPFolding::analyse_chain() {
    /**
     * Determine all the characteristic of the chain (Ho, ho, PHOLinks, etc)
     */
    
    len = (int)chain.length();
    
    He = vector<int>(len); // len is the maximum size, the actual size is in HeSize/HoSize
    Ho = vector<int>(len); // They are juste preallocated.
    
    he = vector<int>(len); // Here len is the real size
    ho = vector<int>(len);
    
    PHELinks = vector<int>(len);
    PHOLinks = vector<int>(len);
    
    PHELinks[0] = 0;
    PHOLinks[0] = 0;
    
    HeSize = 0; HoSize = 0; HISize = 0; PSSize = 0;
    
    for (int i=0; i<len; i++) {
        if (chain[i] == 'H') {
            // H's are categorised in odd/even
            if ((i+1) % 2) {
                Ho[HoSize] = i;
                HoSize++;
            } else {
                He[HeSize] = i;
                HeSize++;
            }
            
            // H's with another H behind them are counted
            if (i<len-1) {
                if (chain[i+1] == 'H') {
                    HISize++;
                }
            }
            
        } else if (chain[i] == 'P') {
            // P-singlets are counted
            if (i>0 && i<len-1) {
                if (chain[i-1] == 'H' && chain[i+1] == 'H') {
                    PSSize++;
                }
            }
        } else {
            throw invalid_argument("Bad chain: only 'H' and 'P' are allowed.");
        }
        
        
        if (i>0 && chain[i-1] != chain[i]) { // If P-H link
            
            
            if (i%2==0) {
                if (chain[i]=='P') {
                    PHELinks[i] = PHELinks[i-1]+1;
                    PHOLinks[i] = PHOLinks[i-1];
                } else {
                    PHOLinks[i] = PHOLinks[i-1]+1;
                    PHELinks[i] = PHELinks[i-1];
                }
            } else {
                if (chain[i]=='P') {
                    PHOLinks[i] = PHOLinks[i-1]+1;
                    PHELinks[i] = PHELinks[i-1];
                } else {
                    PHELinks[i] = PHELinks[i-1]+1;
                    PHOLinks[i] = PHOLinks[i-1];
                }
                
            }
            
        } else {
            if (i>0) {
                PHELinks[i] = PHELinks[i-1];
                PHOLinks[i] = PHOLinks[i-1];
            }
            
        }
    }
    CE = (2*HeSize + (He[HeSize-1]==len-1));
    CO = (2*HoSize+(Ho[0]==0) + (Ho[HoSize-1]==len-1));
    
    
    // Build he & ho once for all
    int lastHeIndex = 0;
    int lastHoIndex = (Ho[0] == 0);
    he[0] = HeSize;
    ho[0] = HoSize - (Ho[0] == 0 && HoSize > 0);
    
    for (int k=1; k<len; k++) {
        if (k >= He[lastHeIndex] && HeSize > 0 && k <= He[HeSize-1]) {
            lastHeIndex++;
            he[k] = he[k-1] - 1;
        } else {
            he[k] = he[k-1];
        }
        
        if (k >= Ho[lastHoIndex] && HoSize > 0 && k <= Ho[HoSize-1]) {
            lastHoIndex++;
            ho[k] = ho[k-1] - 1;
        } else {
            ho[k] = ho[k-1];
        }
    }
    
    CUpperBound = calc_min_contacts();
}
void HPFolding::fold(function<void()> callback) {
    if (len>2) {
        // place the first amino on the center
        grid.setCell(grid.getCenter(), grid.getCenter(), 0);
    
        // place the second amino at its right
        grid.setCell(grid.getCenter() + 1, grid.getCenter(), 1);
    
        int CI = CUpperBound;
        int AI = 0;
        while (finalGrid.getScore() == 0 && CI > 0) {
            if (CI == CArea) {
                int semi_perim = ceil( 2.0*(float)(sqrt( (float)(HoSize+HeSize+PSSize) )));
                AI = floor((float)(semi_perim/2.0) )*ceil((float)(semi_perim/2.0));
            } else {
                AI = grid.getSize()*grid.getSize(); // = infinity
            }
            add(2, grid.getCenter() + 1, grid.getCenter(), CI, AI);
            CI--;
        }
    }
    
    if (callback)
        callback();
}

int HPFolding::calc_min_contacts() {
    /**
     * Calculates the maximum number of H-contacts achievable by the chain.
     * Assumes that the chain has been analysed (analyseChain())
     */
    CParity = min(2*HeSize+(He[HeSize-1]==len-1), 2*HoSize + (Ho[HoSize-1]==len-1) + (Ho[0]==0));
    CArea = 2*(HoSize+HeSize) - HISize - ceil( (float)(2*sqrt( (float)((HoSize+HeSize)+PSSize) )) );
    return min(CParity, CArea);
}


void HPFolding::add(const int i, const int x, const int y, const int CI, const int AI) {
    /**
     *  (x, y) is the position of the last positionned amino
     *  i is the index of the next amino to place (aroud (x,y))
     *  CI is the minimum score we are looking for
     *  AI is the maximum Area the h-core can have
    */
    
    // If all aminos have been placed...
    if (i >= len) {
            // if current score is best...
            if (grid.getScore() > finalGrid.getScore())
                finalGrid = grid;
        return; // WIN
    }
    
    if (finalGrid.getScore() == CI) {
        // We stop at the first optimal conformation
        return; // WIN
    }
    
    // Criterion 1
    int WoMax = CO-CI+PHOLinks[i];
    int WeMax = CE-CI+PHELinks[i];
    if (grid.getWe() > WeMax || grid.getWo() > WoMax) {
        return; // LOOSE
    }
    
    // Criterion 2
    int SE = grid.getSe(1) + grid.getSe(2) + grid.getSe(3) + grid.getSe(4);
    int SO = grid.getSo(1) + grid.getSo(2) + grid.getSo(3) + grid.getSo(4);
    //grid.draw_conformation();
    if (ho[i-2]<SE-(WeMax-grid.getWe()) || he[i-2]<SO-(WoMax-grid.getWo())) {
        return; // LOOSE
    }
    
    // Criterion 3
    int MO3 = min(grid.getSo(3), (int)(He[HeSize-1] == len-1));
    int MO2 = min(grid.getSo(2), he[i-2] - MO3);
    int MO1 = min(grid.getSo(1), he[i-2] - MO2 - MO3);
    int Cleft = 2*ho[i-2] + (Ho[HoSize-1] == len-1) + 3*MO3 + 2*MO2 + MO1;
    if (Cleft < (CI-grid.getScore())) {
        return; // LOOSE
    }
    
    // Criterion 4
    int c_area = (grid.getHCore().xmax - grid.getHCore().xmin + 1)*(grid.getHCore().ymax - grid.getHCore().ymin + 1);
    if (c_area > AI) {
        return; // LOOSE
    }
    
    // Criterion 5
    if (grid.getNSPInHcore()>AI-(HoSize+HeSize+PSSize)) {
        return; // LOOSE
    }
    
    // Criterion 6
    if (CI == CArea) {
        // Checks the current line
        int nh = 0;
        for (int i=0; i<grid.getSize(); i++) {
            nh += ((nh+1)%2)*(grid.getAmino(i, y) == 'H') + (nh%2)*(grid.getAmino(i, y) == 'P');
            if (nh>2) {
                return; // LOOSE
            }
        }
        // Checks the current column
        nh = 0;
        for (int j=0; j<grid.getSize(); j++) {
            nh += ((nh+1)%2)*(grid.getAmino(x, j) == 'H') + (nh%2)*(grid.getAmino(x, j) == 'P');
            if (nh>2) {
                return; // LOOSE
            }
        }
    }
    
    

    
    // LEFT
    if (grid.getCell(x-1, y) == -1) { // if left exists and is empty
        grid.setCell(x-1, y, i);
        add(i+1, x-1, y, CI, AI);
        grid.setCell(x-1, y, -1);
    }
    // RIGHT
    if (grid.getCell(x+1, y) == -1) { // if right exists and is empty
        grid.setCell(x+1, y, i);
        add(i+1, x+1, y, CI, AI);
        grid.setCell(x+1, y, -1);
    }
    
    // TOP
    if (grid.getCell(x, y-1) == -1) { // if top exists and is empty
        grid.setCell(x, y-1, i);
        add(i+1, x, y-1, CI, AI);
        grid.setCell(x, y-1, -1);
    }
    
    // DOWN
    if (grid.getCell(x, y+1) == -1) { // if down exists and is empty
        grid.setCell(x, y+1, i);
        add(i+1, x, y+1, CI, AI);
        grid.setCell(x, y+1, -1);
    }
    return; // LOOSE
    
}

int HPFolding::getScore()  const{
    return finalGrid.getScore();
}

int HPFolding::getLength()  const{
    return len;
}

string HPFolding::getChain()  const{
    return chain;
}
void HPFolding::setChain(const string chain) {
    this->chain = chain;
    analyse_chain();
    finalGrid = Grid(3*len, chain);
    grid = Grid(3*len, chain);
}

void HPFolding::show()  const{
    cout << "Score: " << finalGrid.getScore() << endl;
    finalGrid.draw_conformation();
}

