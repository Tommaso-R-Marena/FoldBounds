//
//  Grid.cpp
//  HPFolding
//
//  Created by Stéphane Sercu on 1/11/16.
//  Copyright © 2016 Stéphane Sercu. All rights reserved.
//

#include "Grid.hpp"


// =================================================
// BasicGrid
// =================================================

BasicGrid::BasicGrid(int size, string prot): size(size), protein(prot), score(0){
    if (size <0)
        size=0;
    cells = vector<int>(size*size, -1);
}




void BasicGrid::draw_conformation(ostream& out)  const {
    for (int j=0; j<size;j++) {
        for (int i=0;i<size; i++) {
            out << getAmino(i, j);
            if (i<size-1) {
                if (abs(getCell(i, j)-getCell(i+1, j))==1 && getCell(i+1, j) != -1 && getCell(i, j) != -1) {
                    out << "-";
                } else {
                    out << " ";
                }
            } else {
                out << " ";
            }
        }
        out << endl;
        
        for (int i=0; i<size;i++) {
            if (j<size-1) {
                if (abs(getCell(i, j)-getCell(i, j+1))==1 && getCell(i, j+1) != -1 && getCell(i, j) != -1) {
                    out << "| ";
                } else {
                    out << "  ";
                }
            }
        }
        out << endl;
    }
    out<<endl;
}


bool BasicGrid::setCell(int x, int y, int amino_ind)
{
    if (amino_ind < (int)protein.size() && amino_ind >= -1) {
        cells[y*size+x] = amino_ind;
        return true;
    }
    return false;
}

int BasicGrid::getCell(int x, int y)  const {
    if (x>=0 && x<size && y>=0 && y< size) {
        return cells[y*size+x];
    }
    return -2; // OUT
}

// =================================================
// Grid
// =================================================


Grid::Grid(int size, string prot) : BasicGrid(size, prot){
    
    exp_sites = vector<int>(size*size, 0);
    
    so[0]=0; so[1]=0; so[2]=0, so[3]=0;
    se[0]=0; se[1]=0; se[2]=0, se[3]=0;
    
    Wo = 0; We = 0;
    
    hcore.xmin = 0; hcore.xmax = 0; hcore.ymax = 0; hcore.ymin = 0;
    hcoreHistory.push_back(hcore);
    
    NSP_in_hcore = 0;
    int nsp = 0;
    if (prot.size()>1) {
        nsp = (prot[0] == 'P' && prot[1] == 'P'); // # of non singlet P
        
        nsp+=(prot[prot.size()-1]=='P' && prot[prot.size()-2]=='P');
        for (int i=1; i<(int)prot.size()-1;i++) {
            nsp += (prot[i] == 'P' && (prot[i-1] == 'P' || prot[i+1]=='P'));
        }
    }
    NSPInHcore = vector<Point>();
    NSPOutHcore = vector<Point>();
    
    first_h_index = (int)prot.size();
    int j = 0;
    while (j<(int)prot.size() && first_h_index==(int)prot.size()) {
        if (prot[j] == 'H') {
            first_h_index = j;
        }
        j++;
    }
}



bool Grid::setCell(int x, int y, int amino_ind) {
    if (amino_ind < (int)protein.size() && amino_ind >= -1) {
        
        bool odd = false; // if the AA is odd
        if (amino_ind==-1) { // Has to be calculated before the cell is changed!
            odd = (getCell(x, y)+1)%2; //
        } else {
            odd = (amino_ind+1)%2;
        }
        
        int h_contacts = (getAmino(x-1, y) == 'H') +
        (getAmino(x, y-1) == 'H') +
        (getAmino(x, y+1) == 'H') +
        (getAmino(x+1, y) == 'H'); // # of H contact
        
        int p_contacts = (getAmino(x-1, y) == 'P') +
        (getAmino(x, y-1) == 'P') +
        (getAmino(x, y+1) == 'P') +
        (getAmino(x+1, y) == 'P'); // # of H contact
        
        if (amino_ind == -1) {
            if (getAmino(x, y) == 'H') { // remove an H
                
                
                
                // dec score
                score -= (h_contacts - (protein[getCell(x, y)-1] == 'H')); // Assumes that beads are placed in the order they are in the protein
                
                // H-core update (back in history)
                hcore = hcoreHistory[hcoreHistory.size()-2]; 
                hcoreHistory.erase(hcoreHistory.end()-1);
                
                NSPInHcore.erase(remove_if(NSPInHcore.begin(), NSPInHcore.end(), [this](Point p){
                    if (!this->hcore.contains(p.x, p.y)) {
                        this->NSPOutHcore.push_back(p);
                        return true;
                    }
                    return false;
                }), NSPInHcore.end());
                
                // Update Weight at (x, y)
                if (!odd) {
                    if (h_contacts>0) {
                        so[h_contacts-1]++;
                    }
                } else {
                    if (h_contacts>0) {
                        se[h_contacts-1]++;
                    }
                }
                setWeight(x, y, (!odd ? -1: 1)*h_contacts);
                
                // update Weights around (x, y)
                if (getCell(x+1, y)==-1) {decrement_weight(x+1, y, odd);}
                if (getCell(x-1, y)==-1) {decrement_weight(x-1, y, odd);}
                if (getCell(x, y+1)==-1) {decrement_weight(x, y+1, odd);}
                if (getCell(x, y-1)==-1) {decrement_weight(x, y-1, odd);}
                
                
                // update wasted contact
                if (odd) {
                    Wo -= p_contacts;
                } else {
                    We -= p_contacts;
                }
                
            } else if (getAmino(x, y) == 'P') { // remove a P
                // Update Weight at (x, y)
                if (!odd) {
                    if (h_contacts>0) {
                        so[h_contacts-1]++;
                    }
                } else {
                    if (h_contacts>0) {
                        se[h_contacts-1]++;
                    }
                }
                setWeight(x, y, (!odd ? -1: 1)*h_contacts);
                
                // update wasted contact
                if (odd) {
                    We -= h_contacts;
                } else {
                    Wo -= h_contacts;
                }
                
                // Is it a non-singlets P beads & will it be removed from the h-core?
                int ind = getCell(x, y);
                if ( ( ((ind>0) && (protein[ind-1] == 'P')) || ((ind<(int)protein.size()-1) && (protein[ind+1] == 'P')) )) {
                    if (this->hcore.contains(x, y)) {
                        NSPInHcore.erase(remove(NSPInHcore.begin(), NSPInHcore.end(), Point(x,y)), NSPInHcore.end());
                    } else {
                        NSPOutHcore.erase(remove(NSPOutHcore.begin(), NSPOutHcore.end(), Point(x,y)), NSPOutHcore.end());
                    }
                }                
                
            }
        } else { // amino_ind is a valid index in protein
            if (protein[amino_ind] == 'H') { // place an H
                
                // inc score
                score += (h_contacts - (protein[amino_ind-1] == 'H'));
                
                // H-core update
                if (amino_ind==first_h_index) {
                    hcore.xmin = x;
                    hcore.xmax = x;
                    hcore.ymax = y;
                    hcore.ymin = y;
                } else {
                    if (x<hcore.xmin) {
                        hcore.xmin = x;
                    }
                    
                    if (x>hcore.xmax) {
                        hcore.xmax = x;
                    }
                    
                    if (y<hcore.ymin) {
                        hcore.ymin = y;
                    }
                    
                    if (y>hcore.ymax) {
                        hcore.ymax = y;
                    }

                }
                
                hcoreHistory.push_back(hcore);
                
                
                NSPOutHcore.erase(remove_if(NSPOutHcore.begin(), NSPOutHcore.end(), [this](Point p){
                    if (this->hcore.contains(p.x, p.y)) {
                        this->NSPInHcore.push_back(p);
                        return true;
                    }
                    return false;
                }), NSPOutHcore.end());
                
                
                // Update Weight at (x, y)
                zero_weight(x, y, !odd);
                
                // update Weights around (x, y)
                if (getCell(x+1, y)==-1) {increment_weight(x+1, y, odd);}
                if (getCell(x-1, y)==-1) {increment_weight(x-1, y, odd);}
                if (getCell(x, y+1)==-1) {increment_weight(x, y+1, odd);}
                if (getCell(x, y-1)==-1) {increment_weight(x, y-1, odd);}
                
                // update wasted contact
                if (odd) {
                    Wo += p_contacts;
                } else {
                    We += p_contacts;
                }
                
            } else if (protein[amino_ind] == 'P') { // place a P
                // Update Weight at (x, y)
                zero_weight(x, y, !odd);
                
                // update wasted contact
                if (odd) {
                    We += h_contacts;
                } else {
                    Wo += h_contacts;
                }
                
                // Is it a non-singlets P beads & will it be placed in the h-core?
                if ((!isPSiglet(amino_ind) && protein[amino_ind] == 'P')) {
                    if (this->hcore.contains(x, y) ) {
                        NSPInHcore.push_back(Point(x,y));
                    } else {
                        NSPOutHcore.push_back(Point(x,y));
                    }
                }
            }
        }
        
        
        
        cells[y*size+x] = amino_ind;
        
        return true;
    } else {
        return false;
    }
}


bool Grid::isPSiglet(int ind) {
    return ( (ind>0) && (ind<(int)protein.size()) && protein[ind] == 'P' ) // is a P
            && ( (ind>0) && protein[ind-1] == 'H' ) // previous is an H
            && ( (ind<(int)protein.size()-1) && (protein[ind+1] == 'H') ); // next is an H
}


char BasicGrid::getAmino(int x, int y)  const {
    int ind = getCell(x, y);
    if (ind > -1) {
        if (ind >= 0 and ind < (int)protein.size()) {
            return protein[ind];
        }
    } else if (ind == -2) {
        return OUT_CELL;
    } else if (ind == -1) {
        return FREE_CELL;
    }
    return '?';
    
}

int Grid::getWeight(int x, int y) {
    return exp_sites[y*size+x];
}

void Grid::setWeight(int x, int y, int w) {
    exp_sites[y*size+x] = w; // TODO: checks if (x, y) are correct and if 0<=w<=4 ?
}
void Grid::increment_weight(int x, int y, bool odd) {
    // Assumes that so (& se) aren't changed by another way than decrement & increment_weight
    if (odd) {
        if (getWeight(x, y) != 0) {
            so[-1*getWeight(x, y)-1] --;
        }
        so[-1*getWeight(x, y)] ++;
    } else {
        if (getWeight(x, y) != 0) {
            se[getWeight(x, y)-1] --;
        }
        se[getWeight(x, y)] ++;
    }
    
    exp_sites[y*size+x] += (odd ? -1: 1); // TODO: checks if (x, y) are correct and if 0<=w<=4 ?
}
void Grid::decrement_weight(int x, int y, bool odd) {
    // Assumes that so (& se) aren't changed by another way than decrement & increment_weight
    if (odd) {
        so[-1*getWeight(x, y)-1] --;
        if (getWeight(x, y)+1 != 0) {
            so[-1*getWeight(x, y)-2] ++;
        }
    } else {
        se[getWeight(x, y)-1] --;
        if (getWeight(x, y)-1 != 0) {
            se[getWeight(x, y)-2] ++;
        }
    }
    
    
    exp_sites[y*size+x] -= (odd ? -1: 1); // TODO: checks if (x, y) are correct and if 0<=w<=4 ?
}

void Grid::zero_weight(int x, int y, bool odd) {
    // Assumes that so (& se) aren't changed by another way than decrement & increment_weight
    if (getWeight(x, y) != 0) { // If the exposed site will disapear...
        if (odd) {
            so[-1*getWeight(x, y)-1] -=1;
        } else {
            se[getWeight(x, y)-1] -=1;
        }
    }
    
    exp_sites[y*size+x] = 0;
}

int BasicGrid::getScore() const {
    return score;
}
int BasicGrid::getSize() const {
    return size;
}

int Grid::getSe(int w) {
    if (w>0 && w<=4) {
        return se[w-1];
    }
    return 0;
}
int Grid::getSo(int w) {
    if (w>0 && w<=4) {
        return so[w-1];
    }
    return 0;
}
Rect Grid::getHCore() {
    return hcore;
}

int BasicGrid::getCenter() const {
    return floor((float)size/2);
}

int Grid::getWe() {
    return We;
}
int Grid::getWo() {
    return Wo;
}

int Grid::getNSPInHcore() {
    return (int)NSPInHcore.size();
}


// DEBUG
void Grid::showAminos() {
    for (int j=0; j<size;j++) {
        for (int i=0;i<size; i++) {
            cout << getAmino(i, j) << " ";
        }
        cout << endl;
    }
    cout<<endl;
}

void Grid::showWeights() {
    for (int j=0; j<size;j++) {
        for (int i=0;i<size; i++) {
            cout << (getWeight(i, j)<0 ? "-" : " ") << abs(getWeight(i, j)) << " ";
        }
        cout << endl;
    }
    cout<<endl;
}