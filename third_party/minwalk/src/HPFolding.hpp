//  HPFolding.hpp
//  HPFolding

#ifndef HPFolding_hpp
#define HPFolding_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <functional>

#include "Grid.hpp"


using namespace std;

class HPFolding {
private:
    string chain;
    int len;
    
    vector<int> He, Ho; // {i|chain[i]=H && i%mod 2==0 (1)}
    vector<int> he, ho; // at index k: # of H in He (Ho) with i>k
    vector<int> PHOLinks, PHELinks; // # of H-P links before index [k] (where H is Odd / Even)
    
    int HeSize; // # of H with even index
    int HoSize; // # of H with odd index
    
    int HISize; // # of H with another H on its right
    int PSSize; // # of P-singlets
    
    int CO; // Maximum contacts achievable by beads of type O
    int CE;
    
    int CUpperBound; // # of expected contacts (upper-bound in MINWALK)
    int CParity;
    int CArea;
    
    
    void add(const int i, const int x, const int y, const int CI, const int AI);
    int calc_min_contacts();
    
    void analyse_chain();
    
public:
    HPFolding();
    HPFolding(const string chain);
    
    
    Grid grid;
    Grid finalGrid; // best grid
    
    void fold(function<void()> callback = NULL);
    
    int getScore() const;
    int getLength() const;
    string getChain() const;
    void setChain(const string chain);
    void show() const;
    
};

#endif /* HPFolding_hpp */
