//
//  main.cpp
//  HPFolding
//
//  Created by Stéphane Sercu on 29/10/16.
//  Copyright © 2016 Stéphane Sercu. All rights reserved.
//


#ifdef GRAPH
bool graphic = true;
#define run(obj) projetgraph(argc, argv, &hp)
#include "HPGUI.h"
#else
bool graphic = false;
#define run(obj) obj.fold()
#endif


#include <iostream>
#include <string>
#include <fstream>
#include <functional>

#include "HPFolding.hpp"


using namespace std;

int main(int argc, char * argv[]) {
    ifstream infile;
    ofstream outfile;
    
    if (argc > 1) {
        infile.open(argv[1]); // First arg is the input file
    }
    if (argc > 2) {
        outfile.open(argv[2]); // Second arg is the output file
    }
    ostream & output = (outfile.is_open() ? outfile : cout);
    
    string chain;
    
    if (infile.is_open()) {
        HPFolding hp;
        
        function<void()> callback = NULL;
        if (graphic)
            callback = [&](){if (hp.getScore()>0) {hp.finalGrid.draw_conformation(output);}};
        
        while (getline(infile, chain)) {
            hp.setChain(chain);
            hp.fold(callback);
            output << chain << "\t" << hp.getScore() << endl;
        }
        infile.close();
    } else {
        // ask for chain
        cout << "Chain: \t";
        cin >> chain;
        
        HPFolding hp(chain);
        run(hp);
        
        output << chain << "\t" << hp.getScore() << endl;
    }
    
    
    return 0;
}
