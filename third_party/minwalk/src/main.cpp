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
#include <cstdlib>

#include "HPFolding.hpp"


using namespace std;

int main(int argc, char * argv[]) {
    ifstream infile;
    ofstream outfile;
    int startCI = 0;
    bool noHCore = false;
    string inPath;
    string outPath;

    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        if (a == "--ci" && i + 1 < argc) {
            startCI = atoi(argv[++i]);
        } else if (a == "--no-hcore") {
            noHCore = true;
        } else if (inPath.empty()) {
            inPath = a;
        } else if (outPath.empty()) {
            outPath = a;
        }
    }

    if (!inPath.empty()) {
        infile.open(inPath);
    }
    if (!outPath.empty()) {
        outfile.open(outPath);
    }
    ostream & output = (outfile.is_open() ? outfile : cout);
    
    string chain;
    
    if (infile.is_open()) {
        HPFolding hp;
        if (startCI > 0) {
            hp.setStartCI(startCI);
        }
        hp.setDisableHCore(noHCore);
        
        function<void()> callback = NULL;
        if (graphic)
            callback = [&](){if (hp.getScore()>0) {hp.finalGrid.draw_conformation(output);}};
        
        while (getline(infile, chain)) {
            hp.setChain(chain);
            if (startCI > 0) {
                hp.setStartCI(startCI);
            }
            hp.setDisableHCore(noHCore);
            hp.fold(callback);
            output << chain << "\t" << hp.getScore() << endl;
            // Machine-readable FoldBounds witness (anchored coordinates).
            hp.finalGrid.dump_path(output);
        }
        infile.close();
    } else {
        // ask for chain
        cout << "Chain: \t";
        cin >> chain;
        
        HPFolding hp(chain);
        if (startCI > 0) {
            hp.setStartCI(startCI);
        }
        hp.setDisableHCore(noHCore);
        run(hp);
        
        output << chain << "\t" << hp.getScore() << endl;
    }
    
    
    return 0;
}
