//
//  testFilenames.h
//  HPFolding
//
//  Created by Stéphane Sercu on 1/11/16.
//  Copyright © 2016 Stéphane Sercu. All rights reserved.
//

//
//  This code is used to generate the test chains.
//  It generates all the possible chains of length 10 and 13
//  and calculates their optimal score.
//  To calculate the score, it uses the HPFolding class. So to
//  be certain that the scores are optimal, the 6 criteria must
//  be commented... yup that's sketchy.

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#include "../src/HPFolding.hpp"

#include "testFilenames.h"


using namespace std;

void generateTestChains(int length, ofstream* testChainsFile, ofstream* testResultsFile);

int main(int argc, char * argv[]) {
    
    remove(test10InputsFilePath.c_str());
    remove(test10ResultsFilePath.c_str());
    remove(test13InputsFilePath.c_str());
    remove(test13ResultsFilePath.c_str());
    
    ofstream testChains10File(test10InputsFilePath, ios_base::app | ios_base::out);
    ofstream testResults10File(test10ResultsFilePath, ios_base::app | ios_base::out);
    ofstream testChains13File(test13InputsFilePath, ios_base::app | ios_base::out);
    ofstream testResults13File(test13ResultsFilePath, ios_base::app | ios_base::out);
    
    
    if (testChains10File.is_open() && testResults10File.is_open()) {
        generateTestChains(10, &testChains10File, &testResults10File);
        cout << "Generated chains of length 10" << endl;
    }
    if (testChains13File.is_open() && testResults13File.is_open()) {
        generateTestChains(13, &testChains13File, &testResults13File);
        cout << "Generated chains of length 13" << endl;
    }
    
    
}


void generateTestChains(int length, ofstream* testChainsFile, ofstream* testResultsFile) {
    HPFolding hp;
    for (int i = 0; i < pow(2, length); i++) {
        string chain(length, 'H');
        for (int chainIndex = 0; chainIndex < length; chainIndex++) {
            if (i & (1<<chainIndex)) {
                chain[chainIndex] = 'P';
            }
        }
        *testChainsFile << chain << endl;
        hp.setChain(chain);
        hp.fold();
        *testResultsFile << hp.getScore() << endl;
        
        // Print progress
        if (i % 100 == 0) {
            cout << "Chain generation: " << i << "/" << pow(2, length) << endl;
        }
    }
}
