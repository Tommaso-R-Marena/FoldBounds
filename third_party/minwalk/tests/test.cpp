#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "../src/HPFolding.hpp"

#include "testFilenames.h"

using namespace std;



const string testInputsFilePath = test13InputsFilePath;
const string testResultsFilePath = test13ResultsFilePath;

int main(int argc, char * argv[]) {
    ifstream testInputsFile;
    ifstream testResultsFile;
    
    
    
    testInputsFile.open(testInputsFilePath);
    testResultsFile.open(testResultsFilePath);
    
    string chain;
    string expectedScoreBuf;
    int expectedScore;
    int failed = 0;
    int passed = 0;
    
    if (testInputsFile.is_open() && testResultsFile.is_open()) {
        
        
        HPFolding hp;

        while (getline(testInputsFile, chain) && getline(testResultsFile, expectedScoreBuf)) {
            stringstream(expectedScoreBuf) >> expectedScore;
            //cout << passed << endl;
            hp.setChain(chain);
            hp.fold();
            if (hp.getScore() == expectedScore) {
                passed++;
            } else {
                cout << "Fail: " << chain << endl;
                cout << "\t Expected score: " << expectedScore << "\t Produced: " << hp.getScore() << endl;
                failed++;
            }
        }
        
        testInputsFile.close();
        testInputsFile.close();
        
        cout << "Result: " << passed << "/" << (passed+failed) << " successfull tests" << endl; // TODO: Time elapsed
    } else {
        cout << "Couldn't open test files." << endl;
    }
    
    
    return 0;
}

