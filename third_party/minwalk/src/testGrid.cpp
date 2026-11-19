//
//  testGrid.cpp
//  HPFolding
//
//  Created by Stéphane Sercu on 25/06/17.
//  Copyright © 2017 Stéphane Sercu. All rights reserved.
//

#include <stdio.h>
#include "Grid.hpp"

using namespace std;

int main(int argc, char * argv[]) {
    Grid g = Grid(10, "PHPPPHHHHH");
    g.setCell(5, 5, 0);
    g.draw_conformation();
    g.setCell(6, 5, 1);
    g.draw_conformation();
    g.setCell(7, 5, 2);
    g.draw_conformation();
    g.setCell(8, 5, 3);
    g.draw_conformation();
    g.setCell(9, 5, 4);
    g.draw_conformation();
    g.setCell(9, 4, 5);
    g.draw_conformation();
    g.setCell(8, 4, 6);
    g.draw_conformation();
    g.setCell(7, 4, 7);
    g.draw_conformation();
    g.setCell(6, 4, 8);
    g.draw_conformation();
}
