//
//  projetgraph.h
//  HPFolding
//
//  Created by Stéphane Sercu on 21/12/16.
//  Copyright © 2016 Stéphane Sercu. All rights reserved.
//

#ifndef projetgraph_h
#define projetgraph_h

#ifdef OSX
#include <GLUT/glut.h>
#elif LINUX
#include <GL/glut.h>
#endif

#include <stdlib.h>
#include <iostream>
#include <thread>
#include <math.h>
#include <sstream>


#include "HPFolding.hpp"


void projetgraph(int argc, char * argv[], HPFolding * hp);

#endif /* projetgraph_h */
