//
//  projectgraph.cpp
//  HPFolding
//
//  Created by Stéphane Sercu on 10/11/16.
//  Copyright © 2016 Stéphane Sercu. All rights reserved.
//

#include "HPGUI.h"

// ================================
// CONFIGURATION
// --------------------------------

const int TIMERMSECS = 40; // actualisation rate
const float FRAME_SIZE = 800.0; // Size of the grid (px)

const bool SHOW_GRID = true;
const bool SHOW_WEIGHTS = true;
const bool SHOW_HCORE = true;

// =================================


int win_w = 0;
int win_h = 0;
int frame_size = 0;

int grid_size = 0;
int cell_size = 0;

bool searching = true;

string chain;

HPFolding *hp = NULL;
Grid *g = NULL;


using namespace std;
void init_window(int argc, char * argv[]) {
    glutInit(&argc,argv);
    
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    
    glutInitWindowSize (win_w, win_h);
    glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH)-win_w)/2, (glutGet(GLUT_SCREEN_HEIGHT)-win_h)/2);
    glutCreateWindow ("HPFolding");
    
    glClearColor(255.0,255.0,255.0,0.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D (0.0, win_w, 0.0,win_h);
    
    

}
void draw_filled_circle(GLfloat x, GLfloat y, GLfloat radius){
    // https://gist.github.com/strife25/803118
    int i;
    int triangleAmount = 20; //# of triangles used to draw circle
    
    GLfloat twicePi = 2.0f * M_PI;
    
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y); // center of circle
    for(i = 0; i <= triangleAmount;i++) {
        glVertex2f(
                   x + (radius * cos(i *  twicePi / triangleAmount)),
                   y + (radius * sin(i * twicePi / triangleAmount))
                   );
    }
    glEnd();
}

void draw_filled_square(GLfloat x, GLfloat y, GLfloat size) {
    // (x, y) = bottom left of the square
    glBegin(GL_QUADS);
    glVertex3f((float)x, (float)y, 0.0f);
    glVertex3f((float)x+size, (float)y, 0.0f);
    glVertex3f((float)x+size,y+size, 0.0f);
    glVertex3f((float)x, (float)y+size, 0.0f);
    glEnd();
}



void draw_grid() {
    for(int i = 0; i <= grid_size; i ++)
    {
        glBegin(GL_LINES);
        glVertex3f((float)i*cell_size+1, 1.0f, 0.0f);
        glVertex3f((float)i*cell_size+1, (float)grid_size*cell_size+1, 0.0f);
        glVertex3f(1.0f, (float)i*cell_size+1, 0.0f);
        glVertex3f((float)grid_size*cell_size+1, (float)i*cell_size+1, 0.0f);
        glEnd();
        
    }
}

void draw_score() {
    string s= "Score: " + to_string(g->getScore());
    //ss<< "Score: " << g->getScore();
    glRasterPos2i(frame_size+10.0,frame_size-25.0);
    for (int i=0; i<(int)s.size(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
}

void draw_hcore() {
    glLineWidth(2.0f);
    glColor3f(0.0,0.0,0.0);
    glBegin(GL_LINES);
    glVertex3f((float)g->getHCore().xmin*cell_size+1, (float)g->getHCore().ymin*cell_size +1, 0.0f);
    glVertex3f(((float)g->getHCore().xmax+1)*cell_size+1, (float)g->getHCore().ymin*cell_size +1, 0.0f);
    
    glVertex3f(((float)g->getHCore().xmax+1)*cell_size+1, (float)g->getHCore().ymin*cell_size +1, 0.0f);
    glVertex3f(((float)g->getHCore().xmax+1)*cell_size+1, ((float)g->getHCore().ymax+1)*cell_size +1, 0.0f);
    
    glVertex3f(((float)g->getHCore().xmax+1)*cell_size+1, ((float)g->getHCore().ymax+1)*cell_size +1, 0.0f);
    glVertex3f((float)g->getHCore().xmin*cell_size+1, ((float)g->getHCore().ymax+1)*cell_size +1, 0.0f);
    
    glVertex3f((float)g->getHCore().xmin*cell_size+1, ((float)g->getHCore().ymax+1)*cell_size +1, 0.0f);
    glVertex3f((float)g->getHCore().xmin*cell_size+1, (float)g->getHCore().ymin*cell_size +1, 0.0f);
    glEnd();
}

void draw_protein() {
    glColor3f(0.0,0.0,0.0);
    glRasterPos2i((win_w-18*chain.size())/2.0,frame_size+45);
    for (int i=0; i<(int)chain.size(); i++)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, chain[i]);

}

void draw_legend() {
    glColor3f(0.0,0.0,0.0);
    string s ="Even Exp. Sites";
    glRasterPos2i(frame_size+20.0+20.0,55.0);
    for (int i=0; i<(int)s.size(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, s[i]);
    
    s = "Odd Exp. Sites";
    glRasterPos2i(frame_size+20.0+20.0,30.0);
    for (int i=0; i<(int)s.size(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, s[i]);
    
    glColor3f(1.0-0.15,1.0,1.0);
    draw_filled_square((float)frame_size+10.0, 50.0, 25.0);
    
    
    glColor3f(1.0,1.0,1.0-0.15);
    draw_filled_square((float)frame_size+10.0, 25.0, 25.0);
}
void draw_latice() {
    for (int j=0; j<grid_size; j++) {
        for (int i=0; i<grid_size; i++) {
            int index = g->getCell(i, j);
            if (index != -1) {
                // beads
                if (chain[index] == 'H') {
                    glColor3f(1.0,0.0,0.0);
                    draw_filled_circle((float)i*cell_size + cell_size/2, (float)j*cell_size + cell_size/2, cell_size/4);
                    glColor3f(1.0,1.0,1.0);
                    glRasterPos2i((float)i*cell_size + cell_size/2 - 5, (float)j*cell_size + cell_size/2 - 5);
                    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, 'H');
                    
                } else {
                    glColor3f(0.0,0.0,1.0);
                    draw_filled_square((float)i*cell_size+ cell_size/4,(float) j*cell_size+ cell_size/4, (float)cell_size/2);
                    glColor3f(1.0,1.0,1.0);
                    glRasterPos2i((float)i*cell_size + cell_size/2 - 5, (float)j*cell_size + cell_size/2 - 5);
                    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, 'P');
                    
                }
                
                // Connectors
                if (index < (int)chain.size()-1) {
                    glBegin(GL_LINES);
                    glColor3f(0.0,0.0,0.0);
                    // Find the next amino
                    if (g->getCell(i+1, j) == index + 1) { // left
                        glVertex3f((float)(i+1)*cell_size + cell_size/4, (float)j*cell_size + cell_size/2, 0.0f);
                        glVertex3f((float)(i)*cell_size + 3*cell_size/4, (float)j*cell_size + cell_size/2, 0.0f);
                    } else if (g->getCell(i-1, j) == index + 1) { // right
                        glVertex3f((float)(i-1)*cell_size + 3*cell_size/4, (float)j*cell_size + cell_size/2, 0.0f);
                        glVertex3f((float)(i)*cell_size + cell_size/4, (float)j*cell_size + cell_size/2, 0.0f);
                    } else if (g->getCell(i, j+1) == index + 1) { // top
                        glVertex3f((float)i*cell_size + cell_size/2, (float)(j+1)*cell_size + cell_size/4, 0.0f);
                        glVertex3f((float)i*cell_size + cell_size/2, (float)j*cell_size + 3*cell_size/4, 0.0f);
                    }else if (g->getCell(i, j-1) == index + 1) { // down
                        glVertex3f((float)i*cell_size + cell_size/2, (float)(j-1)*cell_size + 3*cell_size/4, 0.0f);
                        glVertex3f((float)i*cell_size + cell_size/2, (float)j*cell_size + cell_size/4, 0.0f);
                    }
                    glEnd();
                }
            }
            
            // Weights
            if (SHOW_WEIGHTS) {
                if (g->getWeight(i, j) > 0) { // TODO: legend colors
                    glColor3f(1.0-0.15*g->getWeight(i, j),1.0,1.0);
                    draw_filled_square((float)i*cell_size+1, (float)j*cell_size+1, cell_size-1);
                } else if (g->getWeight(i, j) < 0) {
                    glColor3f(1.0,1.0,1.0+0.15*g->getWeight(i, j));
                    draw_filled_square((float)i*cell_size+1, (float)j*cell_size+1, cell_size-1);
                }
            }
        }
    }

}

void draw_length() {
    string s = " Length: " + to_string(chain.size());
    glRasterPos2i(frame_size+10.0,frame_size-50.0);
    for (int i=0; i<(int)s.size(); i++)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
}


void draw()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glColor3f(0.0, 0.0, 0.0);
    glViewport(20.0,20.0,win_w,win_h);
    glPointSize(1.0f);
    glLineWidth(1.0f);
    
    if (SHOW_GRID)
        draw_grid();
    
    draw_score();
    
    draw_length();
    
    draw_protein();
    
    draw_latice();
    
    
    
    if (SHOW_HCORE) {
        draw_hcore();
        draw_legend();
    }
    
    
    glutSwapBuffers();
}
void animate(int val)
{
    glutTimerFunc(TIMERMSECS,animate,0);
    draw();
}

void projetgraph(int argc, char * argv[], HPFolding * _hp) {
    
    
    //*g = hp.finalGrid;
    hp = _hp;
    chain = hp->getChain();
    g = &hp->grid;
    grid_size = g->getSize();
    cell_size = ceil(FRAME_SIZE/grid_size);
    frame_size = cell_size*grid_size;
    win_w = frame_size+201;
    win_h = frame_size+101;
    
    init_window(argc, argv);
    thread t(&HPFolding::fold, hp, [](){g=&hp->finalGrid; g->draw_conformation();});
    
    glutTimerFunc(TIMERMSECS,animate,0);
    glutDisplayFunc(draw);
    glutMainLoop();
    
    t.join();
}

