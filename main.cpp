// http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf

#include "fluidsolver.h"

#define WINDOW_SIDE 480
#define RECORDING true

#include <GL/gl.h>
#include <GL/glut.h>

#include <utility>

#include <iostream>
#include <fstream>
#include "fluidsolver.h"

using namespace std;

int s = WINDOW_SIDE; // window width and height
int window_id;

bool SIMULATION_PAUSED = false;

void init_gl();

void captureFrame(int framenum);
unsigned char *pRGB; // for capturing framebuffer data

void timerFunc(int);
void keyboardFunc(unsigned char, int, int);
void mouseFunc(int, int, int, int);
void addDensity(int, int);
void displayFunc();
void reshapeFunc(int, int);

int getColorFromDensity(int i, int j);
void displayNearestNeighbour();
int interpolateColors(float x, float y, int left, int top);
void displayBilinear();

void render();

bool vkey = false; // toggles velocity field display

int n = 120; // width and height of fluid box
int dg = WINDOW_SIDE / n; // cell dimensions
int dg_2 = dg / 2;
float dt = 0.2f;

enum UpscalingType
{
	NEAREST_NEIGHBOUR,
	BILINEAR
};

UpscalingType upscaling = NEAREST_NEIGHBOUR;
// UpscalingType upscaling = BILINEAR;

inline int I(int i, int j){ return i + (n + 2) * j; }

FluidSolver fs(n, dt);

int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  render();
  return 0;
}

void render()
{
  glutInitWindowSize(WINDOW_SIDE, WINDOW_SIDE);
  window_id = glutCreateWindow("Flame");
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);

  // register callbacks
  glutDisplayFunc(displayFunc);
  glutReshapeFunc(reshapeFunc);
  glutKeyboardFunc(keyboardFunc);
  glutMouseFunc(mouseFunc);
  glutMotionFunc(addDensity);

  glClearColor(1, 1, 1, 1);
  glColor3f(0, 0, 0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, WINDOW_SIDE, 0, WINDOW_SIDE);
  glutTimerFunc(1000, timerFunc, 0);
  glutMainLoop();
}

void timerFunc(int iter)
{
  if (!SIMULATION_PAUSED)
  {
    fs.step();
    iter++;

  	addDensity(WINDOW_SIDE/2, 3*WINDOW_SIDE/4);

  	if (RECORDING)
	  	captureFrame(iter);

    glutPostRedisplay();
  }

  glutTimerFunc(0, timerFunc, iter);
}

void displayFunc()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  switch(upscaling)
  {
    case NEAREST_NEIGHBOUR:
      displayNearestNeighbour();
      break;
    
    case BILINEAR:
      displayBilinear();
      break;
  }

  glutSwapBuffers();
}

void displayNearestNeighbour()
{
  for (int i = 1 ; i <= n ; i++)
  {
    int dx, dy;
    dx = (int) ((i - 0.5f) * dg);
    for (int j = 1 ; j <= n ; j++)
    {
      dy = (int)((j - 0.5f) * dg);
      if (fs.d[I(i, j)] > 0)
      {
        int c = getColorFromDensity(i, j);
        glBegin(GL_QUADS);
        glColor3ub(c, c, c);
        glVertex2f(dx-dg_2, dy-dg_2);
        glVertex2f(dx-dg_2, dy+dg_2);
        glVertex2f(dx+dg_2, dy+dg_2);
        glVertex2f(dx+dg_2, dy-dg_2);
        glEnd();
      }

      // draw velocity field
      if (vkey && i % 5 == 1 && j % 5 == 1)
      {
        int u = (int)(50 * fs.u[I(i,j)]);
        int v = (int)(50 * fs.v[I(i,j)]);
        glColor3ub(255,0,0);
        glBegin(GL_LINES);
        glVertex2f(dx,dy);
        glVertex2f(dx+u,dy+v);
        glEnd();
      }
    }
  }
}

// i = x coordinate, j = y coordinate
int getColorFromDensity(int i, int j)
{
  int c = (int)((1.0 - fs.d[I(i, j)]) * 255);
  if (c < 0) c = 0;
  return c;
}

void displayBilinear()
{
  for (int i = 0 ; i < WINDOW_SIDE ; i++)
  {
    float x = (float)i / (float)dg;
    int left = (int)x;
    for (int j = 0 ; j < WINDOW_SIDE ; j++)
    {
      float y = (float)j / (float)dg;
      int top = (int)y;
      int c = interpolateColors(x, y, left, top);
      glColor3ub(c,c,c);
      glBegin(GL_POINTS);
      glVertex2f(i, j);
      glEnd();
    }
  }
}

int interpolateColors(float x, float y, int left, int top)
{
  int color_topleft = getColorFromDensity(left, top);
  int color_topright = getColorFromDensity(left + dg, top);
  int color_bottomleft = getColorFromDensity(left, top + dg);
  int color_bottomright = getColorFromDensity(left + dg, top + dg);
  float t = (x - left) / (float)dg;
  float s = (y - top) / (float)dg;
  float result = (1-t)*(1-s)*color_topleft + t*(1-s)*color_topright + t*s*color_bottomright + (1-t)*s*color_bottomleft;
  return (int)result;
}

void keyboardFunc(unsigned char key, int x, int y)
{
  switch(key)
  {
    case 27: // esc
      glutDestroyWindow(window_id);
      break;

    case 32: // space
      SIMULATION_PAUSED = ! SIMULATION_PAUSED;
      break;

    case 'v': // velocity field
      vkey = ! vkey;
      break;
  }

}

void addDensity(int x, int y)
{
  if (x < 0 or y < 0 or x >= WINDOW_SIDE or y >= WINDOW_SIDE)
  {
  	cout << "not adding" << x << " " << y << endl;
    return;
  }

  int cellx = x / dg;
  int celly = (WINDOW_SIDE - y) / dg;
  fs.initFluidSingle(make_pair(cellx, celly));
}

void mouseFunc(int button, int state, int x, int y)
{
  switch(button)
  {
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_UP)
      {
        addDensity(x, y);
      }
  }
}

void reshapeFunc(int w, int h)
{
  if (h == 0) h = 1;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glOrtho(0.0, w, 0.0, h, -1., 1.);
  glViewport(0, 0, w, h);

  glutPostRedisplay();
}

// credit: Harsha
void captureFrame(int framenum)
{
  //global pointer char *pRGB
  pRGB = new unsigned char [3 * (WINDOW_SIDE+1) * (WINDOW_SIDE + 1) ];

  // set the framebuffer to read
  //default for double buffered
  glReadBuffer(GL_FRONT);

  glPixelStoref(GL_PACK_ALIGNMENT, 1); //for word alignment
  
  glReadPixels(0, 0, WINDOW_SIDE, WINDOW_SIDE, GL_RGB, GL_UNSIGNED_BYTE, pRGB);
  char filename[200];
  sprintf(filename, "frames/frame_%04d.ppm", framenum);
  ofstream out(filename, std::ios::out);
  out << "P6" << endl;
  out << WINDOW_SIDE << " " << WINDOW_SIDE << " 255" << endl;
  out.write(reinterpret_cast<char const *>(pRGB), (3 * (WINDOW_SIDE+1) * (WINDOW_SIDE + 1)) * sizeof(int));
  out.close();

  //function to store pRGB in a file named count
  delete pRGB;
}