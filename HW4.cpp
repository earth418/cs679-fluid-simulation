///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "demo.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cstdio>
#ifdef USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "FLUID_2D.h"
#include "FLUID_2D_PERIODIC.h"
#include "FLUID_2D_BOUNDED.h"
#include "FIELD_2D.h"

int mouseDown[3];
int oldMouseX, oldMouseY, mouseX, mouseY;

int   windowX      = 512;
int   windowY      = 512;
bool  drawVelocity = false;
int   res          = 128;
float force        = 5.0;
float source       = 100.0;

// the fluid simulation object
FLUID_2D* fluid;

///////////////////////////////////////////////////////////////////////////////
// Keyboard command processing function
///////////////////////////////////////////////////////////////////////////////
void keyboardCallback(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 'c':
		case 'C':
      fluid->clear();
			break;

		case 'q':
		case 'Q':
			exit(0);
			break;

		case 'v':
		case 'V':
			drawVelocity = !drawVelocity;
			break;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Mouse position and click processing function
///////////////////////////////////////////////////////////////////////////////
void mouseCallback(int button, int state, int x, int y)
{
	oldMouseX = mouseX = x;
	oldMouseX = mouseY = y;

	mouseDown[button] = state == GLUT_DOWN;
}

///////////////////////////////////////////////////////////////////////////////
// Mouse movement processing function
///////////////////////////////////////////////////////////////////////////////
void motionCallback(int x, int y)
{
	mouseX = x;
	mouseY = y;
}

///////////////////////////////////////////////////////////////////////////////
// Window shaping function
///////////////////////////////////////////////////////////////////////////////
void reshapeCallback(int width, int height)
{
	glutReshapeWindow(width, height);

	windowX = width;
	windowY = height;
}

///////////////////////////////////////////////////////////////////////////////
// Idle command processing function
///////////////////////////////////////////////////////////////////////////////
void idleCallback()
{
  // add a packet of smoke at the bottom
  fluid->addSource();

  // add a buoyancy force from the smoke packet
  fluid->addBuoyancy();

  // clear the old field so we can store the input forces and
  // sources in them
  fluid->clearOlds();

  // if one of the mouse buttons are pressed, process their input
	if (mouseDown[0] || mouseDown[2])
  {
    int i = (int)((           mouseX /(float)windowX) * res + 1);
    int j = (int)(((windowY - mouseY)/(float)windowY) * res + 1);

    // make sure the click is in bounds
    if (i > 1 && i < res - 2 && j > 1 && j < res - 2)
    {
      // left mouse -- add a force
      if (mouseDown[0]) 
        fluid->addForce(i, j, 
                       force * (mouseX - oldMouseX), 
                       force * (oldMouseY - mouseY));

      // right mouse -- add a density
      if (mouseDown[2]) 
        fluid->addDensity(i,j, source);

      oldMouseX = mouseX;
      oldMouseY = mouseY;
    }
  }

  // step the fluid simulation
  fluid->step();

	glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////////////
// The drawing function
///////////////////////////////////////////////////////////////////////////////
void displayCallback()
{
	glViewport(0, 0, windowX, windowY);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity ();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

  if (drawVelocity) 
    fluid->drawVelocity();
  else 
    fluid->drawDensity();

	glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv)
{
	glutInit(&argc, argv);
	
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowSize(windowX, windowY);
	glutCreateWindow(argv[0]);

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers ();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers ();

	glutKeyboardFunc(keyboardCallback);
	glutMouseFunc(mouseCallback);
	glutMotionFunc(motionCallback);
	glutReshapeFunc(reshapeCallback);
	glutIdleFunc(idleCallback);
	glutDisplayFunc(displayCallback);
  
  fluid = new FLUID_2D_BOUNDED(res, res, 0.1);
//   fluid = new FLUID_2D_PERIODIC(res, res, 0.1);

	glutMainLoop ();

  return 0;
}
