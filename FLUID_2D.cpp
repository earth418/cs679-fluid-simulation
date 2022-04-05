///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include "FLUID_2D.h"

///////////////////////////////////////////////////////////////////////
// Constructor / Destructor
///////////////////////////////////////////////////////////////////////
FLUID_2D::FLUID_2D(int xRes, int yRes, float dt) :
  _xRes(xRes), _yRes(yRes), _dt(dt),
  _density(xRes, yRes),
  _densityOld(xRes, yRes),
  _xVelocity(xRes, yRes),
  _xVelocityOld(xRes, yRes),
  _yVelocity(xRes, yRes),
  _yVelocityOld(xRes, yRes)
{
}

///////////////////////////////////////////////////////////////////////
// swap the pointers for fields 'left' and 'right'
///////////////////////////////////////////////////////////////////////
void FLUID_2D::swapFields(FIELD_2D& left, FIELD_2D& right)
{
  float*& leftData = left.data();
  float*& rightData = right.data();

  float* temp = leftData;
  leftData = rightData;
  rightData = temp;
}

///////////////////////////////////////////////////////////////////////
// add the contents of 'source' to 'field'
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addSource(FIELD_2D& field, FIELD_2D& source)
{
  for (int y = 0; y < _xRes; y++)
    for (int x = 0; x < _yRes; x++)
      field(x,y) += _dt * source(x,y);
}

///////////////////////////////////////////////////////////////////////
// user input force
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addForce(int x, int y, float xForce, float yForce)
{
  _xVelocityOld(x,y) = xForce;
  _yVelocityOld(x,y) = yForce;
}

///////////////////////////////////////////////////////////////////////
// clear all fields
///////////////////////////////////////////////////////////////////////
void FLUID_2D::clear()
{
  _xVelocity.clear();
  _yVelocity.clear();
  _density.clear();

  _xVelocityOld.clear();
  _yVelocityOld.clear();
  _densityOld.clear();
}

///////////////////////////////////////////////////////////////////////
// clear old fields
///////////////////////////////////////////////////////////////////////
void FLUID_2D::clearOlds()
{
  _xVelocityOld.clear();
  _yVelocityOld.clear();
  _densityOld.clear();
}

///////////////////////////////////////////////////////////////////////
// draw the density field to GL
///////////////////////////////////////////////////////////////////////
void FLUID_2D::drawDensity()
{
	float h = 1.0f / _xRes;

	glBegin(GL_QUADS);
		for (int i = 1; i < _xRes; i++)  
    {
			float x = (i - 0.5f) * h;
			for (int j = 1; j < _yRes; j++) 
      {
				float y = (j - 0.5f) * h;

  		  float density = _density(i,j);

				glColor3f(density, density, density);
        glVertex2f(x, y);
        glVertex2f(x + h, y);
        glVertex2f(x + h, y + h);
        glVertex2f(x, y + h);
			}
		}
	glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the velocity field to GL
///////////////////////////////////////////////////////////////////////
void FLUID_2D::drawVelocity()
{
	float h = 1.0f / _xRes;

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

		for(int i = 1; i < _xRes; i++) 
    {
			float x = (i - 0.5f) * h;
			for(int j = 1; j < _yRes; j++) 
      {
				float y = (j - 0.5f) * h;
        float xDiff = _xVelocity(i,j);
        float yDiff = _yVelocity(i,j);

				glVertex2f(x, y);
				glVertex2f(x + xDiff, y + yDiff);
			}
		}

	glEnd();
}

///////////////////////////////////////////////////////////////////////
// timestep the fluid simulation with periodic boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D::step()
{
  stepVelocity();
  stepDensity();
}

///////////////////////////////////////////////////////////////////////
// add a smoke packet to the center
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addSource()
{
  float dx = 1.0 / _xRes;
  float dy = 1.0 / _yRes;
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
    {
      float xReal = x * dx;
      float yReal = y * dy;

      if (xReal > 0.475 && xReal < 0.525 &&
          yReal > 0.1 && yReal < 0.15)
      {
        _density(x,y) = 1.0;
      }
    }
}

///////////////////////////////////////////////////////////////////////
// add a buoyancy force from the density
///////////////////////////////////////////////////////////////////////
void FLUID_2D::addBuoyancy()
{
  for (int y = 0; y < _yRes; y++)
    for (int x = 0; x < _xRes; x++)
      _yVelocity(x,y) += 0.001 * _density(x,y);
}
