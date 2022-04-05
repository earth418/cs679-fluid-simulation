///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include "FLUID_2D_PERIODIC.h"

///////////////////////////////////////////////////////////////////////
// Constructor / Destructor
///////////////////////////////////////////////////////////////////////
FLUID_2D_PERIODIC::FLUID_2D_PERIODIC(int xRes, int yRes, float dt) :
  FLUID_2D(xRes, yRes, dt)
{
}

///////////////////////////////////////////////////////////////////////
// copy contents of 'field' to have periodic boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_PERIODIC::setPeriodicBoundary(FIELD_2D& field)
{
  int N = _xRes - 2;
	for (int i = 1; i <= N; i++)
  {
		field(0,     i) = field(N, i);
		field(N + 1, i) = field(1, i);
		field(i,     0) = field(i, N);
		field(i, N + 1) = field(i, 1);
	}
	field(0, 0    ) = field(N, N);
	field(0, N + 1) = field(N, 1);
	field(N + 1, 0  ) = field(1, N);
	field(N + 1, N + 1) = field(1, 1);
}

///////////////////////////////////////////////////////////////////////
// solve linear system with Gauss-Seidel iteration, with periodic 
// boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_PERIODIC::gaussSeidel(FIELD_2D& current, FIELD_2D& old)
{
	for (int k = 0; k < 10; k++) 
  {
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
			  current(x,y) = (old(x,y) + current(x-1,y) + current(x+1,y) + current(x,y-1) + current(x,y+1)) * 0.25;

    setPeriodicBoundary(current);
	}
}

///////////////////////////////////////////////////////////////////////
// advect field 'old' into 'current' using velocity field
// 'xVelocity' and 'yVelocity' and periodic boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_PERIODIC::advect(FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity)
{
  int N = _xRes - 2;
	float dt0 = _dt * N;

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      // trace backwards through the velocity field
		  float tempX = x - dt0 * xVelocity(x,y); 
      float tempY = y - dt0 * yVelocity(x,y);

      // wrap around for periodic boundary conditions in x direction
      int   xInt    = (int)tempX;
      float xFrac = tempX - xInt;
      if (xInt > N) 
        xInt = xInt % N;
      if (xInt < 1) {
        xInt = xInt % N;
        xInt += N - 1;
      }
      tempX = xInt + xFrac;

      // wrap around for periodic boundary conditions in y direction
      int   yInt    = (int)tempY;
      float yFrac = tempY - yInt;
      if (yInt > N) 
        yInt = yInt % N;
      if (yInt < 1) {
        yInt = yInt % N;
        yInt += N - 1;
      }
      tempY = yInt + yFrac;

      // retrieve the coordinates of the grid cells to interpolate
      int x0 = (int) tempX; 
      int y0 = (int) tempY; 
      int x1 = x0 + 1;
      int y1 = y0 + 1;

      // compute the interpolation weights
		  float s1 = tempX - x0; 
      float s0 = 1 - s1; 
      float t1 = tempY - y0; 
      float t0 = 1 - t1;

      // compute the final interpolation
		  current(x,y) = s0 * (t0 * old(x0,y0) + t1 * old(x0, y1)) +
					           s1 * (t0 * old(x1,y0) + t1 * old(x1, y1));
    }
}

///////////////////////////////////////////////////////////////////////
// perform projection using periodic boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_PERIODIC::project()
{
  int N = _xRes - 2;
  FIELD_2D& pressure = _xVelocityOld;
  FIELD_2D& divergence = _yVelocityOld;

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      divergence(x,y) = -0.5f * (_xVelocity(x + 1,y)  - _xVelocity(x - 1,y) + 
                                 _yVelocity(x, y + 1) - _yVelocity(x, y - 1)) / N;
      pressure(x,y) = 0;
    }
  setPeriodicBoundary(divergence);
  setPeriodicBoundary(pressure);
	gaussSeidel(pressure, divergence);

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      _xVelocity(x, y) -= 0.5f * N * (pressure(x + 1, y) - pressure(x - 1, y));
      _yVelocity(x, y) -= 0.5f * N * (pressure(x, y + 1) - pressure(x, y - 1));
    }

  setPeriodicBoundary(_xVelocity);
  setPeriodicBoundary(_yVelocity);
}

///////////////////////////////////////////////////////////////////////
// step density field using periodic boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_PERIODIC::stepDensity()
{
  addSource(_density, _densityOld);
  swapFields(_density, _densityOld);
  advect(_density, _densityOld, _xVelocity, _yVelocity);
  setPeriodicBoundary(_density);
}

///////////////////////////////////////////////////////////////////////
// step velocity field using periodic boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_PERIODIC::stepVelocity()
{
	addSource(_xVelocity, _xVelocityOld); 
  addSource(_yVelocity, _yVelocityOld);
	project();

	swapFields(_xVelocityOld, _xVelocity); 
  swapFields(_yVelocityOld, _yVelocity);

	advect(_xVelocity, _xVelocityOld, _xVelocityOld, _yVelocityOld); 
  setPeriodicBoundary(_xVelocity);

  advect(_yVelocity, _yVelocityOld, _xVelocityOld, _yVelocityOld);
  setPeriodicBoundary(_yVelocity);
	project();
}
