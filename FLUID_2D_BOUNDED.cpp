///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include "FLUID_2D_BOUNDED.h"
#include <cmath>


///////////////////////////////////////////////////////////////////////////
// Sets the edges such that (alongX) d/dy or (alongY) d/dx is zero
///////////////////////////////////////////////////////////////////////////
void setNeumannBoundary(FIELD_2D& field, bool alongX, bool alongY) {
  const int xr = field.xRes();
  const int yr = field.yRes();

  if (alongY)
    for (int j = 0; j < yr; ++j) { 
      field(0, j) = field(2, j); // copy thing two over from edge into the edge
      field(xr - 1, j) = field(xr - 3, j);
  }

  if (alongX)
    for (int i = 0; i < xr; ++i) {
      field(i, 0) = field(i, 2); // copy thing two over from edge into the edge
      field(i, yr - 1) = field(i, yr - 3);
  }
}

void ZeroBoundary(FIELD_2D& field, bool alongX, bool alongY) {
  const int xr = field.xRes();
  const int yr = field.yRes();

  if (alongY)
    for (int j = 0; j < yr; ++j) { 
        field(0, j) = 0; // copy thing two over from edge into the edge
        field(xr - 1, j) = 0.0;
    }

  if (alongX)
    for (int i = 0; i < xr; ++i) {
      field(i, 0) = 0; // copy thing two over from edge into the edge
      field(i, yr - 1) = 0.0;
    }
}

void setNeumann1away(FIELD_2D& field) {
  const int xr = field.xRes();
  const int yr = field.yRes();

  for (int j = 0; j < yr; ++j) { 
      field(0, j) = field(1, j); // copy thing two over from edge into the edge
      field(xr - 1, j) = field(xr - 2, j);
  }

  for (int i = 0; i < xr; ++i) {
    field(i, 0) = field(i, 1); // copy thing two over from edge into the edge
    field(i, yr - 1) = field(i, yr - 2);
  }
}

///////////////////////////////////////////////////////////////////////
// Constructor / Destructor
///////////////////////////////////////////////////////////////////////
FLUID_2D_BOUNDED::FLUID_2D_BOUNDED(int xRes, int yRes, float dt) :
  FLUID_2D(xRes, yRes, dt)
{
}

///////////////////////////////////////////////////////////////////////
// solve linear system with Gauss-Seidel iteration, with Dirichlet 
// boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::gaussSeidel(FIELD_2D& current, FIELD_2D& old)
{
  // IMPLEMENT ME
  for (int k = 0; k < 10; k++) 
  {
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++)
			  current(x,y) = (old(x,y) + current(x-1,y) + current(x+1,y) + current(x,y-1) + current(x,y+1)) * 0.25;

    // setPeriodicBoundary(current);
  // ZeroBoundary(current, true, true);
    // setNeumann1away(current);
	}
}



///////////////////////////////////////////////////////////////////////
// advect field 'old' into 'current' using velocity field
// 'xVelocity' and 'yVelocity' and Dirichlet boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::advect(FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity)
{
  // IMPLEMENT ME
  float dt0 = _dt * (_xRes - 2);

  for (int i = 1; i < _xRes - 1; ++i)
    for (int j = 1; j < _yRes - 1; ++j) {

      // Clamp stuff to the edges
      const float sampleX = max(0.0f, min((float) _xRes - 2, i - dt0 * xVelocity(i, j)));
      const float sampleY = max(0.0f, min((float) _yRes - 2, j - dt0 * yVelocity(i, j)));

      const int gridX = (int) sampleX;
      const int gridY = (int) sampleY;
      const float alphaX1 = sampleX - gridX; // interp x axis
      const float alphaY1 = sampleY - gridY;
      const float alphaX0 = 1 - alphaX1;
      const float alphaY0 = 1 - alphaY1;
    
      current(i, j) = alphaX0 * (alphaY0 * old(gridX    , gridY) + alphaY1 * old(gridX    , gridY + 1)) + 
                      alphaX1 * (alphaY0 * old(gridX + 1, gridY) + alphaY1 * old(gridX + 1, gridY + 1));

  }
}

///////////////////////////////////////////////////////////////////////
// perform projection using Neumann boundary conditions
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::project()
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

  // setNeumannBoundary(divergence, true, true);
  // setNeumannBoundary(pressure, true, true);
  setNeumann1away(divergence);
  setNeumann1away(pressure);

	gaussSeidel(pressure, divergence);
  // conjugateGradient(pressure, divergence);

  for (int y = 1; y < _yRes - 1; y++)
    for (int x = 1; x < _xRes - 1; x++)
    {
      _xVelocity(x, y) -= 0.5f * N * (pressure(x + 1, y) - pressure(x - 1, y));
      _yVelocity(x, y) -= 0.5f * N * (pressure(x, y + 1) - pressure(x, y - 1));
    }

  setNeumannBoundary(_xVelocity, false, true);
  setNeumannBoundary(_yVelocity, true, false);
}



///////////////////////////////////////////////////////////////////////
// Calculate and add vorticity force
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::addVorticity() {

    const float N = _xRes - 2;
    FIELD_2D vorticity(_xRes, _yRes);
    int i, j;

    for (i = 1; i < _xRes - 1; ++i)
        for (j = 1; j < _yRes - 1; ++j)
            vorticity(i, j) = (_yVelocity(i + 1, j) - _yVelocity(i - 1, j) - _xVelocity(i, j + 1) + _xVelocity(i, j - 1)) * (2 * N); 

    for (i = 1; i < _xRes - 1; ++i)
        for (j = 1; j < _yRes - 1; ++j) {

            const float etaX = (vorticity(i + 1, j) - vorticity(i - 1, j)) * (2 * N);
            const float etaY = (vorticity(i, j + 1) - vorticity(i, j - 1)) * (2 * N);

            const float n = sqrtf(etaX * etaX + etaY * etaY);

            const float Nx = etaX / n;
            const float Ny = etaY / n;

            // N x w == w is just in the z direction, and N is in the x-y plane, so N x w gets us something in the x-y plane of the form
            // (Ny * W, -Nx * w)
            const float tempMul = vorticity(i, j) * vorticityEpsilon * N;
            // if (tempMul > 0.01)
            //   std::cout << tempMul << '\n';
            // std::cout << "(" << Nx << ", " << Ny << ")\n";

            // std::cout << tempMul;
            // addForce(i, j, tempMul * Ny, -tempMul * Nx);
            _xVelocityOld(i, j) = tempMul * Ny;
            _yVelocityOld(i, j) = -tempMul * Nx;
    }

}


///////////////////////////////////////////////////////////////////////
// step density field using Neumann boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::stepDensity()
{
  // ZeroBoundary(_density, true, true);
  ZeroBoundary(_xVelocity, false, true);
  ZeroBoundary(_yVelocity, true, false);

  // for (int i = 0; i < _xRes - 0; ++i)
  //       for (int j = 0; j < _yRes - 0; ++j)
  //           _xVelocityOld(i, j) = _yVelocityOld(i, j) = 0.0;

  // addVorticity();


  addSource(_density, _densityOld);
  swapFields(_density, _densityOld);
  advect(_density, _densityOld, _xVelocity, _yVelocity);
  setNeumannBoundary(_density, true, true);
}

///////////////////////////////////////////////////////////////////////
// step velocity field using Neumann boundaries
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::stepVelocity()
{
  // ZeroBoundary(_xVelocityOld, true, true);
  // ZeroBoundary(_yVelocityOld, true, true);
  // addVorticity();

  addSource(_xVelocity, _xVelocityOld); 
  addSource(_yVelocity, _yVelocityOld);
	project();

	swapFields(_xVelocityOld, _xVelocity); 
  swapFields(_yVelocityOld, _yVelocity);

	advect(_xVelocity, _xVelocityOld, _xVelocityOld, _yVelocityOld);
  setNeumannBoundary(_xVelocity, true, false);

  advect(_yVelocity, _yVelocityOld, _xVelocityOld, _yVelocityOld);
  setNeumannBoundary(_yVelocity, false, true);
	
  project();
}

///////////////////////////////////////////////////////////////////////
// solve the system using conjugate gradients
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::conjugateGradient(FIELD_2D& pressure, FIELD_2D& divergence)
{
  // IMPLEMENT ME
  int x, y;
  FIELD_2D r(_xRes, _yRes);
  FIELD_2D d(_xRes, _yRes);
  float delta_new = 0.0;

  for (x = 1; x < _xRes - 1; ++x)
    for (y = 1; y < _yRes - 1; ++y) {
      const float curr_r = divergence(x, y) - 
      (5 * pressure(x, y) - pressure(x + 1, y) - pressure(x - 1, y) - pressure(x, y - 1) - pressure(x, y + 1));

      r(x, y) = curr_r;
      d(x, y) = curr_r;
      delta_new += curr_r * curr_r;
  }

  const float delta_0 = delta_new;

  int i = 0;
  const int i_max = 50;
  const float eps = 1e-6;
  const float maxD = delta_0 * eps * eps;
  FIELD_2D q(_xRes, _yRes);


  while (i < i_max && delta_new > maxD) {

    // q = Ad
    // alpha = dNew / (d.transpose * q)
    float alpha = 0.0;
    for (x = 1; x < _xRes - 1; ++x)
      for (y = 1; y < _yRes - 1; ++y) {
        const float curr_d = d(x, y);
        const float curr_q = 4 * curr_d - d(x + 1, y) - d(x - 1, y) - d(x, y - 1) - d(x, y + 1);
        q(x, y) = curr_q;

        alpha += curr_q * curr_d;
    }
    
    if (fabs(alpha) > 0.0f)
      alpha = delta_new / alpha;
    
    float delta_old = delta_new;
    delta_new = 0.0;

    // x = x + alpha * d;
    // r = r - alpha * q;
    // deltaNew = transpose(r) * r
    for (x = 1; x < _xRes - 1; ++x)
      for (y = 1; y < _yRes - 1; ++y) {
        pressure(x, y) += alpha * d(x, y);
        r(x, y) -= alpha * q(x, y);

        delta_new += r(x, y) * r(x, y);
    }

    float beta = delta_new / delta_old;
    // d = r + beta * d
    for (x = 1; x < _xRes - 1; ++x)
      for (y = 1; y < _yRes - 1; ++y) {
        d(x, y) = r(x, y) + beta * d(x, y);
    }

    ++i;
  }

  cout << "Converged in " << i << " iterations to " << delta_new << endl;
  
}

///////////////////////////////////////////////////////////////////////
// solve the system using conjugate gradients
///////////////////////////////////////////////////////////////////////
void FLUID_2D_BOUNDED::preconditionedConjugateGradient(FIELD_2D& pressure, FIELD_2D& divergence)
{
  // IMPLEMENT ME
}
