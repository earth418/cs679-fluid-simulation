///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FLUID_2D_PERIODIC_H
#define FLUID_2D_PERIODIC_H

#ifdef USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "FIELD_2D.h"
#include "FLUID_2D.h"

class FLUID_2D_PERIODIC : public FLUID_2D {

public:
  FLUID_2D_PERIODIC(int xRes, int yRes, float dt);
  virtual ~FLUID_2D_PERIODIC() {};

protected:
  // solve linear system with Gauss-Seidel iteration, with periodic boundaries
  virtual void gaussSeidel(FIELD_2D& current, FIELD_2D& old);

  // advect field 'old' into 'current' using velocity field
  // 'xVelocity' and 'yVelocity' and periodic boundary conditions
  virtual void advect(FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity);

  // perform projection using periodic boundary conditions
  virtual void project();

  // step density field using periodic boundaries
  virtual void stepDensity();

  // step velocity field using periodic boundaries
  virtual void stepVelocity();

  // copy contents of 'field' to have periodic boundary conditions
  void setPeriodicBoundary(FIELD_2D& field);
};

#endif
