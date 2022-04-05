///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "solver.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FLUID_2D_BOUNDED_H
#define FLUID_2D_BOUNDED_H

#ifdef USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "FIELD_2D.h"
#include "FLUID_2D.h"

class FLUID_2D_BOUNDED : public FLUID_2D {

public:
  FLUID_2D_BOUNDED(int xRes, int yRes, float dt);
  virtual ~FLUID_2D_BOUNDED() {};

  float vorticityEpsilon = 5.0f;

private:
  // solve linear system with Gauss-Seidel iteration
  virtual void gaussSeidel(FIELD_2D& current, FIELD_2D& old);

  // advect field 'old' into 'current' using velocity field
  // 'xVelocity' and 'yVelocity'
  virtual void advect(FIELD_2D& current, FIELD_2D& old, FIELD_2D& xVelocity, FIELD_2D& yVelocity);

  // perform projection
  virtual void project();

  // step density field
  virtual void stepDensity();

  // step velocity field
  virtual void stepVelocity();

  void addVorticity();

  // solve the system using conjugate gradients
  void conjugateGradient(FIELD_2D& pressure, FIELD_2D& divergence);

  // solve the system using preconditioned conjugate gradients
  void preconditionedConjugateGradient(FIELD_2D& pressure, FIELD_2D& divergence);
};

#endif
