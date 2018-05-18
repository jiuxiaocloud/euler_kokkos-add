#ifndef HYDRO_INIT_FUNCTORS_3D_H_
#define HYDRO_INIT_FUNCTORS_3D_H_

#include <limits> // for std::numeric_limits
#ifdef __CUDA_ARCH__
#include <math_constants.h> // for cuda math constants, e.g. CUDART_INF
#endif // __CUDA_ARCH__

#include "shared/kokkos_shared.h"
#include "HydroBaseFunctor3D.h"

// init conditions
#include "shared/problems/BlastParams.h"

namespace euler_kokkos { namespace muscl {
    
/*************************************************/
/*************************************************/
/*************************************************/
class InitFakeFunctor3D : public HydroBaseFunctor3D {
  
public:
  InitFakeFunctor3D(HydroParams params,
		    DataArray3d Udata) :
    HydroBaseFunctor3D(params), Udata(Udata)  {};
  
  // static method which does it all: create and execute functor
  static void apply(HydroParams params,
                    DataArray3d Udata,
		    int         nbCells)
  {
    InitFakeFunctor3D functor(params, Udata);
    Kokkos::parallel_for(nbCells, functor);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& index) const
  {

    const int isize = params.isize;
    const int jsize = params.jsize;
    const int ksize = params.ksize;
    
    int i,j,k;

    index2coord(index,i,j,k,isize,jsize,ksize);
        
    Udata(i  ,j  ,k  , ID) = 0.0;
    Udata(i  ,j  ,k  , IP) = 0.0;
    Udata(i  ,j  ,k  , IU) = 0.0;
    Udata(i  ,j  ,k  , IV) = 0.0;
    Udata(i  ,j  ,k  , IW) = 0.0;
    
  } // end operator ()

  DataArray3d Udata;

}; // InitFakeFunctor3D

/*************************************************/
/*************************************************/
/*************************************************/
class InitImplodeFunctor3D : public HydroBaseFunctor3D {
  
public:
  InitImplodeFunctor3D(HydroParams params,
		       DataArray3d Udata) :
    HydroBaseFunctor3D(params), Udata(Udata)  {};
  
  // static method which does it all: create and execute functor
  static void apply(HydroParams params,
                    DataArray3d Udata,
		    int         nbCells)
  {
    InitImplodeFunctor3D functor(params, Udata);
    Kokkos::parallel_for(nbCells, functor);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& index) const
  {

    const int isize = params.isize;
    const int jsize = params.jsize;
    const int ksize = params.ksize;
    const int ghostWidth = params.ghostWidth;
    
#ifdef USE_MPI
    const int i_mpi = params.myMpiPos[IX];
    const int j_mpi = params.myMpiPos[IY];
    const int k_mpi = params.myMpiPos[IZ];
#else
    const int i_mpi = 0;
    const int j_mpi = 0;
    const int k_mpi = 0;
#endif

    const int nx = params.nx;
    const int ny = params.ny;
    const int nz = params.nz;

    const real_t xmin = params.xmin;
    const real_t ymin = params.ymin;
    const real_t zmin = params.zmin;

    const real_t dx = params.dx;
    const real_t dy = params.dy;
    const real_t dz = params.dz;
    
    const real_t gamma0 = params.settings.gamma0;
    
    int i,j,k;
    index2coord(index,i,j,k,isize,jsize,ksize);
    
    real_t x = xmin + dx/2 + (i+nx*i_mpi-ghostWidth)*dx;
    real_t y = ymin + dy/2 + (j+ny*j_mpi-ghostWidth)*dy;
    real_t z = zmin + dz/2 + (k+nz*k_mpi-ghostWidth)*dz;
    
    real_t tmp = x + y + z;
    if (tmp > 0.5 && tmp < 2.5) {
      Udata(i  ,j  ,k  , ID) = 1.0;
      Udata(i  ,j  ,k  , IP) = 1.0/(gamma0-1.0);
      Udata(i  ,j  ,k  , IU) = 0.0;
      Udata(i  ,j  ,k  , IV) = 0.0;
      Udata(i  ,j  ,k  , IW) = 0.0;
    } else {
      Udata(i  ,j  ,k  , ID) = 0.125;
      Udata(i  ,j  ,k  , IP) = 0.14/(gamma0-1.0);
      Udata(i  ,j  ,k  , IU) = 0.0;
      Udata(i  ,j  ,k  , IV) = 0.0;
      Udata(i  ,j  ,k  , IW) = 0.0;	    
    }
    
  } // end operator ()

  DataArray3d Udata;

}; // InitImplodeFunctor3D
  
/*************************************************/
/*************************************************/
/*************************************************/
class InitBlastFunctor3D : public HydroBaseFunctor3D {

public:
  InitBlastFunctor3D(HydroParams params,
		     BlastParams bParams,
		     DataArray3d Udata) :
    HydroBaseFunctor3D(params), bParams(bParams), Udata(Udata)  {};
  
  // static method which does it all: create and execute functor
  static void apply(HydroParams params,
		    BlastParams bParams,
                    DataArray3d Udata,
		    int         nbCells)
  {
    InitBlastFunctor3D functor(params, bParams, Udata);
    Kokkos::parallel_for(nbCells, functor);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& index) const
  {

    const int isize = params.isize;
    const int jsize = params.jsize;
    const int ksize = params.ksize;
    const int ghostWidth = params.ghostWidth;
    
#ifdef USE_MPI
    const int i_mpi = params.myMpiPos[IX];
    const int j_mpi = params.myMpiPos[IY];
    const int k_mpi = params.myMpiPos[IZ];
#else
    const int i_mpi = 0;
    const int j_mpi = 0;
    const int k_mpi = 0;
#endif

    const int nx = params.nx;
    const int ny = params.ny;
    const int nz = params.nz;

    const real_t xmin = params.xmin;
    const real_t ymin = params.ymin;
    const real_t zmin = params.zmin;

    const real_t dx = params.dx;
    const real_t dy = params.dy;
    const real_t dz = params.dz;
    
    const real_t gamma0 = params.settings.gamma0;

    // blast problem parameters
    const real_t blast_radius      = bParams.blast_radius;
    const real_t radius2           = blast_radius*blast_radius;
    const real_t blast_center_x    = bParams.blast_center_x;
    const real_t blast_center_y    = bParams.blast_center_y;
    const real_t blast_center_z    = bParams.blast_center_z;
    const real_t blast_density_in  = bParams.blast_density_in;
    const real_t blast_density_out = bParams.blast_density_out;
    const real_t blast_pressure_in = bParams.blast_pressure_in;
    const real_t blast_pressure_out= bParams.blast_pressure_out;
  

    int i,j,k;
    index2coord(index,i,j,k,isize,jsize,ksize);
    
    real_t x = xmin + dx/2 + (i+nx*i_mpi-ghostWidth)*dx;
    real_t y = ymin + dy/2 + (j+ny*j_mpi-ghostWidth)*dy;
    real_t z = zmin + dz/2 + (k+nz*k_mpi-ghostWidth)*dz;

    real_t d2 = 
      (x-blast_center_x)*(x-blast_center_x)+
      (y-blast_center_y)*(y-blast_center_y)+
      (z-blast_center_z)*(z-blast_center_z);    
    
    if (d2 < radius2) {
      Udata(i  ,j  ,k  , ID) = blast_density_in;
      Udata(i  ,j  ,k  , IP) = blast_pressure_in/(gamma0-1.0);
      Udata(i  ,j  ,k  , IU) = 0.0;
      Udata(i  ,j  ,k  , IV) = 0.0;
      Udata(i  ,j  ,k  , IW) = 0.0;
    } else {
      Udata(i  ,j  ,k  , ID) = blast_density_out;
      Udata(i  ,j  ,k  , IP) = blast_pressure_out/(gamma0-1.0);
      Udata(i  ,j  ,k  , IU) = 0.0;
      Udata(i  ,j  ,k  , IV) = 0.0;
      Udata(i  ,j  ,k  , IW) = 0.0;
    }
    
  } // end operator ()
  
  BlastParams bParams;
  DataArray3d Udata;
  
}; // InitBlastFunctor3D

/*************************************************/
/*************************************************/
/*************************************************/
/**
 * Test of the Rayleigh-Taylor instability.
 * See
 * http://www.astro.princeton.edu/~jstone/Athena/tests/rt/rt.html
 * for a description of such initial conditions
 */
class RayleighTaylorInstabilityFunctor3D : public HydroBaseFunctor3D {

public:
  RayleighTaylorInstabilityFunctor3D(HydroParams params,
				     RayleighTaylorInstabilityParams rtiparams,
				     DataArray3d Udata,
				     VectorField3d gravity) :
    HydroBaseFunctor3D(params),
    rtiparams(rtiparams),
    Udata(Udata),
    gravity(gravity)
  {};

  // static method which does it all: create and execute functor
  static void apply(HydroParams params,
		    RayleighTaylorInstabilityParams rtiparams,
                    DataArray3d Udata,
		    VectorField3d gravity)
  {
    uint64_t nbCells = params.isize * params.jsize * params.ksize;
    RayleighTaylorInstabilityFunctor3D functor(params, rtiparams, Udata, gravity);
    Kokkos::parallel_for(nbCells, functor);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& index) const
  {

    const int isize = params.isize;
    const int jsize = params.jsize;
    const int ksize = params.ksize;
    const int ghostWidth = params.ghostWidth;
    
#ifdef USE_MPI
    const int i_mpi = params.myMpiPos[IX];
    const int j_mpi = params.myMpiPos[IY];
    const int k_mpi = params.myMpiPos[IZ];
#else
    const int i_mpi = 0;
    const int j_mpi = 0;
    const int k_mpi = 0;
#endif

    const int nx = params.nx;
    const int ny = params.ny;
    const int nz = params.nz;

    const real_t xmin = params.xmin;
    const real_t ymin = params.ymin;
    const real_t zmin = params.zmin;

    const real_t xmax = params.xmax;
    const real_t ymax = params.ymax;
    const real_t zmax = params.zmax;

    const real_t Lx = xmax-xmin;
    const real_t Ly = ymax-ymin;
    const real_t Lz = zmax-zmin;

    const real_t dx = params.dx;
    const real_t dy = params.dy;
    const real_t dz = params.dz;
    
    const real_t gamma0 = params.settings.gamma0;
  
    int i,j,k;
    index2coord(index,i,j,k,isize,jsize,ksize);
    
    real_t x = xmin + dx/2 + (i+nx*i_mpi-ghostWidth)*dx;
    real_t y = ymin + dy/2 + (j+ny*j_mpi-ghostWidth)*dy;
    real_t z = zmin + dz/2 + (k+nz*k_mpi-ghostWidth)*dz;

    /* initialize perturbation amplitude */
    real_t amplitude = rtiparams.amplitude;
    real_t        d0 = rtiparams.d0;
    real_t        d1 = rtiparams.d1;
    
    //bool  randomEnabled = rti.randomEnabled;


    /* uniform static gravity field */
    const real_t gravity_x = rtiparams.gx;
    const real_t gravity_y = rtiparams.gy;
    const real_t gravity_z = rtiparams.gz;
    
    real_t         P0 = 1.0;
      
  
    // the initial condition must ensure the condition of
    // hydrostatic equilibrium for pressure P = P0 - 0.1*\rho*y
	
    // Athena initial conditions are
    // if ( y > 0.0 ) {
    //   h_U(i,j,ID) = 2.0;
    // } else {
    //   h_U(i,j,ID) = 1.0;
    // }
    // h_U(i,j,IP) = P0 + gravity_x*x + gravity_y*y;
    // h_U(i,j,IU) = 0.0;
    // h_U(i,j,IV) = amplitude*(1+cos(2*M_PI*x))*(1+cos(0.5*M_PI*y))/4;
    
    if ( z > (zmin+zmax)/2 ) {
      Udata(i,j,k,ID) = d1;
    } else {
      Udata(i,j,k,ID) = d0;
    }
    Udata(i,j,k,IU) = 0.0;
    Udata(i,j,k,IV) = 0.0;
    // if (randomEnabled)
    //   Udata(i,j,IV) = amplitude * ( rand() * 1.0 / RAND_MAX - 0.5);
    // else
    Udata(i,j,k,IW) = amplitude * 
      (1+cos(2*M_PI*x/Lx))*
      (1+cos(2*M_PI*y/Ly))/4;

    Udata(i,j,k,IE) = (P0 + Udata(i,j,k,ID)*(gravity_x*x +
					     gravity_y*y +
					     gravity_z*z))/(gamma0-1);

    // init gravity field
    gravity(i,j,k,IX) = gravity_x;
    gravity(i,j,k,IY) = gravity_y;
    gravity(i,j,k,IZ) = gravity_z;
    
  } // end operator ()

  RayleighTaylorInstabilityParams rtiparams;
  DataArray3d Udata;
  VectorField3d gravity;
  
}; // class RayleighTaylorInstabilityFunctor3D

} // namespace  muscl

} // namespace euler_kokkos

#endif // HYDRO_INIT_FUNCTORS_3D_H_
