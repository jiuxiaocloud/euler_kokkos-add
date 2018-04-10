#ifndef MHD_INIT_FUNCTORS_3D_H_
#define MHD_INIT_FUNCTORS_3D_H_

#include <limits> // for std::numeric_limits
#ifdef __CUDA_ARCH__
#include <math_constants.h> // for cuda math constants, e.g. CUDART_INF
#endif // __CUDA_ARCH__

#include "shared/kokkos_shared.h"
#include "MHDBaseFunctor3D.h"

// some utils
#include "shared/utils.h"

// init conditions
#include "shared/BlastParams.h"
#include "muscl/OrszagTangInit.h"
#include "shared/RotorParams.h"

namespace euler_kokkos { namespace muscl {

/*************************************************/
/*************************************************/
/*************************************************/
class InitImplodeFunctor3D_MHD : public MHDBaseFunctor3D {

public:
  InitImplodeFunctor3D_MHD(HydroParams params,
			   DataArray3d Udata) :
    MHDBaseFunctor3D(params), Udata(Udata)  {};
  
  // static method which does it all: create and execute functor
  static void apply(HydroParams params,
                    DataArray3d Udata,
		    int         nbCells)
  {
    InitImplodeFunctor3D_MHD functor(params, Udata);
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
    
    real_t tmp = x+y+z;
    if (tmp > 0.5 && tmp < 2.5) {
      Udata(i,j,k , ID)  = 1.0;
      Udata(i,j,k , IU)  = 0.0;
      Udata(i,j,k , IV)  = 0.0;
      Udata(i,j,k , IW)  = 0.0;
      Udata(i,j,k , IBX) = 0.5;
      Udata(i,j,k , IBY) = 0.0;
      Udata(i,j,k , IBZ) = 0.0;
      Udata(i,j,k , IP)  = 1.0/(gamma0-1.0) +
	0.5* ( SQR(Udata(i,j,k , IBX)) +
	       SQR(Udata(i,j,k , IBY)) +
	       SQR(Udata(i,j,k , IBZ)) );
    } else {
      Udata(i,j,k , ID)  = 0.125;
      Udata(i,j,k , IU)  = 0.0;
      Udata(i,j,k , IV)  = 0.0;
      Udata(i,j,k , IW)  = 0.0;
      Udata(i,j,k , IBX) = 0.5;
      Udata(i,j,k , IBY) = 0.0;
      Udata(i,j,k , IBZ) = 0.0;
      Udata(i,j,k , IP)  = 0.14/(gamma0-1.0)  +
	0.5* ( SQR(Udata(i,j,k , IBX)) +
	       SQR(Udata(i,j,k , IBY)) +
	       SQR(Udata(i,j,k , IBZ)) );
    }
    
  } // end operator ()

  DataArray3d Udata;

}; // InitImplodeFunctor3D_MHD

/*************************************************/
/*************************************************/
/*************************************************/
class InitBlastFunctor3D_MHD : public MHDBaseFunctor3D {

public:
  InitBlastFunctor3D_MHD(HydroParams params,
			 BlastParams bParams,
			 DataArray3d Udata) :
    MHDBaseFunctor3D(params), bParams(bParams), Udata(Udata)  {};
  
  // static method which does it all: create and execute functor
  static void apply(HydroParams params,
		    BlastParams bParams,
                    DataArray3d Udata,
		    int         nbCells)
  {
    InitBlastFunctor3D_MHD functor(params, bParams, Udata);
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
      Udata(i,j,k , ID) = blast_density_in;
      Udata(i,j,k , IU) = 0.0;
      Udata(i,j,k , IV) = 0.0;
      Udata(i,j,k , IW) = 0.0;
      Udata(i,j,k , IA) = 0.5;
      Udata(i,j,k , IB) = 0.5;
      Udata(i,j,k , IC) = 0.5;
      Udata(i,j,k , IP) = blast_pressure_in/(gamma0-1.0) +
	0.5* ( SQR(Udata(i,j,k , IA)) +
	       SQR(Udata(i,j,k , IB)) +
	       SQR(Udata(i,j,k , IC)) );
    } else {
      Udata(i,j,k , ID) = blast_density_out;
      Udata(i,j,k , IU) = 0.0;
      Udata(i,j,k , IV) = 0.0;
      Udata(i,j,k , IW) = 0.0;
      Udata(i,j,k , IA) = 0.5;
      Udata(i,j,k , IB) = 0.5;
      Udata(i,j,k , IC) = 0.5;
      Udata(i,j,k , IP) = blast_pressure_out/(gamma0-1.0) +
	0.5* ( SQR(Udata(i,j,k , IA)) +
	       SQR(Udata(i,j,k , IB)) +
	       SQR(Udata(i,j,k , IC)) );
    }
    
  } // end operator ()
  
  BlastParams bParams;
  DataArray3d Udata;
  
}; // InitBlastFunctor3D_MHD

/*************************************************/
/*************************************************/
/*************************************************/
template<OrszagTang_init_type ot_type>
class InitOrszagTangFunctor3D : public MHDBaseFunctor3D {

public:
  InitOrszagTangFunctor3D(HydroParams params,
			  DataArray3d Udata) :
    MHDBaseFunctor3D(params), Udata(Udata)  {};
  
  // static method which does it all: create and execute functor
  static void apply(HydroParams params,
                    DataArray3d Udata,
		    int         nbCells)
  {
    InitOrszagTangFunctor3D<ot_type> functor(params, Udata);
    Kokkos::parallel_for(nbCells, functor);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& index) const
  {

    if (ot_type == INIT_ALL_VAR_BUT_ENERGY)
      init_all_var_but_energy(index);
    else if(ot_type == INIT_ENERGY)
      init_energy(index);

  } // end operator ()

  KOKKOS_INLINE_FUNCTION
  void init_all_var_but_energy(const int index) const
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
    UNUSED(k_mpi);

    const int nx = params.nx;
    const int ny = params.ny;
    const int nz = params.nz;
    UNUSED(nz);
    
    const real_t xmin = params.xmin;
    const real_t ymin = params.ymin;
    const real_t zmin = params.zmin;
    UNUSED(zmin);
        
    const double dx = params.dx;
    const double dy = params.dy;
    const double dz = params.dz;
    UNUSED(dz);
    
    const real_t gamma0 = params.settings.gamma0;
    
    const double TwoPi = 4.0*asin(1.0);
    const double B0    = 1.0/sqrt(2.0*TwoPi);
    const double p0    = gamma0/(2.0*TwoPi);
    const double d0    = gamma0*p0;
    const double v0    = 1.0;

    int i,j,k;
    index2coord(index,i,j,k,isize,jsize,ksize);
    
    double xPos = xmin + dx/2 + (i+nx*i_mpi-ghostWidth)*dx;
    double yPos = ymin + dy/2 + (j+ny*j_mpi-ghostWidth)*dy;
    //double zPos = zmin + dz/2 + (k+nz*k_mpi-ghostWidth)*dz;
        
    // density
    Udata(i,j,k,ID) = d0;
    
    // rho*vx
    Udata(i,j,k,IU)  = static_cast<real_t>(-d0*v0*sin(yPos*TwoPi));
    
    // rho*vy
    Udata(i,j,k,IV)  = static_cast<real_t>( d0*v0*sin(xPos*TwoPi));
    
    // rho*vz
    Udata(i,j,k,IW) =  ZERO_F;

    // bx, by, bz
    Udata(i,j,k, IBX) = -B0*sin(    yPos*TwoPi);
    Udata(i,j,k, IBY) =  B0*sin(2.0*xPos*TwoPi);
    Udata(i,j,k, IBZ) =  0.0;

  } // init_all_var_but_energy

  KOKKOS_INLINE_FUNCTION
  void init_energy(const int index) const
  {

    const int isize = params.isize;
    const int jsize = params.jsize;
    const int ksize = params.ksize;
    const int ghostWidth = params.ghostWidth;
    
    //const double xmin = params.xmin;
    //const double ymin = params.ymin;
    //const double zmin = params.zmin;
    
    //const double dx = params.dx;
    //const double dy = params.dy;
    //const double dz = params.dz;
    
    const real_t gamma0 = params.settings.gamma0;
    
    const double TwoPi = 4.0*asin(1.0);
    //const double B0    = 1.0/sqrt(2.0*TwoPi);
    const double p0    = gamma0/(2.0*TwoPi);
    //const double d0    = gamma0*p0;
    //const double v0    = 1.0;

    int i,j,k;
    index2coord(index,i,j,k,isize,jsize,ksize);
    
    //double xPos = xmin + dx/2 + (i-ghostWidth)*dx;
    //double yPos = ymin + dy/2 + (j-ghostWidth)*dy;
    //double zPos = zmin + dz/2 + (k-ghostWidth)*dz;
        
    if (i<isize-1 and j<jsize-1) {

      Udata(i,j,k,IP)  = p0 / (gamma0-1.0) +
	0.5 * ( SQR(Udata(i,j,k,IU)) / Udata(i,j,k,ID) +
		SQR(Udata(i,j,k,IV)) / Udata(i,j,k,ID) +
		0.25*SQR(Udata(i,j,k,IBX) + Udata(i+1,j  ,k  ,IBX)) + 
		0.25*SQR(Udata(i,j,k,IBY) + Udata(i  ,j+1,k  ,IBY)) );

    } else if ( (i <isize-1) and (j==jsize-1)) {

      Udata(i,j,k,IP)  = p0 / (gamma0-1.0) +
	0.5 * ( SQR(Udata(i,j,k,IU)) / Udata(i,j,k,ID) +
		SQR(Udata(i,j,k,IV)) / Udata(i,j,k,ID) +
		0.25*SQR(Udata(i,j,k,IBX) + Udata(i+1,j           ,k,IBX)) + 
		0.25*SQR(Udata(i,j,k,IBY) + Udata(i  ,2*ghostWidth,k,IBY)) );

    } else if ( (i==isize-1) and (j <jsize-1)) {

      Udata(i,j,k,IP)  = p0 / (gamma0-1.0) +
	0.5 * ( SQR(Udata(i,j,k,IU)) / Udata(i,j,k,ID) +
		SQR(Udata(i,j,k,IV)) / Udata(i,j,k,ID) +
		0.25*SQR(Udata(i,j,k,IBX) + Udata(2*ghostWidth,j  ,k  ,IBX)) + 
		0.25*SQR(Udata(i,j,k,IBY) + Udata(i           ,j+1,k  ,IBY)) );

    } else if ( (i==isize-1) and (j==jsize-1) ) {

      Udata(i,j,k,IP)  = p0 / (gamma0-1.0) +
	0.5 * ( SQR(Udata(i,j,k,IU)) / Udata(i,j,k,ID) +
		SQR(Udata(i,j,k,IV)) / Udata(i,j,k,ID) +
		0.25*SQR(Udata(i,j,k,IBX) + Udata(2*ghostWidth,j,k ,IBX)) + 
		0.25*SQR(Udata(i,j,k,IBY) + Udata(i,2*ghostWidth,k ,IBY)) );

    }
    
  } // init_energy
  
  DataArray3d Udata;
  
}; // InitOrszagTangFunctor3D

/*************************************************/
/*************************************************/
/*************************************************/
/**
 * The two-dimensional MHD rotor problem. Some references:
 * - Balsara and Spicer, 1999, JCP, 149, 270.
 * - G. Toth, "The div(B)=0 constraint in shock-capturing MHD codes",
 *   JCP, 161, 605 (2000)
 * 
 * Initial conditions are taken from Toth's paper.
 *
 */
class InitRotorFunctor3D_MHD : public MHDBaseFunctor3D {

public:
  InitRotorFunctor3D_MHD(HydroParams params,
			 RotorParams rParams,
			 DataArray3d Udata) :
    MHDBaseFunctor3D(params), rParams(rParams), Udata(Udata)  {};
  
  // static method which does it all: create and execute functor
  static void apply(HydroParams params,
		    RotorParams rParams,
                    DataArray3d Udata,
		    int         nbCells)
  {
    InitRotorFunctor3D_MHD functor(params, rParams, Udata);
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
    UNUSED(k_mpi);
    
    const int nx = params.nx;
    const int ny = params.ny;
    const int nz = params.nz;
    UNUSED(nz);

    const real_t xmin = params.xmin;
    const real_t ymin = params.ymin;
    const real_t xmax = params.xmax;
    const real_t ymax = params.ymax;

    const real_t dx = params.dx;
    const real_t dy = params.dy;
    
    const real_t gamma0 = params.settings.gamma0;

    // rotor problem parameters
    const real_t r0      = rParams.r0;
    const real_t r1      = rParams.r1;
    const real_t u0      = rParams.u0;
    const real_t p0      = rParams.p0;
    const real_t b0      = rParams.b0;

    const real_t xCenter = (xmax + xmin)/2;
    const real_t yCenter = (ymax + ymin)/2;

    int i,j,k;
    index2coord(index,i,j,k,isize,jsize,ksize);
    
    real_t x = xmin + dx/2 + (i+nx*i_mpi-ghostWidth)*dx;
    real_t y = ymin + dy/2 + (j+ny*j_mpi-ghostWidth)*dy;
    
    real_t r = SQRT( (x-xCenter)*(x-xCenter) +
		     (y-yCenter)*(y-yCenter) );
    real_t f_r = (r1-r)/(r1-r0);

    if (r<=r0) {
      Udata(i,j,k,ID) = 10.0;
      Udata(i,j,k,IU) = -u0*(y-yCenter)/r0;
      Udata(i,j,k,IV) =  u0*(x-xCenter)/r0;
    } else if (r<=r1) {
      Udata(i,j,k,ID) = 1+9*f_r;
      Udata(i,j,k,IU) = -f_r*u0*(y-yCenter)/r;
      Udata(i,j,k,IV) =  f_r*u0*(x-xCenter)/r;
    } else {
      Udata(i,j,k,ID) = 1.0;
      Udata(i,j,k,IU) = 0.0;
      Udata(i,j,k,IV) = 0.0;
    }

    Udata(i,j,k,IW) = 0.0;
    Udata(i,j,k,IA) = b0; //5.0/SQRT(FourPi);
    Udata(i,j,k,IB) = 0.0;
    Udata(i,j,k,IC) = 0.0;
    Udata(i,j,k,IP) = p0/(gamma0-1.0) + 
      ( Udata(i,j,k,IU)*Udata(i,j,k,IU) + 
	Udata(i,j,k,IV)*Udata(i,j,k,IV) +
	Udata(i,j,k,IW)*Udata(i,j,k,IW) )/2/Udata(i,j,k,ID) +
      ( Udata(i,j,k,IA)*Udata(i,j,k,IA) )/2;
    
    
  } // end operator ()
  
  RotorParams rParams;
  DataArray3d Udata;
  
}; // InitRotorFunctor3D_MHD

} // namespace muscl

} // namespace euler_kokkos

#endif // MHD_INIT_FUNCTORS_3D_H_
