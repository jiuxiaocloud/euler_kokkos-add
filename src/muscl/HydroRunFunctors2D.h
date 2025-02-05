#ifndef HYDRO_RUN_FUNCTORS_2D_H_
#define HYDRO_RUN_FUNCTORS_2D_H_

#include <limits> // for std::numeric_limits
#ifdef __CUDA_ARCH__
#include <math_constants.h> // for cuda math constants, e.g. CUDART_INF
#endif                      // __CUDA_ARCH__

#include "shared/kokkos_shared.h"
#include "HydroBaseFunctor2D.h"
#include "shared/RiemannSolvers.h"

namespace euler_kokkos
{
  namespace muscl
  {

    /*************************************************/
    /*************************************************/
    /*************************************************/
    class ComputeDtFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Compute time step satisfying CFL constraint.
       *
       * \param[in] params
       * \param[in] Udata
       */
      ComputeDtFunctor2D(HydroParams params,
                         DataArray2d Udata) : HydroBaseFunctor2D(params),
                                              Udata(Udata){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        DataArray2d Udata,
                        int nbCells,
                        real_t &invDt)
      {
        ComputeDtFunctor2D functor(params, Udata);
        Kokkos::parallel_reduce(nbCells, functor, invDt);
      }

      // Tell each thread how to initialize its reduction result.
      KOKKOS_INLINE_FUNCTION
      void init(real_t &dst) const
      {
        // The identity under max is -Inf.
        // Kokkos does not come with a portable way to access
        // floating-point Inf and NaN.
#ifdef __CUDA_ARCH__
        dst = -CUDART_INF;
#else
        dst = std::numeric_limits<real_t>::min();
#endif  // __CUDA_ARCH__
      } // init

      /* this is a reduce (max) functor */
      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index, real_t &invDt) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;
        // const int nbvar = params.nbvar;
        const real_t dx = params.dx;
        const real_t dy = params.dy;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth && j < jsize - ghostWidth &&
            i >= ghostWidth && i < isize - ghostWidth)
        {

          HydroState uLoc; // conservative    variables in current cell
          HydroState qLoc; // primitive    variables in current cell
          real_t c = 0.0;
          real_t vx, vy;

          // get local conservative variable
          uLoc[ID] = Udata(i, j, ID);
          uLoc[IP] = Udata(i, j, IP);
          uLoc[IU] = Udata(i, j, IU);
          uLoc[IV] = Udata(i, j, IV);
#if COP
          uLoc[IDY1] = Udata(i, j, IDY1);
#endif
          // get primitive variables in current cell
          computePrimitives(uLoc, &c, qLoc);
// NOTE:时间组分计算
#if COP
          real_t T = 298.15;
          get_CopC2(T, params, qLoc, &c);
          c = sqrt(c);
#endif
          vx = c + FABS(qLoc[IU]);
          vy = c + FABS(qLoc[IV]);
#if COP
#if ADD_VIS
          // NOTE:添加粘性时时间要多一个对比，使时间步变小，以避免算爆，算爆一般是时间间步太大导致
          // real_t T = qLoc[IP] / qLoc[ID] / Ru / ((1 - y1) / params.species(0, 6) + y1 / params.species(1, 6)); //单位K
          HydroSpecies Yi = {1 - qLoc[IDY1], qLoc[IDY1]};
          real_t miu_ = get_CopMiu(params.species, Yi, T); //单位Pa.s
          invDt = FMAX(invDt, 14.0 / 3 * miu_ / (qLoc[ID] * FMIN(dx * dx, dy * dy)));
#endif
#endif
          invDt = FMAX(invDt, vx / dx + vy / dy);
        }

      } // operator ()

      // "Join" intermediate results from different threads.
      // This should normally implement the same reduction
      // operation as operator() above. Note that both input
      // arguments MUST be declared volatile.
      KOKKOS_INLINE_FUNCTION
      void join(volatile real_t &dst,
                const volatile real_t &src) const
      {
        // max reduce
        if (dst < src)
        {
          dst = src;
        }
      } // join

      DataArray2d Udata;

    }; // ComputeDtFunctor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    /**
     * A specialized functor to compute CFL dt constraint when gravity source
     * term is activated.
     */
    class ComputeDtGravityFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Compute time step satisfying CFL constraint.
       *
       * \param[in] params
       * \param[in] Udata
       */
      // #if COP
      //       ComputeDtGravityFunctor2D(HydroParams params,
      //                                 real_t cfl,
      //                                 VectorField2d gravity,
      //                                 DataArray2d Udata) : HydroBaseFunctor2D(params),
      //                                                      cfl(cfl),
      //                                                      gravity(gravity),
      //                                                      Udata(Udata){};

      //       // static method which does it all: create and execute functor
      //       static void apply(HydroParams params,
      //                         real_t cfl,
      //                         VectorField2d gravity,
      //                         DataArray2d Udata,
      //                         int nbCells,
      //                         real_t &invDt)
      //       {
      //         ComputeDtGravityFunctor2D functor(params, cfl, gravity, Udata);
      //         Kokkos::parallel_reduce(nbCells, functor, invDt);
      //       }
      // #else
      ComputeDtGravityFunctor2D(HydroParams params,
                                real_t cfl,
                                VectorField2d gravity,
                                DataArray2d Udata) : HydroBaseFunctor2D(params),
                                                     cfl(cfl),
                                                     gravity(gravity),
                                                     Udata(Udata){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        real_t cfl,
                        VectorField2d gravity,
                        DataArray2d Udata,
                        int nbCells,
                        real_t &invDt)
      {
        ComputeDtGravityFunctor2D functor(params, cfl, gravity, Udata);
        Kokkos::parallel_reduce(nbCells, functor, invDt);
      }
      // #endif
      // Tell each thread how to initialize its reduction result.
      KOKKOS_INLINE_FUNCTION
      void init(real_t &dst) const
      {
        // The identity under max is -Inf.
        // Kokkos does not come with a portable way to access
        // floating-point Inf and NaN.
#ifdef __CUDA_ARCH__
        dst = -CUDART_INF;
#else
        dst = std::numeric_limits<real_t>::min();
#endif  // __CUDA_ARCH__
      } // init

      /* this is a reduce (max) functor */
      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index, real_t &invDt) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;
        // const int nbvar = params.nbvar;
        const real_t dx = fmin(params.dx, params.dy);

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth && j < jsize - ghostWidth &&
            i >= ghostWidth && i < isize - ghostWidth)
        {

          HydroState uLoc; // conservative    variables in current cell
          HydroState qLoc; // primitive    variables in current cell
          real_t c = 0.0;

          // get local conservative variable
          uLoc[ID] = Udata(i, j, ID);
          uLoc[IP] = Udata(i, j, IP);
          uLoc[IU] = Udata(i, j, IU);
          uLoc[IV] = Udata(i, j, IV);

          // get primitive variables in current cell
          computePrimitives(uLoc, &c, qLoc);
          real_t velocity = 0.0;
          velocity += c + FABS(qLoc[IU]);
          velocity += c + FABS(qLoc[IV]);

          /* Due to the gravitational acceleration, the CFL condition
           * can be written as
           * g dt^2 / (2 dx) + u dt / dx <= cfl
           * where u = sum(|v_i| + c_s) and g = sum(|g_i|)
           *
           * u / dx has to be corrected by a factor k / (sqrt(1 + 2k) - 1)
           * in order to satisfy the new CFL, where k = g dx cfl / u^2
           */
          double k = fabs(gravity(i, j, IX)) + fabs(gravity(i, j, IY));

          k *= cfl * dx / (velocity * velocity);

          /* prevent numerical errors due to very low gravity */
          k = fmax(k, 1e-4);

          velocity *= k / (sqrt(1.0 + 2.0 * k) - 1.0);

          invDt = fmax(invDt, velocity / dx);
        }

      } // operator ()

      // "Join" intermediate results from different threads.
      // This should normally implement the same reduction
      // operation as operator() above. Note that both input
      // arguments MUST be declared volatile.
      KOKKOS_INLINE_FUNCTION
      void join(volatile real_t &dst,
                const volatile real_t &src) const
      {
        // max reduce
        if (dst < src)
        {
          dst = src;
        }
      } // join

      real_t cfl;
      VectorField2d gravity;
      DataArray2d Udata;

    }; // ComputeDtGravityFunctor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    class ConvertToPrimitivesFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Convert conservative variables to primitive ones using equation of state.
       *
       * \param[in] params
       * \param[in] Udata conservative variables
       * \param[out] Qdata primitive variables
       */
      ConvertToPrimitivesFunctor2D(HydroParams params,
                                   DataArray2d Udata,
                                   DataArray2d Qdata) : HydroBaseFunctor2D(params), Udata(Udata), Qdata(Qdata){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        DataArray2d Udata,
                        DataArray2d Qdata)
      {
        int nbCells = params.isize * params.jsize;
        ConvertToPrimitivesFunctor2D functor(params, Udata, Qdata);
        Kokkos::parallel_for(nbCells, functor);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        // const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= 0 && j < jsize &&
            i >= 0 && i < isize)
        {

          HydroState uLoc; // conservative    variables in current cell
          HydroState qLoc; // primitive    variables in current cell
          real_t c;

          // get local conservative variable
          uLoc[ID] = Udata(i, j, ID);
          uLoc[IP] = Udata(i, j, IP);
          uLoc[IU] = Udata(i, j, IU);
          uLoc[IV] = Udata(i, j, IV);
#if COP
          uLoc[IDY1] = Udata(i, j, IDY1);
#endif
          // get primitive variables in current cell
          computePrimitives(uLoc, &c, qLoc);

          // copy q state in q global
          Qdata(i, j, ID) = qLoc[ID];
          Qdata(i, j, IP) = qLoc[IP];
          Qdata(i, j, IU) = qLoc[IU];
          Qdata(i, j, IV) = qLoc[IV];
#if COP
          Qdata(i, j, IDY1) = qLoc[IDY1];
#endif
        }
      }

      DataArray2d Udata;
      DataArray2d Qdata;

    }; // ConvertToPrimitivesFunctor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    class ComputeFluxesAndUpdateFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Perform time update by computing Riemann fluxes at cell interfaces.
       *
       * \param[in] params
       * \param[in,out] Udata conservative variables
       * \param[in] Qm_x primitive variables reconstructed on face -X
       * \param[in] Qm_y primitive variables reconstructed on face -Y
       * \param[in] Qp_x primitive variables reconstructed on face +X
       * \param[in] Qp_y primitive variables reconstructed on face +Y
       */
      ComputeFluxesAndUpdateFunctor2D(HydroParams params,
                                      DataArray2d Udata,
                                      DataArray2d Qm_x,
                                      DataArray2d Qm_y,
                                      DataArray2d Qp_x,
                                      DataArray2d Qp_y,
                                      real_t dtdx,
                                      real_t dtdy) : HydroBaseFunctor2D(params), Udata(Udata),
                                                     Qm_x(Qm_x), Qm_y(Qm_y), Qp_x(Qp_x), Qp_y(Qp_y),
                                                     dtdx(dtdx), dtdy(dtdy){};

      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth && j <= jsize - ghostWidth &&
            i >= ghostWidth && i <= isize - ghostWidth)
        {

          HydroState qleft, qright;
          HydroState flux_x, flux_y;
          HydroState qgdnv;

          //
          // Solve Riemann problem at X-interfaces and compute
          // X-fluxes
          //
          qleft[ID] = Qm_x(i - 1, j, ID);
          qleft[IP] = Qm_x(i - 1, j, IP);
          qleft[IU] = Qm_x(i - 1, j, IU);
          qleft[IV] = Qm_x(i - 1, j, IV);

          qright[ID] = Qp_x(i, j, ID);
          qright[IP] = Qp_x(i, j, IP);
          qright[IU] = Qp_x(i, j, IU);
          qright[IV] = Qp_x(i, j, IV);

          // compute hydro flux_x
          // riemann_hllc(qleft,qright,qgdnv,flux_x);
          riemann_hydro(qleft, qright, qgdnv, flux_x, params);

          //
          // Solve Riemann problem at Y-interfaces and compute Y-fluxes
          //
          qleft[ID] = Qm_y(i, j - 1, ID);
          qleft[IP] = Qm_y(i, j - 1, IP);
          qleft[IU] = Qm_y(i, j - 1, IV); // watchout IU, IV permutation
          qleft[IV] = Qm_y(i, j - 1, IU); // watchout IU, IV permutation

          qright[ID] = Qp_y(i, j, ID);
          qright[IP] = Qp_y(i, j, IP);
          qright[IU] = Qp_y(i, j, IV); // watchout IU, IV permutation
          qright[IV] = Qp_y(i, j, IU); // watchout IU, IV permutation

          // compute hydro flux_y
          // riemann_hllc(qleft,qright,qgdnv,flux_y);
          riemann_hydro(qleft, qright, qgdnv, flux_y, params);

          //
          // update hydro array
          //
          Udata(i - 1, j, ID) += -flux_x[ID] * dtdx;
          Udata(i - 1, j, IP) += -flux_x[IP] * dtdx;
          Udata(i - 1, j, IU) += -flux_x[IU] * dtdx;
          Udata(i - 1, j, IV) += -flux_x[IV] * dtdx;

          Udata(i, j, ID) += flux_x[ID] * dtdx;
          Udata(i, j, IP) += flux_x[IP] * dtdx;
          Udata(i, j, IU) += flux_x[IU] * dtdx;
          Udata(i, j, IV) += flux_x[IV] * dtdx;

          Udata(i, j - 1, ID) += -flux_y[ID] * dtdy;
          Udata(i, j - 1, IP) += -flux_y[IP] * dtdy;
          Udata(i, j - 1, IU) += -flux_y[IV] * dtdy; // watchout IU and IV swapped
          Udata(i, j - 1, IV) += -flux_y[IU] * dtdy; // watchout IU and IV swapped

          Udata(i, j, ID) += flux_y[ID] * dtdy;
          Udata(i, j, IP) += flux_y[IP] * dtdy;
          Udata(i, j, IU) += flux_y[IV] * dtdy; // watchout IU and IV swapped
          Udata(i, j, IV) += flux_y[IU] * dtdy; // watchout IU and IV swapped
        }
      }

      DataArray2d Udata;
      DataArray2d Qm_x, Qm_y, Qp_x, Qp_y;
      real_t dtdx, dtdy;

    }; // ComputeFluxesAndUpdateFunctor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    class ComputeTraceFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Compute (slope extrapolated) reconstructed states at face centers for all faces.
       *
       * \param[in] params
       * \param[in] Qdata primitive variables
       * \param[out] Qm_x primitive variables reconstructed at center face -X
       * \param[out] Qm_y primitive variables reconstructed at center face -Y
       * \param[out] Qp_x primitive variables reconstructed at center face +X
       * \param[out] Qp_y primitive variables reconstructed at center face +Y
       */
      ComputeTraceFunctor2D(HydroParams params,
                            DataArray2d Qdata,
                            DataArray2d Qm_x,
                            DataArray2d Qm_y,
                            DataArray2d Qp_x,
                            DataArray2d Qp_y,
                            real_t dtdx,
                            real_t dtdy) : HydroBaseFunctor2D(params),
                                           Qdata(Qdata),
                                           Qm_x(Qm_x), Qm_y(Qm_y), Qp_x(Qp_x), Qp_y(Qp_y),
                                           dtdx(dtdx), dtdy(dtdy){};

      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= 1 && j <= jsize - ghostWidth &&
            i >= 1 && i <= isize - ghostWidth)
        {

          HydroState qLoc; // local primitive variables
          HydroState qPlusX;
          HydroState qMinusX;
          HydroState qPlusY;
          HydroState qMinusY;

          HydroState dqX;
          HydroState dqY;

          HydroState qmX;
          HydroState qmY;
          HydroState qpX;
          HydroState qpY;

          // get primitive variables state vector
          {
            qLoc[ID] = Qdata(i, j, ID);
            qPlusX[ID] = Qdata(i + 1, j, ID);
            qMinusX[ID] = Qdata(i - 1, j, ID);
            qPlusY[ID] = Qdata(i, j + 1, ID);
            qMinusY[ID] = Qdata(i, j - 1, ID);

            qLoc[IP] = Qdata(i, j, IP);
            qPlusX[IP] = Qdata(i + 1, j, IP);
            qMinusX[IP] = Qdata(i - 1, j, IP);
            qPlusY[IP] = Qdata(i, j + 1, IP);
            qMinusY[IP] = Qdata(i, j - 1, IP);

            qLoc[IU] = Qdata(i, j, IU);
            qPlusX[IU] = Qdata(i + 1, j, IU);
            qMinusX[IU] = Qdata(i - 1, j, IU);
            qPlusY[IU] = Qdata(i, j + 1, IU);
            qMinusY[IU] = Qdata(i, j - 1, IU);

            qLoc[IV] = Qdata(i, j, IV);
            qPlusX[IV] = Qdata(i + 1, j, IV);
            qMinusX[IV] = Qdata(i - 1, j, IV);
            qPlusY[IV] = Qdata(i, j + 1, IV);
            qMinusY[IV] = Qdata(i, j - 1, IV);

          } //

          // get hydro slopes dq
          slope_unsplit_hydro_2d(qLoc,
                                 qPlusX, qMinusX,
                                 qPlusY, qMinusY,
                                 dqX, dqY);

          // compute qm, qp
          trace_unsplit_hydro_2d(qLoc,
                                 dqX, dqY,
                                 dtdx, dtdy,
                                 qmX, qmY,
                                 qpX, qpY);

          // store qm, qp : only what is really needed
          Qm_x(i, j, ID) = qmX[ID];
          Qp_x(i, j, ID) = qpX[ID];
          Qm_y(i, j, ID) = qmY[ID];
          Qp_y(i, j, ID) = qpY[ID];

          Qm_x(i, j, IP) = qmX[IP];
          Qp_x(i, j, IP) = qpX[IP];
          Qm_y(i, j, IP) = qmY[IP];
          Qp_y(i, j, IP) = qpY[IP];

          Qm_x(i, j, IU) = qmX[IU];
          Qp_x(i, j, IU) = qpX[IU];
          Qm_y(i, j, IU) = qmY[IU];
          Qp_y(i, j, IU) = qpY[IU];

          Qm_x(i, j, IV) = qmX[IV];
          Qp_x(i, j, IV) = qpX[IV];
          Qm_y(i, j, IV) = qmY[IV];
          Qp_y(i, j, IV) = qpY[IV];
        }
      }

      DataArray2d Qdata;
      DataArray2d Qm_x, Qm_y, Qp_x, Qp_y;
      real_t dtdx, dtdy;

    }; // ComputeTraceFunctor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    class ComputeAndStoreFluxesFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Compute (all-in-one) reconstructed states on faces, then compute Riemann fluxes and store them.
       *
       * \note All-in-one here means the stencil of this operator is larger (need to
       * fetch data in neighbor of neighbor).
       *
       * \param[in] Qdata primitive variables (at cell center)
       * \param[out] FluxData_x flux coming from the left neighbor along X
       * \param[out] FluxData_y flux coming from the left neighbor along Y
       * \param[in] gravity_enabled boolean value to activate static gravity
       * \param[in] gravity is a vector field
       */
      // #if COP
      //       ComputeAndStoreFluxesFunctor2D(HydroParams params,
      //                                      Add_Params add_params,
      //                                      DataArray2d Qdata,
      //                                      DataArray2d FluxData_x,
      //                                      DataArray2d FluxData_y,
      //                                      real_t dt,
      //                                      bool gravity_enabled,
      //                                      VectorField2d gravity) : HydroBaseFunctor2D(params, add_params),
      //                                                               Qdata(Qdata),
      //                                                               FluxData_x(FluxData_x),
      //                                                               FluxData_y(FluxData_y),
      //                                                               dt(dt),
      //                                                               dtdx(dt / params.dx),
      //                                                               dtdy(dt / params.dy),
      //                                                               gravity_enabled(gravity_enabled),
      //                                                               gravity(gravity){};

      //       // static method which does it all: create and execute functor
      //       static void apply(HydroParams params,
      //                         Add_Params add_params,
      //                         DataArray2d Qdata,
      //                         DataArray2d FluxData_x,
      //                         DataArray2d FluxData_y,
      //                         real_t dt,
      //                         bool gravity_enabled,
      //                         VectorField2d gravity)
      //       {
      //         int nbCells = params.isize * params.jsize;
      //         ComputeAndStoreFluxesFunctor2D functor(params, add_params, Qdata,
      //                                                FluxData_x, FluxData_y,
      //                                                dt,
      //                                                gravity_enabled,
      //                                                gravity);
      //         Kokkos::parallel_for(nbCells, functor);
      //       }
      // #else
      ComputeAndStoreFluxesFunctor2D(HydroParams params,
                                     DataArray2d Qdata,
                                     DataArray2d FluxData_x,
                                     DataArray2d FluxData_y,
                                     real_t dt,
                                     bool gravity_enabled,
                                     VectorField2d gravity) : HydroBaseFunctor2D(params),
                                                              Qdata(Qdata),
                                                              FluxData_x(FluxData_x),
                                                              FluxData_y(FluxData_y),
                                                              dt(dt),
                                                              dtdx(dt / params.dx),
                                                              dtdy(dt / params.dy),
                                                              gravity_enabled(gravity_enabled),
                                                              gravity(gravity){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        DataArray2d Qdata,
                        DataArray2d FluxData_x,
                        DataArray2d FluxData_y,
                        real_t dt,
                        bool gravity_enabled,
                        VectorField2d gravity)
      {
        int nbCells = params.isize * params.jsize;
        ComputeAndStoreFluxesFunctor2D functor(params, Qdata,
                                               FluxData_x, FluxData_y,
                                               dt,
                                               gravity_enabled,
                                               gravity);
        Kokkos::parallel_for(nbCells, functor);
      }
      // #endif

      KOKKOS_INLINE_FUNCTION void
      operator()(const int &index) const
      { // NOTE： 调用的计算和储存通量核函数
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth && j <= jsize - ghostWidth &&
            i >= ghostWidth && i <= isize - ghostWidth)
        {

          // local primitive variables
          HydroState qLoc; // local primitive variables

          // local primitive variables in neighbor cell
          HydroState qLocNeighbor;

          // local primitive variables in neighborbood
          HydroState qNeighbors_0;
          HydroState qNeighbors_1;
          HydroState qNeighbors_2;
          HydroState qNeighbors_3;

          // Local slopes and neighbor slopes
          HydroState dqX;
          HydroState dqY;
          HydroState dqX_neighbor;
          HydroState dqY_neighbor;

          // Local variables for Riemann problems solving
          HydroState qleft;
          HydroState qright;
          HydroState qgdnv;
          HydroState flux_x;
          HydroState flux_y;

          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          // deal with left interface along X !
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          // get primitive variables state vector
          qLoc[ID] = Qdata(i, j, ID);
          qNeighbors_0[ID] = Qdata(i + 1, j, ID);
          qNeighbors_1[ID] = Qdata(i - 1, j, ID);
          qNeighbors_2[ID] = Qdata(i, j + 1, ID);
          qNeighbors_3[ID] = Qdata(i, j - 1, ID);

          qLoc[IP] = Qdata(i, j, IP);
          qNeighbors_0[IP] = Qdata(i + 1, j, IP);
          qNeighbors_1[IP] = Qdata(i - 1, j, IP);
          qNeighbors_2[IP] = Qdata(i, j + 1, IP);
          qNeighbors_3[IP] = Qdata(i, j - 1, IP);

          qLoc[IU] = Qdata(i, j, IU);
          qNeighbors_0[IU] = Qdata(i + 1, j, IU);
          qNeighbors_1[IU] = Qdata(i - 1, j, IU);
          qNeighbors_2[IU] = Qdata(i, j + 1, IU);
          qNeighbors_3[IU] = Qdata(i, j - 1, IU);

          qLoc[IV] = Qdata(i, j, IV);
          qNeighbors_0[IV] = Qdata(i + 1, j, IV);
          qNeighbors_1[IV] = Qdata(i - 1, j, IV);
          qNeighbors_2[IV] = Qdata(i, j + 1, IV);
          qNeighbors_3[IV] = Qdata(i, j - 1, IV);

#if COP
          qLoc[IDY1] = Qdata(i, j, IDY1);
          qNeighbors_0[IDY1] = Qdata(i + 1, j, IDY1);
          qNeighbors_1[IDY1] = Qdata(i - 1, j, IDY1);
          qNeighbors_2[IDY1] = Qdata(i, j + 1, IDY1);
          qNeighbors_3[IDY1] = Qdata(i, j - 1, IDY1);
#endif

          slope_unsplit_hydro_2d(qLoc,
                                 qNeighbors_0, qNeighbors_1,
                                 qNeighbors_2, qNeighbors_3,
                                 dqX, dqY);
          // NOTE：这个函数计算了（i，j）点的ux=dqX[IU],uy=dqX[IV],vx=dqY[IU],vy=dqY[IV]
          // slopes at left neighbor along X
          qLocNeighbor[ID] = Qdata(i - 1, j, ID);
          qNeighbors_0[ID] = Qdata(i, j, ID);
          qNeighbors_1[ID] = Qdata(i - 2, j, ID);
          qNeighbors_2[ID] = Qdata(i - 1, j + 1, ID);
          qNeighbors_3[ID] = Qdata(i - 1, j - 1, ID);

          qLocNeighbor[IP] = Qdata(i - 1, j, IP);
          qNeighbors_0[IP] = Qdata(i, j, IP);
          qNeighbors_1[IP] = Qdata(i - 2, j, IP);
          qNeighbors_2[IP] = Qdata(i - 1, j + 1, IP);
          qNeighbors_3[IP] = Qdata(i - 1, j - 1, IP);

          qLocNeighbor[IU] = Qdata(i - 1, j, IU);
          qNeighbors_0[IU] = Qdata(i, j, IU);
          qNeighbors_1[IU] = Qdata(i - 2, j, IU);
          qNeighbors_2[IU] = Qdata(i - 1, j + 1, IU);
          qNeighbors_3[IU] = Qdata(i - 1, j - 1, IU);

          qLocNeighbor[IV] = Qdata(i - 1, j, IV);
          qNeighbors_0[IV] = Qdata(i, j, IV);
          qNeighbors_1[IV] = Qdata(i - 2, j, IV);
          qNeighbors_2[IV] = Qdata(i - 1, j + 1, IV);
          qNeighbors_3[IV] = Qdata(i - 1, j - 1, IV);

#if COP
          qLocNeighbor[IDY1] = Qdata(i - 1, j, IDY1);
          qNeighbors_0[IDY1] = Qdata(i, j, IDY1);
          qNeighbors_1[IDY1] = Qdata(i - 2, j, IDY1);
          qNeighbors_2[IDY1] = Qdata(i - 1, j + 1, IDY1);
          qNeighbors_3[IDY1] = Qdata(i - 1, j - 1, IDY1);
#endif
          // NOTE:使用一阶差分算出了斜率并进行了斜率限制，//FAQ：可以用来构造通量但是能用来构造粘性项？
          slope_unsplit_hydro_2d(qLocNeighbor,
                                 qNeighbors_0, qNeighbors_1,
                                 qNeighbors_2, qNeighbors_3,
                                 dqX_neighbor, dqY_neighbor);

          //
          // compute reconstructed states at left interface along X
          //

          // left interface : right state
          trace_unsplit_2d_along_dir(qLoc,
                                     dqX, dqY,
                                     dtdx, dtdy, FACE_XMIN, qright);

          // left interface : left state
          trace_unsplit_2d_along_dir(qLocNeighbor,
                                     dqX_neighbor, dqY_neighbor,
                                     dtdx, dtdy, FACE_XMAX, qleft);

          if (gravity_enabled)
          {
            // we need to modify input to flux computation with
            // gravity predictor (half time step)

            qleft[IU] += 0.5 * dt * gravity(i - 1, j, IX);
            qleft[IV] += 0.5 * dt * gravity(i - 1, j, IY);

            qright[IU] += 0.5 * dt * gravity(i, j, IX);
            qright[IV] += 0.5 * dt * gravity(i, j, IY);
          }

          // Solve Riemann problem at X-interfaces and compute X-fluxes
          // riemann_2d(qleft,qright,qgdnv,flux_x);

#if COP
          // NOTE:flux_x 计算函数调用,在这个函数中添加Roe求解器
          riemann_hydro(qleft, qright, qgdnv, flux_x, params);
#if ADD_VIS
          // NOTE:给通量添加粘性，粘性项有主方向，需要交换速度
          DataArrayGPUf uy, vy, uy0, uy1, uy2, uy3, vy0, vy1, vy2, vy3;
          HydroState2d qNeighbors_r1, qNeighbors_r2; // i-1=qLocNeighbor,i=qLoc;

          uy0[0] = Qdata(i - 1, j - 2, IU);
          uy0[1] = Qdata(i - 1, j - 1, IU);
          uy0[2] = Qdata(i - 1, j + 1, IU);
          uy0[3] = Qdata(i - 1, j + 2, IU);

          uy1[0] = Qdata(i, j - 2, IU);
          uy1[1] = Qdata(i, j - 1, IU);
          uy1[2] = Qdata(i, j + 1, IU);
          uy1[3] = Qdata(i, j + 2, IU);

          uy2[0] = Qdata(i + 1, j - 2, IU);
          uy2[1] = Qdata(i + 1, j - 1, IU);
          uy2[2] = Qdata(i + 1, j + 1, IU);
          uy2[3] = Qdata(i + 1, j + 2, IU);

          uy3[0] = Qdata(i + 2, j - 2, IU);
          uy3[1] = Qdata(i + 2, j - 1, IU);
          uy3[2] = Qdata(i + 2, j + 1, IU);
          uy3[3] = Qdata(i + 2, j + 2, IU);

          uy[0] = 1.0 / 12.0 / params.dy * (uy0[0] - 8.0 * uy0[1] + 8.0 * uy0[2] - uy0[3]); // i-1
          uy[1] = 1.0 / 12.0 / params.dy * (uy1[0] - 8.0 * uy1[1] + 8.0 * uy1[2] - uy1[3]); // i
          uy[2] = 1.0 / 12.0 / params.dy * (uy2[0] - 8.0 * uy2[1] + 8.0 * uy2[2] - uy2[3]); // i+1
          uy[3] = 1.0 / 12.0 / params.dy * (uy3[0] - 8.0 * uy3[1] + 8.0 * uy3[2] - uy3[3]); // i+2

          vy0[0] = Qdata(i - 1, j - 2, IV);
          vy0[1] = Qdata(i - 1, j - 1, IV);
          vy0[2] = Qdata(i - 1, j + 1, IV);
          vy0[3] = Qdata(i - 1, j + 2, IV);

          vy1[0] = Qdata(i, j - 2, IV);
          vy1[1] = Qdata(i, j + 1, IV);
          vy1[2] = Qdata(i, j + 1, IV);
          vy1[3] = Qdata(i, j - 2, IV);

          vy2[0] = Qdata(i + 1, j - 2, IV);
          vy2[1] = Qdata(i + 1, j + 1, IV);
          vy2[2] = Qdata(i + 1, j + 1, IV);
          vy2[3] = Qdata(i + 1, j - 2, IV);

          vy3[0] = Qdata(i + 2, j - 2, IV);
          vy3[1] = Qdata(i + 2, j + 1, IV);
          vy3[2] = Qdata(i + 2, j + 1, IV);
          vy3[3] = Qdata(i + 2, j - 2, IV);

          vy[0] = 1.0 / 12.0 / params.dy * (vy0[0] - 8.0 * vy0[1] + 8.0 * vy0[2] - vy0[3]); // i-1
          vy[1] = 1.0 / 12.0 / params.dy * (vy1[0] - 8.0 * vy1[1] + 8.0 * vy1[2] - vy1[3]); // i
          vy[2] = 1.0 / 12.0 / params.dy * (vy2[0] - 8.0 * vy2[1] + 8.0 * vy2[2] - vy2[3]); // i+1
          vy[3] = 1.0 / 12.0 / params.dy * (vy3[0] - 8.0 * vy3[1] + 8.0 * vy3[2] - vy3[3]); // i+2

          qNeighbors_r1[ID] = Qdata(i + 1, j, ID);
          qNeighbors_r1[IP] = Qdata(i + 1, j, IP);
          qNeighbors_r1[IU] = Qdata(i + 1, j, IU);
          qNeighbors_r1[IV] = Qdata(i + 1, j, IV);
          qNeighbors_r1[IDY1] = Qdata(i + 1, j, IDY1);

          qNeighbors_r2[ID] = Qdata(i + 2, j, ID);
          qNeighbors_r2[IP] = Qdata(i + 2, j, IP);
          qNeighbors_r2[IU] = Qdata(i + 2, j, IU);
          qNeighbors_r2[IV] = Qdata(i + 2, j, IV);
          qNeighbors_r2[IDY1] = Qdata(i + 2, j, IDY1);

          if (params.if_debug)
          {

            printf("输出网格点位置(%d,%d),X方向上网格点的变量值,Yi=%lf，%lf,%lf,%lf \n", i, j, qLocNeighbor[IDY1], qLoc[IDY1], qNeighbors_r1[IDY1], qNeighbors_r2[IDY1]);
          }
          add_vflux(qLocNeighbor, qLoc, qNeighbors_r1, qNeighbors_r2, params, uy, vy, flux_x, params.dx, i, j);

#endif
#endif
          // store fluxes X
          /*NOTE:下面计算的就是用于计算U的半点通量，根据Godunov格式方程离散
          U_{i}^{n+1}=U_{i}^{n}-dt/dx*(Flux_{i+1/2}-Flux_{i-1/2})+dt*S_{i}
          所以可以推断下边的flux_x的组成
          */
          FluxData_x(i, j, ID) = flux_x[ID] * dtdx;
          FluxData_x(i, j, IP) = flux_x[IP] * dtdx;
          FluxData_x(i, j, IU) = flux_x[IU] * dtdx;
          FluxData_x(i, j, IV) = flux_x[IV] * dtdx;
#if COP
          FluxData_x(i, j, IDY1) = flux_x[IDY1] * dtdx;
#endif
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          // deal with left interface along Y !
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          // slopes at left neighbor along Y
          qLocNeighbor[ID] = Qdata(i, j - 1, ID);
          qNeighbors_0[ID] = Qdata(i + 1, j - 1, ID);
          qNeighbors_1[ID] = Qdata(i - 1, j - 1, ID);
          qNeighbors_2[ID] = Qdata(i, j, ID);
          qNeighbors_3[ID] = Qdata(i, j - 2, ID);

          qLocNeighbor[IP] = Qdata(i, j - 1, IP);
          qNeighbors_0[IP] = Qdata(i + 1, j - 1, IP);
          qNeighbors_1[IP] = Qdata(i - 1, j - 1, IP);
          qNeighbors_2[IP] = Qdata(i, j, IP);
          qNeighbors_3[IP] = Qdata(i, j - 2, IP);

          qLocNeighbor[IU] = Qdata(i, j - 1, IU);
          qNeighbors_0[IU] = Qdata(i + 1, j - 1, IU);
          qNeighbors_1[IU] = Qdata(i - 1, j - 1, IU);
          qNeighbors_2[IU] = Qdata(i, j, IU);
          qNeighbors_3[IU] = Qdata(i, j - 2, IU);

          qLocNeighbor[IV] = Qdata(i, j - 1, IV);
          qNeighbors_0[IV] = Qdata(i + 1, j - 1, IV);
          qNeighbors_1[IV] = Qdata(i - 1, j - 1, IV);
          qNeighbors_2[IV] = Qdata(i, j, IV);
          qNeighbors_3[IV] = Qdata(i, j - 2, IV);

#if COP
          qLocNeighbor[IDY1] = Qdata(i, j - 1, IDY1);
          qNeighbors_0[IDY1] = Qdata(i + 1, j - 1, IDY1);
          qNeighbors_1[IDY1] = Qdata(i - 1, j - 1, IDY1);
          qNeighbors_2[IDY1] = Qdata(i, j, IDY1);
          qNeighbors_3[IDY1] = Qdata(i, j - 2, IDY1);
#endif

          slope_unsplit_hydro_2d(qLocNeighbor,
                                 qNeighbors_0, qNeighbors_1,
                                 qNeighbors_2, qNeighbors_3,
                                 dqX_neighbor, dqY_neighbor);

          //
          // compute reconstructed states at left interface along Y
          //

          // left interface : right state
          trace_unsplit_2d_along_dir(qLoc,
                                     dqX, dqY,
                                     dtdx, dtdy, FACE_YMIN, qright);

          // left interface : left state
          trace_unsplit_2d_along_dir(qLocNeighbor,
                                     dqX_neighbor, dqY_neighbor,
                                     dtdx, dtdy, FACE_YMAX, qleft);

          if (gravity_enabled)
          {
            // we need to modify input to flux computation with
            // gravity predictor (half time step)

            qleft[IU] += 0.5 * dt * gravity(i, j - 1, IX);
            qleft[IV] += 0.5 * dt * gravity(i, j - 1, IY);

            qright[IU] += 0.5 * dt * gravity(i, j, IX);
            qright[IV] += 0.5 * dt * gravity(i, j, IY);
          }

          // Solve Riemann problem at Y-interfaces and compute Y-fluxes
          swapValues(&(qleft[IU]), &(qleft[IV]));
          swapValues(&(qright[IU]), &(qright[IV]));
          // riemann_2d(qleft,qright,qgdnv,flux_y);
#if COP
          riemann_hydro(qleft, qright, qgdnv, flux_y, params);
#if ADD_VIS
          uy0[0] = Qdata(i - 2, j - 1, IV);
          uy0[1] = Qdata(i - 1, j - 1, IV);
          uy0[2] = Qdata(i + 1, j - 1, IV);
          uy0[3] = Qdata(i + 2, j - 1, IV);

          uy1[0] = Qdata(i - 2, j, IV);
          uy1[1] = Qdata(i - 1, j, IV);
          uy1[2] = Qdata(i + 1, j, IV);
          uy1[3] = Qdata(i + 2, j, IV);

          uy2[0] = Qdata(i - 2, j + 1, IV);
          uy2[1] = Qdata(i - 1, j + 1, IV);
          uy2[2] = Qdata(i + 1, j + 1, IV);
          uy2[3] = Qdata(i + 2, j + 1, IV);

          uy3[0] = Qdata(i - 2, j + 2, IV);
          uy3[1] = Qdata(i - 1, j + 2, IV);
          uy3[2] = Qdata(i + 1, j + 2, IV);
          uy3[3] = Qdata(i + 2, j + 2, IV);

          uy[0] = 1.0 / 12.0 / params.dx * (uy0[0] - 8.0 * uy0[1] + 8.0 * uy0[2] - uy0[3]); // i-1//在y方向上对应的物理量其实是vx
          uy[1] = 1.0 / 12.0 / params.dx * (uy1[0] - 8.0 * uy1[1] + 8.0 * uy1[2] - uy1[3]); // i
          uy[2] = 1.0 / 12.0 / params.dx * (uy2[0] - 8.0 * uy2[1] + 8.0 * uy2[2] - uy2[3]); // i+1
          uy[3] = 1.0 / 12.0 / params.dx * (uy3[0] - 8.0 * uy3[1] + 8.0 * uy3[2] - uy3[3]); // i+2

          vy0[0] = Qdata(i - 2, j - 1, IU);
          vy0[1] = Qdata(i - 1, j - 1, IU);
          vy0[2] = Qdata(i + 1, j - 1, IU);
          vy0[3] = Qdata(i + 2, j - 1, IU);

          vy1[0] = Qdata(i - 2, j, IU);
          vy1[1] = Qdata(i - 1, j, IU);
          vy1[2] = Qdata(i + 1, j, IU);
          vy1[3] = Qdata(i + 2, j, IU);

          vy2[0] = Qdata(i - 2, j + 1, IU);
          vy2[1] = Qdata(i - 1, j + 1, IU);
          vy2[2] = Qdata(i + 1, j + 1, IU);
          vy2[3] = Qdata(i + 2, j + 1, IU);

          vy3[0] = Qdata(i - 2, j + 2, IU);
          vy3[1] = Qdata(i - 1, j + 2, IU);
          vy3[2] = Qdata(i + 1, j + 2, IU);
          vy3[3] = Qdata(i + 2, j + 2, IU);

          vy[0] = 1.0 / 12.0 / params.dx * (vy0[0] - 8.0 * vy0[1] + 8.0 * vy0[2] - vy0[3]); // i-1
          vy[1] = 1.0 / 12.0 / params.dx * (vy1[0] - 8.0 * vy1[1] + 8.0 * vy1[2] - vy1[3]); // i
          vy[2] = 1.0 / 12.0 / params.dx * (vy2[0] - 8.0 * vy2[1] + 8.0 * vy2[2] - vy2[3]); // i+1
          vy[3] = 1.0 / 12.0 / params.dx * (vy3[0] - 8.0 * vy3[1] + 8.0 * vy3[2] - vy3[3]); // i+2

          qNeighbors_r1[ID] = Qdata(i, j + 1, ID);
          qNeighbors_r1[IP] = Qdata(i, j + 1, IP);
          qNeighbors_r1[IU] = Qdata(i, j + 1, IU);
          qNeighbors_r1[IV] = Qdata(i, j + 1, IV);
          qNeighbors_r1[IDY1] = Qdata(i, j + 1, IDY1);

          qNeighbors_r2[ID] = Qdata(i, j + 2, ID);
          qNeighbors_r2[IP] = Qdata(i, j + 2, IP);
          qNeighbors_r2[IU] = Qdata(i, j + 2, IU);
          qNeighbors_r2[IV] = Qdata(i, j + 2, IV);
          qNeighbors_r2[IDY1] = Qdata(i, j + 2, IDY1);

          swapValues(&(qLocNeighbor[IU]), &(qLocNeighbor[IV]));
          swapValues(&(qLoc[IU]), &(qLoc[IV]));
          swapValues(&(qNeighbors_r1[IU]), &(qNeighbors_r1[IV]));
          swapValues(&(qNeighbors_r2[IU]), &(qNeighbors_r2[IV]));
          swapValues(&(flux_y[IU]), &(flux_y[IV]));
          if (params.if_debug)
          {
            printf("输出网格点位置(%d,%d),Y方向上网格点的变量值,Yi=%lf，%lf,%lf,%lf \n", i, j, qLocNeighbor[IDY1], qLoc[IDY1], qNeighbors_r1[IDY1], qNeighbors_r2[IDY1]);
          }
          add_vflux(qLocNeighbor, qLoc, qNeighbors_r1, qNeighbors_r2, params, uy, vy, flux_y, params.dy, i, j);
          swapValues(&(flux_y[IU]), &(flux_y[IV]));
#endif
#endif
          //
          // store fluxes Y
          //
          FluxData_y(i, j, ID) = flux_y[ID] * dtdy;
          FluxData_y(i, j, IP) = flux_y[IP] * dtdy;
          FluxData_y(i, j, IU) = flux_y[IV] * dtdy; //
          FluxData_y(i, j, IV) = flux_y[IU] * dtdy; //
#if COP
          FluxData_y(i, j, IDY1) = flux_y[IDY1] * dtdy;
#endif

        } // end if

      } // end operator ()

      DataArray2d Qdata;
      DataArray2d FluxData_x;
      DataArray2d FluxData_y;
      real_t dt, dtdx, dtdy;
      bool gravity_enabled;
      VectorField2d gravity;
    }; // ComputeAndStoreFluxesFunctor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    class UpdateFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Perform time update using the stored fluxes.
       *
       * \note this functor must be called after ComputeAndStoreFluxesFunctor2D
       *
       * \param[in,out] Udata
       * \param[in] FluxData_x flux coming from the left neighbor along X
       * \param[in] FluxData_y flux coming from the left neighbor along Y
       */
      UpdateFunctor2D(HydroParams params,
                      DataArray2d Udata,
                      DataArray2d FluxData_x,
                      DataArray2d FluxData_y) : HydroBaseFunctor2D(params),
                                                Udata(Udata),
                                                FluxData_x(FluxData_x),
                                                FluxData_y(FluxData_y){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        DataArray2d Udata,
                        DataArray2d FluxData_x,
                        DataArray2d FluxData_y)
      {
        int nbCells = params.isize * params.jsize;
        UpdateFunctor2D functor(params, Udata, FluxData_x, FluxData_y);
        Kokkos::parallel_for(nbCells, functor);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth && j < jsize - ghostWidth &&
            i >= ghostWidth && i < isize - ghostWidth)
        {

          Udata(i, j, ID) += FluxData_x(i, j, ID);
          Udata(i, j, IP) += FluxData_x(i, j, IP);
          Udata(i, j, IU) += FluxData_x(i, j, IU);
          Udata(i, j, IV) += FluxData_x(i, j, IV);

          Udata(i, j, ID) -= FluxData_x(i + 1, j, ID);
          Udata(i, j, IP) -= FluxData_x(i + 1, j, IP);
          Udata(i, j, IU) -= FluxData_x(i + 1, j, IU);
          Udata(i, j, IV) -= FluxData_x(i + 1, j, IV);

          Udata(i, j, ID) += FluxData_y(i, j, ID);
          Udata(i, j, IP) += FluxData_y(i, j, IP);
          Udata(i, j, IU) += FluxData_y(i, j, IU);
          Udata(i, j, IV) += FluxData_y(i, j, IV);

          Udata(i, j, ID) -= FluxData_y(i, j + 1, ID);
          Udata(i, j, IP) -= FluxData_y(i, j + 1, IP);
          Udata(i, j, IU) -= FluxData_y(i, j + 1, IU);
          Udata(i, j, IV) -= FluxData_y(i, j + 1, IV);

#if COP
          Udata(i, j, IDY1) += FluxData_x(i, j, IDY1);
          Udata(i, j, IDY1) -= FluxData_x(i + 1, j, IDY1);
          Udata(i, j, IDY1) += FluxData_y(i, j, IDY1);
          Udata(i, j, IDY1) -= FluxData_y(i, j + 1, IDY1);
#endif
        } // end if

      } // end operator ()

      DataArray2d Udata;
      DataArray2d FluxData_x;
      DataArray2d FluxData_y;

    }; // UpdateFunctor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    template <Direction dir>
    class UpdateDirFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Perform time update using the stored fluxes along direction dir.
       *
       * \param[in,out] Udata
       * \param[in] FluxData flux coming from the left neighbor along direction dir
       *
       */
      UpdateDirFunctor2D(HydroParams params,
                         DataArray2d Udata,
                         DataArray2d FluxData) : HydroBaseFunctor2D(params),
                                                 Udata(Udata),
                                                 FluxData(FluxData){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        DataArray2d Udata,
                        DataArray2d FluxData)
      {
        int nbCells = params.isize * params.jsize;
        UpdateDirFunctor2D<dir> functor(params, Udata, FluxData);
        Kokkos::parallel_for(nbCells, functor);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth && j < jsize - ghostWidth &&
            i >= ghostWidth && i < isize - ghostWidth)
        {

          if (dir == XDIR)
          {

            Udata(i, j, ID) += FluxData(i, j, ID);
            Udata(i, j, IP) += FluxData(i, j, IP);
            Udata(i, j, IU) += FluxData(i, j, IU);
            Udata(i, j, IV) += FluxData(i, j, IV);

            Udata(i, j, ID) -= FluxData(i + 1, j, ID);
            Udata(i, j, IP) -= FluxData(i + 1, j, IP);
            Udata(i, j, IU) -= FluxData(i + 1, j, IU);
            Udata(i, j, IV) -= FluxData(i + 1, j, IV);
          }
          else if (dir == YDIR)
          {

            Udata(i, j, ID) += FluxData(i, j, ID);
            Udata(i, j, IP) += FluxData(i, j, IP);
            Udata(i, j, IU) += FluxData(i, j, IU);
            Udata(i, j, IV) += FluxData(i, j, IV);

            Udata(i, j, ID) -= FluxData(i, j + 1, ID);
            Udata(i, j, IP) -= FluxData(i, j + 1, IP);
            Udata(i, j, IU) -= FluxData(i, j + 1, IU);
            Udata(i, j, IV) -= FluxData(i, j + 1, IV);
          }

        } // end if

      } // end operator ()

      DataArray2d Udata;
      DataArray2d FluxData;

    }; // UpdateDirFunctor

    /*************************************************/
    /*************************************************/
    /*************************************************/
    class ComputeSlopesFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Compute limited slopes.
       *
       * \param[in] Qdata primitive variables
       * \param[out] Slopes_x limited slopes along direction X
       * \param[out] Slopes_y limited slopes along direction Y
       */
      ComputeSlopesFunctor2D(HydroParams params,
                             DataArray2d Qdata,
                             DataArray2d Slopes_x,
                             DataArray2d Slopes_y) : HydroBaseFunctor2D(params), Qdata(Qdata),
                                                     Slopes_x(Slopes_x), Slopes_y(Slopes_y){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        DataArray2d Qdata,
                        DataArray2d Slopes_x,
                        DataArray2d Slopes_y)
      {
        int nbCells = params.isize * params.jsize;
        ComputeSlopesFunctor2D functor(params, Qdata, Slopes_x, Slopes_y);
        Kokkos::parallel_for(nbCells, functor);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth - 1 && j <= jsize - ghostWidth &&
            i >= ghostWidth - 1 && i <= isize - ghostWidth)
        {

          // local primitive variables
          HydroState qLoc; // local primitive variables

          // local primitive variables in neighborbood
          HydroState qNeighbors_0;
          HydroState qNeighbors_1;
          HydroState qNeighbors_2;
          HydroState qNeighbors_3;

          // Local slopes and neighbor slopes
          HydroState dqX{};
          HydroState dqY{};

          // get primitive variables state vector
          qLoc[ID] = Qdata(i, j, ID);
          qNeighbors_0[ID] = Qdata(i + 1, j, ID);
          qNeighbors_1[ID] = Qdata(i - 1, j, ID);
          qNeighbors_2[ID] = Qdata(i, j + 1, ID);
          qNeighbors_3[ID] = Qdata(i, j - 1, ID);

          qLoc[IP] = Qdata(i, j, IP);
          qNeighbors_0[IP] = Qdata(i + 1, j, IP);
          qNeighbors_1[IP] = Qdata(i - 1, j, IP);
          qNeighbors_2[IP] = Qdata(i, j + 1, IP);
          qNeighbors_3[IP] = Qdata(i, j - 1, IP);

          qLoc[IU] = Qdata(i, j, IU);
          qNeighbors_0[IU] = Qdata(i + 1, j, IU);
          qNeighbors_1[IU] = Qdata(i - 1, j, IU);
          qNeighbors_2[IU] = Qdata(i, j + 1, IU);
          qNeighbors_3[IU] = Qdata(i, j - 1, IU);

          qLoc[IV] = Qdata(i, j, IV);
          qNeighbors_0[IV] = Qdata(i + 1, j, IV);
          qNeighbors_1[IV] = Qdata(i - 1, j, IV);
          qNeighbors_2[IV] = Qdata(i, j + 1, IV);
          qNeighbors_3[IV] = Qdata(i, j - 1, IV);

          slope_unsplit_hydro_2d(qLoc,
                                 qNeighbors_0, qNeighbors_1,
                                 qNeighbors_2, qNeighbors_3,
                                 dqX, dqY);

          // copy back slopes in global arrays
          Slopes_x(i, j, ID) = dqX[ID];
          Slopes_y(i, j, ID) = dqY[ID];

          Slopes_x(i, j, IP) = dqX[IP];
          Slopes_y(i, j, IP) = dqY[IP];

          Slopes_x(i, j, IU) = dqX[IU];
          Slopes_y(i, j, IU) = dqY[IU];

          Slopes_x(i, j, IV) = dqX[IV];
          Slopes_y(i, j, IV) = dqY[IV];

        } // end if

      } // end operator ()

      DataArray2d Qdata;
      DataArray2d Slopes_x, Slopes_y;

    }; // ComputeSlopesFunctor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    template <Direction dir>
    class ComputeTraceAndFluxes_Functor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Compute reconstructed states on faces (not stored), and fluxes (stored).
       *
       * \param[in] Qdata primitive variables
       * \param[in] Slopes_x limited slopes along direction X
       * \param[in] Slopes_y limited slopes along direction Y
       * \param[out] Fluxes along direction dir
       *
       * \tparam dir direction along which fluxes are computed.
       */
      ComputeTraceAndFluxes_Functor2D(HydroParams params,
                                      DataArray2d Qdata,
                                      DataArray2d Slopes_x,
                                      DataArray2d Slopes_y,
                                      DataArray2d Fluxes,
                                      real_t dt,
                                      bool gravity_enabled,
                                      VectorField2d gravity) : HydroBaseFunctor2D(params), Qdata(Qdata),
                                                               Slopes_x(Slopes_x), Slopes_y(Slopes_y),
                                                               Fluxes(Fluxes),
                                                               dt(dt),
                                                               dtdx(dt / params.dx), dtdy(dt / params.dy),
                                                               gravity_enabled(gravity_enabled),
                                                               gravity(gravity){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        DataArray2d Qdata,
                        DataArray2d Slopes_x,
                        DataArray2d Slopes_y,
                        DataArray2d Fluxes,
                        real_t dt,
                        bool gravity_enabled,
                        VectorField2d gravity)
      {
        int nbCells = params.isize * params.jsize;
        ComputeTraceAndFluxes_Functor2D<dir> functor(params, Qdata,
                                                     Slopes_x, Slopes_y,
                                                     Fluxes,
                                                     dt,
                                                     gravity_enabled,
                                                     gravity);
        Kokkos::parallel_for(nbCells, functor);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth && j <= jsize - ghostWidth &&
            i >= ghostWidth && i <= isize - ghostWidth)
        {

          // local primitive variables
          HydroState qLoc; // local primitive variables

          // local primitive variables in neighbor cell
          HydroState qLocNeighbor;

          // Local slopes and neighbor slopes
          HydroState dqX;
          HydroState dqY;
          HydroState dqX_neighbor;
          HydroState dqY_neighbor;

          // Local variables for Riemann problems solving
          HydroState qleft;
          HydroState qright;
          HydroState qgdnv;
          HydroState flux;

          //
          // compute reconstructed states at left interface along X
          //
          qLoc[ID] = Qdata(i, j, ID);
          dqX[ID] = Slopes_x(i, j, ID);
          dqY[ID] = Slopes_y(i, j, ID);

          qLoc[IP] = Qdata(i, j, IP);
          dqX[IP] = Slopes_x(i, j, IP);
          dqY[IP] = Slopes_y(i, j, IP);

          qLoc[IU] = Qdata(i, j, IU);
          dqX[IU] = Slopes_x(i, j, IU);
          dqY[IU] = Slopes_y(i, j, IU);

          qLoc[IV] = Qdata(i, j, IV);
          dqX[IV] = Slopes_x(i, j, IV);
          dqY[IV] = Slopes_y(i, j, IV);

          if (dir == XDIR)
          {

            // left interface : right state
            trace_unsplit_2d_along_dir(qLoc,
                                       dqX, dqY,
                                       dtdx, dtdy, FACE_XMIN, qright);

            if (gravity_enabled)
            {
              // we need to modify input to flux computation with
              // gravity predictor (half time step)

              qright[IU] += 0.5 * dt * gravity(i, j, IX);
              qright[IV] += 0.5 * dt * gravity(i, j, IY);
            }

            qLocNeighbor[ID] = Qdata(i - 1, j, ID);
            dqX_neighbor[ID] = Slopes_x(i - 1, j, ID);
            dqY_neighbor[ID] = Slopes_y(i - 1, j, ID);

            qLocNeighbor[IP] = Qdata(i - 1, j, IP);
            dqX_neighbor[IP] = Slopes_x(i - 1, j, IP);
            dqY_neighbor[IP] = Slopes_y(i - 1, j, IP);

            qLocNeighbor[IU] = Qdata(i - 1, j, IU);
            dqX_neighbor[IU] = Slopes_x(i - 1, j, IU);
            dqY_neighbor[IU] = Slopes_y(i - 1, j, IU);

            qLocNeighbor[IV] = Qdata(i - 1, j, IV);
            dqX_neighbor[IV] = Slopes_x(i - 1, j, IV);
            dqY_neighbor[IV] = Slopes_y(i - 1, j, IV);

            // left interface : left state
            trace_unsplit_2d_along_dir(qLocNeighbor,
                                       dqX_neighbor, dqY_neighbor,
                                       dtdx, dtdy, FACE_XMAX, qleft);

            if (gravity_enabled)
            {
              // we need to modify input to flux computation with
              // gravity predictor (half time step)

              qleft[IU] += 0.5 * dt * gravity(i - 1, j, IX);
              qleft[IV] += 0.5 * dt * gravity(i - 1, j, IY);
            }

            // Solve Riemann problem at X-interfaces and compute X-fluxes
            riemann_hydro(qleft, qright, qgdnv, flux, params);

            //
            // store fluxes
            //
            Fluxes(i, j, ID) = flux[ID] * dtdx;
            Fluxes(i, j, IP) = flux[IP] * dtdx;
            Fluxes(i, j, IU) = flux[IU] * dtdx;
            Fluxes(i, j, IV) = flux[IV] * dtdx;
          }
          else if (dir == YDIR)
          {

            // left interface : right state
            trace_unsplit_2d_along_dir(qLoc,
                                       dqX, dqY,
                                       dtdx, dtdy, FACE_YMIN, qright);

            if (gravity_enabled)
            {
              // we need to modify input to flux computation with
              // gravity predictor (half time step)

              qright[IU] += 0.5 * dt * gravity(i, j, IX);
              qright[IV] += 0.5 * dt * gravity(i, j, IY);
            }

            qLocNeighbor[ID] = Qdata(i, j - 1, ID);
            dqX_neighbor[ID] = Slopes_x(i, j - 1, ID);
            dqY_neighbor[ID] = Slopes_y(i, j - 1, ID);

            qLocNeighbor[IP] = Qdata(i, j - 1, IP);
            dqX_neighbor[IP] = Slopes_x(i, j - 1, IP);
            dqY_neighbor[IP] = Slopes_y(i, j - 1, IP);

            qLocNeighbor[IU] = Qdata(i, j - 1, IU);
            dqX_neighbor[IU] = Slopes_x(i, j - 1, IU);
            dqY_neighbor[IU] = Slopes_y(i, j - 1, IU);

            qLocNeighbor[IV] = Qdata(i, j - 1, IV);
            dqX_neighbor[IV] = Slopes_x(i, j - 1, IV);
            dqY_neighbor[IV] = Slopes_y(i, j - 1, IV);

            // left interface : left state
            trace_unsplit_2d_along_dir(qLocNeighbor,
                                       dqX_neighbor, dqY_neighbor,
                                       dtdx, dtdy, FACE_YMAX, qleft);

            if (gravity_enabled)
            {
              // we need to modify input to flux computation with
              // gravity predictor (half time step)

              qleft[IU] += 0.5 * dt * gravity(i, j - 1, IX);
              qleft[IV] += 0.5 * dt * gravity(i, j - 1, IY);
            }

            // Solve Riemann problem at Y-interfaces and compute Y-fluxes
            swapValues(&(qleft[IU]), &(qleft[IV]));
            swapValues(&(qright[IU]), &(qright[IV]));
            riemann_hydro(qleft, qright, qgdnv, flux, params);

            //
            // update hydro array
            //
            Fluxes(i, j, ID) = flux[ID] * dtdy;
            Fluxes(i, j, IP) = flux[IP] * dtdy;
            Fluxes(i, j, IU) = flux[IV] * dtdy; // IU/IV swapped
            Fluxes(i, j, IV) = flux[IU] * dtdy; // IU/IV swapped
          }

        } // end if

      } // end operator ()

      DataArray2d Qdata;
      DataArray2d Slopes_x, Slopes_y;
      DataArray2d Fluxes;
      real_t dt, dtdx, dtdy;
      bool gravity_enabled;
      VectorField2d gravity;

    }; // ComputeTraceAndFluxes_Functor2D

    /*************************************************/
    /*************************************************/
    /*************************************************/
    class GravitySourceTermFunctor2D : public HydroBaseFunctor2D
    {

    public:
      /**
       * Update with gravity source term.
       *
       * \param[in] Udata_in conservative variables at t(n)
       * \param[in,out] Udata_out conservative variables at t(n+1)
       * \param[in] gravity is a vector field
       */
      GravitySourceTermFunctor2D(HydroParams params,
                                 DataArray2d Udata_in,
                                 DataArray2d Udata_out,
                                 VectorField2d gravity,
                                 real_t dt) : HydroBaseFunctor2D(params),
                                              Udata_in(Udata_in),
                                              Udata_out(Udata_out),
                                              gravity(gravity),
                                              dt(dt){};

      // static method which does it all: create and execute functor
      static void apply(HydroParams params,
                        DataArray2d Udata_in,
                        DataArray2d Udata_out,
                        VectorField2d gravity,
                        real_t dt)
      {
        int nbCells = params.isize * params.jsize;
        GravitySourceTermFunctor2D functor(params, Udata_in, Udata_out, gravity, dt);
        Kokkos::parallel_for(nbCells, functor);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const int &index) const
      {
        const int isize = params.isize;
        const int jsize = params.jsize;
        const int ghostWidth = params.ghostWidth;

        int i, j;
        index2coord(index, i, j, isize, jsize);

        if (j >= ghostWidth && j < jsize - ghostWidth &&
            i >= ghostWidth && i < isize - ghostWidth)
        {

          real_t rhoOld = Udata_in(i, j, ID);
          real_t rhoNew = fmax(params.settings.smallr, Udata_out(i, j, ID));

          real_t rhou = Udata_out(i, j, IU);
          real_t rhov = Udata_out(i, j, IV);

          // compute kinetic energy before updating momentum
          real_t ekin_old = 0.5 * (rhou * rhou + rhov * rhov) / rhoNew;

          // update momentum
          rhou += 0.5 * dt * gravity(i, j, IX) * (rhoOld + rhoNew);
          rhov += 0.5 * dt * gravity(i, j, IY) * (rhoOld + rhoNew);
          Udata_out(i, j, IU) = rhou;
          Udata_out(i, j, IV) = rhov;

          // compute kinetic energy after updating momentum
          real_t ekin_new = 0.5 * (rhou * rhou + rhov * rhov) / rhoNew;

          // update total energy
          Udata_out(i, j, IE) += (ekin_new - ekin_old);
        }

      } // end operator ()

      DataArray2d Udata_in, Udata_out;
      VectorField2d gravity;
      real_t dt;

    }; // GravitySourceTermFunctor2D

  } // namespace muscl

} // namespace euler_kokkos

#endif // HYDRO_RUN_FUNCTORS_2D_H_
