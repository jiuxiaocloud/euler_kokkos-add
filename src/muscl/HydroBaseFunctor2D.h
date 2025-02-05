#ifndef HYDRO_BASE_FUNCTOR_2D_H_
#define HYDRO_BASE_FUNCTOR_2D_H_

#include "shared/kokkos_shared.h"

#include "shared/HydroParams.h"
#include "shared/HydroState.h"

namespace euler_kokkos
{
  namespace muscl
  {

    /**
     * Base class to derive actual kokkos functor for hydro 2D.
     * params is passed by copy.
     */
    class HydroBaseFunctor2D
    {

    public:
      using HydroState = HydroState2d;
      using DataArray = DataArray2d;

      HydroBaseFunctor2D(HydroParams params) : params(params){};
      virtual ~HydroBaseFunctor2D(){};

      HydroParams params;
#if COP
      // Add_Params add_params;
      // HydroBaseFunctor2D(HydroParams params, Add_Params add_params) : params(params), add_params(add_params){};
#endif
      const int nbvar = params.nbvar;

      // utility routines used in various computational kernels

      KOKKOS_INLINE_FUNCTION
      void swapValues(real_t *a, real_t *b) const
      {

        real_t tmp = *a;

        *a = *b;
        *b = tmp;

      } // swapValues

      /**
       * Equation of state:
       * compute pressure p and speed of sound c, from density rho and
       * internal energy eint using the "calorically perfect gas" equation
       * of state : \f$ eint=\frac{p}{\rho (\gamma-1)} \f$
       * Recall that \f$ \gamma \f$ is equal to the ratio of specific heats
       *  \f$ \left[ c_p/c_v \right] \f$.
       *
       * @param[in]  rho  density
       * @param[in]  eint internal energy
       * @param[out] p    pressure
       * @param[out] c    speed of sound
       */
      KOKKOS_INLINE_FUNCTION
      void eos(real_t rho,
               real_t eint,
               real_t *p,
               real_t *c) const
      {
        real_t gamma0 = params.settings.gamma0;
        real_t smallp = params.settings.smallp;

        *p = FMAX((gamma0 - ONE_F) * rho * eint, rho * smallp);
        *c = SQRT(gamma0 * (*p) / rho);

      } // eos

      /**
       * Convert conservative variables (rho, rho*u, rho*v, e) to
       * primitive variables (rho,u,v,p)
       * @param[in]  u  conservative variables array
       * @param[out] q  primitive    variables array (allocated in calling routine, size is constant nbvar)
       * @param[out] c  local speed of sound
       */
      KOKKOS_INLINE_FUNCTION
      void computePrimitives(const HydroState &u,
                             real_t *c,
                             HydroState &q) const
      {
        real_t gamma0 = params.settings.gamma0;
        real_t smallr = params.settings.smallr;
        real_t smallp = params.settings.smallp;

        real_t d, p, ux, uy;

        d = fmax(u[ID], smallr);
        ux = u[IU] / d;
        uy = u[IV] / d;

        real_t eken = HALF_F * (ux * ux + uy * uy);
        real_t e = u[IP] / d - eken;

        // compute pressure and speed of sound
        p = fmax((gamma0 - 1.0) * d * e, d * smallp);
        *c = sqrt(gamma0 * (p) / d);

        q[ID] = d;
        q[IP] = p;
        q[IU] = ux;
        q[IV] = uy;
#if COP
        real_t dy1 = (u[IDY1] / d);
        q[IDY1] = fmin(fmax(dy1, 1.0e-10), 1.0);
#endif
      } // NOTE:computePrimitive

      /**
       * Trace computations for unsplit Godunov scheme.
       *
       * \param[in] q          : Primitive variables state.
       * \param[in] qNeighbors : state in the neighbor cells (2 neighbors
       * per dimension, in the following order x+, x-, y+, y-, z+, z-)
       * \param[in] c          : local sound speed.
       * \param[in] dtdx       : dt over dx
       * \param[out] qm        : qm state (one per dimension)
       * \param[out] qp        : qp state (one per dimension)
       */
      KOKKOS_INLINE_FUNCTION
      void trace_unsplit_2d(const HydroState &q,
                            const HydroState &qNeighbors_0,
                            const HydroState &qNeighbors_1,
                            const HydroState &qNeighbors_2,
                            const HydroState &qNeighbors_3,
                            real_t c,
                            real_t dtdx,
                            real_t dtdy,
                            HydroState &qm_x,
                            HydroState &qm_y,
                            HydroState &qp_x,
                            HydroState &qp_y) const
      {

        real_t gamma0 = params.settings.gamma0;
        real_t smallr = params.settings.smallr;

        // first compute slopes
        HydroState dqX, dqY;
        dqX[ID] = 0.0;
        dqX[IP] = 0.0;
        dqX[IU] = 0.0;
        dqX[IV] = 0.0;
        dqY[ID] = 0.0;
        dqY[IP] = 0.0;
        dqY[IU] = 0.0;
        dqY[IV] = 0.0;

        slope_unsplit_hydro_2d(q,
                               qNeighbors_0, qNeighbors_1,
                               qNeighbors_2, qNeighbors_3,
                               dqX, dqY);

        // Cell centered values
        real_t r = q[ID];
        real_t p = q[IP];
        real_t u = q[IU];
        real_t v = q[IV];

        // TVD slopes in all directions
        real_t drx = dqX[ID];
        real_t dpx = dqX[IP];
        real_t dux = dqX[IU];
        real_t dvx = dqX[IV];

        real_t dry = dqY[ID];
        real_t dpy = dqY[IP];
        real_t duy = dqY[IU];
        real_t dvy = dqY[IV];

        // source terms (with transverse derivatives)
        real_t sr0 = -u * drx - v * dry - (dux + dvy) * r;
        real_t sp0 = -u * dpx - v * dpy - (dux + dvy) * gamma0 * p;
        real_t su0 = -u * dux - v * duy - (dpx) / r;
        real_t sv0 = -u * dvx - v * dvy - (dpy) / r;

        // Right state at left interface
        qp_x[ID] = r - HALF_F * drx + sr0 * dtdx * HALF_F;
        qp_x[IP] = p - HALF_F * dpx + sp0 * dtdx * HALF_F;
        qp_x[IU] = u - HALF_F * dux + su0 * dtdx * HALF_F;
        qp_x[IV] = v - HALF_F * dvx + sv0 * dtdx * HALF_F;
        qp_x[ID] = fmax(smallr, qp_x[ID]);

        // Left state at right interface
        qm_x[ID] = r + HALF_F * drx + sr0 * dtdx * HALF_F;
        qm_x[IP] = p + HALF_F * dpx + sp0 * dtdx * HALF_F;
        qm_x[IU] = u + HALF_F * dux + su0 * dtdx * HALF_F;
        qm_x[IV] = v + HALF_F * dvx + sv0 * dtdx * HALF_F;
        qm_x[ID] = fmax(smallr, qm_x[ID]);

        // Top state at bottom interface
        qp_y[ID] = r - HALF_F * dry + sr0 * dtdy * HALF_F;
        qp_y[IP] = p - HALF_F * dpy + sp0 * dtdy * HALF_F;
        qp_y[IU] = u - HALF_F * duy + su0 * dtdy * HALF_F;
        qp_y[IV] = v - HALF_F * dvy + sv0 * dtdy * HALF_F;
        qp_y[ID] = fmax(smallr, qp_y[ID]);

        // Bottom state at top interface
        qm_y[ID] = r + HALF_F * dry + sr0 * dtdy * HALF_F;
        qm_y[IP] = p + HALF_F * dpy + sp0 * dtdy * HALF_F;
        qm_y[IU] = u + HALF_F * duy + su0 * dtdy * HALF_F;
        qm_y[IV] = v + HALF_F * dvy + sv0 * dtdy * HALF_F;
        qm_y[ID] = fmax(smallr, qm_y[ID]);

      } // trace_unsplit_2d

      /**
       * Trace computations for unsplit Godunov scheme.
       *
       * \param[in] q          : Primitive variables state.
       * \param[in] dqX        : slope along X
       * \param[in] dqY        : slope along Y
       * \param[in] dtdx       : dt over dx
       * \param[in] dtdy       : dt over dy
       * \param[in] faceId     : which face will be reconstructed
       * \param[out] qface     : q reconstructed state at cell interface
       */
      KOKKOS_INLINE_FUNCTION
      void trace_unsplit_2d_along_dir(const HydroState &q,
                                      const HydroState &dqX,
                                      const HydroState &dqY,
                                      real_t dtdx,
                                      real_t dtdy,
                                      int faceId,
                                      HydroState &qface) const
      {

        real_t gamma0 = params.settings.gamma0;
        real_t smallr = params.settings.smallr;

        // Cell centered values
        real_t r = q[ID];
        real_t p = q[IP];
        real_t u = q[IU];
        real_t v = q[IV];

        // TVD slopes in all directions
        real_t drx = dqX[ID];
        real_t dpx = dqX[IP];
        real_t dux = dqX[IU];
        real_t dvx = dqX[IV];

        real_t dry = dqY[ID];
        real_t dpy = dqY[IP];
        real_t duy = dqY[IU];
        real_t dvy = dqY[IV];
#if COP
        real_t y1 = q[IDY1];
        real_t dy1x = dqX[IDY1];
        real_t dy1y = dqY[IDY1];
#endif

        // source terms (with transverse derivatives)
        real_t sr0 = -u * drx - v * dry - (dux + dvy) * r;
        real_t sp0 = -u * dpx - v * dpy - (dux + dvy) * gamma0 * p;
        real_t su0 = -u * dux - v * duy - (dpx) / r;
        real_t sv0 = -u * dvx - v * dvy - (dpy) / r;
#if COP
        real_t sy10 = -u * dy1x - v * dy1y - (dux + dvy) * r;
#endif

        if (faceId == FACE_XMIN)
        {
          // Right state at left interface
          qface[ID] = r - HALF_F * drx + sr0 * dtdx * HALF_F;
          qface[IP] = p - HALF_F * dpx + sp0 * dtdx * HALF_F;
          qface[IU] = u - HALF_F * dux + su0 * dtdx * HALF_F;
          qface[IV] = v - HALF_F * dvx + sv0 * dtdx * HALF_F;
          qface[ID] = fmax(smallr, qface[ID]);
#if COP
          qface[IDY1] = y1 - HALF_F * dy1x + sy10 * dtdx * HALF_F;
#endif
        }

        if (faceId == FACE_XMAX)
        {
          // Left state at right interface
          qface[ID] = r + HALF_F * drx + sr0 * dtdx * HALF_F;
          qface[IP] = p + HALF_F * dpx + sp0 * dtdx * HALF_F;
          qface[IU] = u + HALF_F * dux + su0 * dtdx * HALF_F;
          qface[IV] = v + HALF_F * dvx + sv0 * dtdx * HALF_F;
          qface[ID] = fmax(smallr, qface[ID]);
#if COP
          qface[IDY1] = y1 + HALF_F * dy1x + sy10 * dtdx * HALF_F;
#endif
        }

        if (faceId == FACE_YMIN)
        {
          // Top state at bottom interface
          qface[ID] = r - HALF_F * dry + sr0 * dtdy * HALF_F;
          qface[IP] = p - HALF_F * dpy + sp0 * dtdy * HALF_F;
          qface[IU] = u - HALF_F * duy + su0 * dtdy * HALF_F;
          qface[IV] = v - HALF_F * dvy + sv0 * dtdy * HALF_F;
          qface[ID] = fmax(smallr, qface[ID]);
#if COP
          qface[IDY1] = y1 - HALF_F * dy1y + sy10 * dtdy * HALF_F;
#endif
        }

        if (faceId == FACE_YMAX)
        {
          // Bottom state at top interface
          qface[ID] = r + HALF_F * dry + sr0 * dtdy * HALF_F;
          qface[IP] = p + HALF_F * dpy + sp0 * dtdy * HALF_F;
          qface[IU] = u + HALF_F * duy + su0 * dtdy * HALF_F;
          qface[IV] = v + HALF_F * dvy + sv0 * dtdy * HALF_F;
          qface[ID] = fmax(smallr, qface[ID]);
#if COP
          qface[IDY1] = y1 + HALF_F * dy1y + sy10 * dtdy * HALF_F;
#endif
        }

      } // trace_unsplit_2d_along_dir

      /**
       * This another implementation of trace computations for 2D data; it
       * is used when unsplitVersion = 1
       *
       * Note that :
       * - hydro slopes computations are done outside this routine
       *
       * \param[in]  q  primitive variable state vector
       * \param[in]  dq primitive variable slopes
       * \param[in]  dtdx dt divided by dx
       * \param[in]  dtdy dt divided by dy
       * \param[out] qm
       * \param[out] qp
       *
       */
      KOKKOS_INLINE_FUNCTION
      void trace_unsplit_hydro_2d(const HydroState &q,
                                  const HydroState &dqX,
                                  const HydroState &dqY,
                                  real_t dtdx,
                                  real_t dtdy,
                                  HydroState &qm_x,
                                  HydroState &qm_y,
                                  HydroState &qp_x,
                                  HydroState &qp_y) const
      {

        real_t gamma0 = params.settings.gamma0;
        real_t smallr = params.settings.smallr;
        real_t smallp = params.settings.smallp;

        // Cell centered values
        real_t r = q[ID];
        real_t p = q[IP];
        real_t u = q[IU];
        real_t v = q[IV];

        // Cell centered TVD slopes in X direction
        real_t drx = dqX[ID];
        drx *= HALF_F;
        real_t dpx = dqX[IP];
        dpx *= HALF_F;
        real_t dux = dqX[IU];
        dux *= HALF_F;
        real_t dvx = dqX[IV];
        dvx *= HALF_F;

        // Cell centered TVD slopes in Y direction
        real_t dry = dqY[ID];
        dry *= HALF_F;
        real_t dpy = dqY[IP];
        dpy *= HALF_F;
        real_t duy = dqY[IU];
        duy *= HALF_F;
        real_t dvy = dqY[IV];
        dvy *= HALF_F;

        // Source terms (including transverse derivatives)
        real_t sr0, su0, sv0, sp0;

        /*only true for cartesian grid */
        {
          sr0 = (-u * drx - dux * r) * dtdx + (-v * dry - dvy * r) * dtdy;
          su0 = (-u * dux - dpx / r) * dtdx + (-v * duy) * dtdy;
          sv0 = (-u * dvx) * dtdx + (-v * dvy - dpy / r) * dtdy;
          sp0 = (-u * dpx - dux * gamma0 * p) * dtdx + (-v * dpy - dvy * gamma0 * p) * dtdy;
        } // end cartesian

        // Update in time the  primitive variables
        r = r + sr0;
        u = u + su0;
        v = v + sv0;
        p = p + sp0;

        // Face averaged right state at left interface
        qp_x[ID] = r - drx;
        qp_x[IU] = u - dux;
        qp_x[IV] = v - dvx;
        qp_x[IP] = p - dpx;
        qp_x[ID] = fmax(smallr, qp_x[ID]);
        qp_x[IP] = fmax(smallp * qp_x[ID], qp_x[IP]);

        // Face averaged left state at right interface
        qm_x[ID] = r + drx;
        qm_x[IU] = u + dux;
        qm_x[IV] = v + dvx;
        qm_x[IP] = p + dpx;
        qm_x[ID] = fmax(smallr, qm_x[ID]);
        qm_x[IP] = fmax(smallp * qm_x[ID], qm_x[IP]);

        // Face averaged top state at bottom interface
        qp_y[ID] = r - dry;
        qp_y[IU] = u - duy;
        qp_y[IV] = v - dvy;
        qp_y[IP] = p - dpy;
        qp_y[ID] = fmax(smallr, qp_y[ID]);
        qp_y[IP] = fmax(smallp * qp_y[ID], qp_y[IP]);

        // Face averaged bottom state at top interface
        qm_y[ID] = r + dry;
        qm_y[IU] = u + duy;
        qm_y[IV] = v + dvy;
        qm_y[IP] = p + dpy;
        qm_y[ID] = fmax(smallr, qm_y[ID]);
        qm_y[IP] = fmax(smallp * qm_y[ID], qm_y[IP]);

      } // trace_unsplit_hydro_2d

      /**
       * Compute primitive variables slopes (dqX,dqY) for one component from q and its neighbors.
       * This routine is only used in the 2D UNSPLIT integration and slope_type = 0,1 and 2.
       *
       * Only slope_type 1 and 2 are supported.
       *
       * \param[in]  q       : current primitive variable
       * \param[in]  qPlusX  : value in the next neighbor cell along XDIR
       * \param[in]  qMinusX : value in the previous neighbor cell along XDIR
       * \param[in]  qPlusY  : value in the next neighbor cell along YDIR
       * \param[in]  qMinusY : value in the previous neighbor cell along YDIR
       * \param[out] dqX     : reference to an array returning the X slopes
       * \param[out] dqY     : reference to an array returning the Y slopes
       *
       */
      KOKKOS_INLINE_FUNCTION
      void slope_unsplit_hydro_2d_scalar(real_t q,
                                         real_t qPlusX,
                                         real_t qMinusX,
                                         real_t qPlusY,
                                         real_t qMinusY,
                                         real_t *dqX,
                                         real_t *dqY) const
      {
        real_t slope_type = params.settings.slope_type;

        real_t dlft, drgt, dcen, dsgn, slop, dlim;

        // slopes in first coordinate direction
        dlft = slope_type * (q - qMinusX);
        drgt = slope_type * (qPlusX - q);
        dcen = HALF_F * (qPlusX - qMinusX);
        dsgn = (dcen >= ZERO_F) ? ONE_F : -ONE_F;
        slop = fmin(FABS(dlft), FABS(drgt));
        dlim = slop;
        if ((dlft * drgt) <= ZERO_F)
          dlim = ZERO_F;
        *dqX = dsgn * fmin(dlim, FABS(dcen));

        // slopes in second coordinate direction
        dlft = slope_type * (q - qMinusY);
        drgt = slope_type * (qPlusY - q);
        dcen = HALF_F * (qPlusY - qMinusY);
        dsgn = (dcen >= ZERO_F) ? ONE_F : -ONE_F;
        slop = fmin(FABS(dlft), FABS(drgt));
        dlim = slop;
        if ((dlft * drgt) <= ZERO_F)
          dlim = ZERO_F;
        *dqY = dsgn * fmin(dlim, FABS(dcen));

      } // slope_unsplit_hydro_2d_scalar
      /**
       * Compute primitive variables slope (vector dq) from q and its neighbors.
       * This routine is only used in the 2D UNSPLIT integration and slope_type = 0,1 and 2.
       *
       * Only slope_type 1 and 2 are supported.
       *
       * \param[in]  q       : current primitive variable state
       * \param[in]  qPlusX  : state in the next neighbor cell along XDIR
       * \param[in]  qMinusX : state in the previous neighbor cell along XDIR
       * \param[in]  qPlusY  : state in the next neighbor cell along YDIR
       * \param[in]  qMinusY : state in the previous neighbor cell along YDIR
       * \param[out] dqX     : reference to an array returning the X slopes
       * \param[out] dqY     : reference to an array returning the Y slopes
       *
       */
      KOKKOS_INLINE_FUNCTION
      void slope_unsplit_hydro_2d(const HydroState &q,
                                  const HydroState &qPlusX,
                                  const HydroState &qMinusX,
                                  const HydroState &qPlusY,
                                  const HydroState &qMinusY,
                                  HydroState &dqX,
                                  HydroState &dqY) const
      {

        real_t slope_type = params.settings.slope_type;

        if (slope_type == 0)
        {

          dqX[ID] = ZERO_F;
          dqX[IP] = ZERO_F;
          dqX[IU] = ZERO_F;
          dqX[IV] = ZERO_F;

          dqY[ID] = ZERO_F;
          dqY[IP] = ZERO_F;
          dqY[IU] = ZERO_F;
          dqY[IV] = ZERO_F;

#if COP
          dqX[IDY1] = ZERO_F;
          dqY[IDY1] = ZERO_F;
#endif

          return;
        }

        if (slope_type == 1 || slope_type == 2)
        { // minmod or average

          slope_unsplit_hydro_2d_scalar(q[ID], qPlusX[ID], qMinusX[ID], qPlusY[ID], qMinusY[ID], &(dqX[ID]), &(dqY[ID]));
          slope_unsplit_hydro_2d_scalar(q[IP], qPlusX[IP], qMinusX[IP], qPlusY[IP], qMinusY[IP], &(dqX[IP]), &(dqY[IP]));
          slope_unsplit_hydro_2d_scalar(q[IU], qPlusX[IU], qMinusX[IU], qPlusY[IU], qMinusY[IU], &(dqX[IU]), &(dqY[IU]));
          slope_unsplit_hydro_2d_scalar(q[IV], qPlusX[IV], qMinusX[IV], qPlusY[IV], qMinusY[IV], &(dqX[IV]), &(dqY[IV]));
#if COP
          slope_unsplit_hydro_2d_scalar(q[IDY1], qPlusX[IDY1], qMinusX[IDY1], qPlusY[IDY1], qMinusY[IDY1], &(dqX[IDY1]), &(dqY[IDY1]));
#endif
        } // end slope_type == 1 or 2

      } // slope_unsplit_hydro_2d

    }; // class HydroBaseFunctor2D

  } // namespace muscl

} // namespace euler_kokkos

#endif // HYDRO_BASE_FUNCTOR_2D_H_
