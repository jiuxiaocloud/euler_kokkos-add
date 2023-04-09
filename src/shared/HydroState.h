#ifndef HYDRO_STATE_H_
#define HYDRO_STATE_H_

#include "real_type.h"

constexpr int HYDRO_2D_NBVAR = COP ? 5 : 4;
#if COP
constexpr int NUM_SPECIES = HYDRO_2D_NBVAR - 4 + 1;
constexpr int HYDRO_NBVAR_EXCOP = HYDRO_2D_NBVAR - NUM_SPECIES + 1;
#endif
constexpr int HYDRO_3D_NBVAR = 5;
constexpr int MHD_2D_NBVAR = 8;
constexpr int MHD_3D_NBVAR = 8;
constexpr int MHD_NBVAR = 8;

using HydroState2d = Kokkos::Array<real_t, HYDRO_2D_NBVAR>;
using DataArrayGPUt = Kokkos::Array<real_t, 2>;
using DataArrayGPUf = Kokkos::Array<real_t, 4>;
using HydroState3d = Kokkos::Array<real_t, HYDRO_3D_NBVAR>;
using MHDState = Kokkos::Array<real_t, MHD_NBVAR>;
using BField = Kokkos::Array<real_t, 3>;
using DataArrayspecies = Kokkos::Array<real_t, 8>;
using HydroEigen2d = Kokkos::Array<real_t, HYDRO_2D_NBVAR * HYDRO_2D_NBVAR>;
using HydroSpecies = Kokkos::Array<real_t, NUM_SPECIES>;
using HydroSpecies2 = Kokkos::Array<real_t, NUM_SPECIES * NUM_SPECIES>;

#endif // HYDRO_STATE_H_
