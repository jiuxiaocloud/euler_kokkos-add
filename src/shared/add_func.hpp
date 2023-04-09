#pragma once

#include "shared/kokkos_shared.h"
//#include "shared/HydroParams.h"
//#include "shared/HydroState.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using HydroState = HydroState2d;

const real_t kb = 1.3806 * 1.0e-23;                  //玻尔兹曼常数
const real_t pi = 3.1415926535897932384626433832795; //圆周率
const real_t NA = 6.02214128 * 1.0e23;
const real_t Ru = 8.314; //

const double Reference_params[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
#define Thermo 1
const real_t Tref = 1.0;

KOKKOS_INLINE_FUNCTION
real_t get_OmegaMiui(real_t speciesi1, const real_t T)
{
    real_t A = 1.16145;
    real_t B = -0.14874;
    real_t C = 0.52487;
    real_t D = -0.7732;
    real_t E = 2.16178;
    real_t F = -2.43787;
    real_t Tstar = T * 1.0 / speciesi1;
    return A * std::pow(Tstar, B) + C * std::exp(D * Tstar) + E * std::exp(F * Tstar);
    // printf("Tstar=%lf\n", Tstar);
}

KOKKOS_INLINE_FUNCTION
void get_Miui(HydroSpecies &Miui, const DataArrayScalar &species, const real_t T)
{
    HydroSpecies OmegaMiui;
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {
        OmegaMiui[n] = get_OmegaMiui(species(n, 1), T);
        Miui[n] = 2.6693 * 1.0e-6 * std::sqrt(species(n, 6) * 1000 * T) / (OmegaMiui[n] * species(n, 2) * species(n, 2)); //国际制单位结果
        // printf("OmegaMiui[n]=%lf,Miui[n]=%lf,d=%lf\n", OmegaMiui[n], Miui[n], species(n, 2));
    }
}

KOKKOS_INLINE_FUNCTION
real_t get_CopMiu(const DataArrayScalar &species, const HydroSpecies &Yi, const real_t T)
{
    real_t BelowMiu = 0.0;
    real_t OverMiu = 0.0;
    HydroSpecies Miui;
    get_Miui(Miui, species, T);
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {
        BelowMiu += Yi[n] * std::pow(species(n, 6), 0.5);
        OverMiu += Miui[n] * Yi[n] * std::pow(species(n, 6), 0.5);
        // printf("Miui[%d]=%lf,Yi[%d]=%lf,species(%d,6)=%lf\n", n, Miui[n], n, Yi[n], n, species(n, 6));
    }
    return OverMiu / BelowMiu;
}

KOKKOS_INLINE_FUNCTION
real_t get_OmegaDij(real_t &speciesi1, real_t &speciesj1, const real_t T)
{
    real_t A = 1.06036;
    real_t B = -0.1561;
    real_t C = 0.19300;
    real_t D = -0.47635;
    real_t E = 1.03587;
    real_t F = -1.52996;
    real_t G = 1.76474;
    real_t H = -3.89411;
    real_t Tstar = T / std::sqrt(speciesi1 * speciesj1);
    return A * std::pow(Tstar, B) + C * std::exp(D * Tstar) + E * std::exp(F * Tstar) + G * exp(H * Tstar);
}

KOKKOS_INLINE_FUNCTION
void get_CopDij(HydroSpecies2 &CopDij, const DataArrayScalar &species, real_t P, const real_t T)
{
    for (size_t i = 0; i < NUM_SPECIES; i++)
    {
        for (size_t j = 0; j < NUM_SPECIES; j++)
        {
            real_t Mij = 2.0 / (1.0 / species(i, 6) + 1.0 / species(j, 6));
            real_t sigmaij = 0.5 * (species(i, 2) + species(j, 2));
            real_t OmegaDij = get_OmegaDij(species(i, 1), species(j, 1), T);
            CopDij[i * NUM_SPECIES + j] = 0.0266 / OmegaDij * std::pow(T, 1.5) / P / std::sqrt(Mij * 1000) / sigmaij / sigmaij;
            // NOTE:P是压力
            // printf("CopDij=%lf\n", CopDij[i * NUM_SPECIES + j]);
        }
    }
}

KOKKOS_INLINE_FUNCTION
void get_CopDim(HydroSpecies &Dim, const DataArrayScalar &species, const HydroSpecies &Yi, real_t P, const real_t T)
{
    HydroSpecies2 CopDij;
    get_CopDij(CopDij, species, P, T);
    HydroSpecies Xij;
    real_t BelowXij = 0.0;
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {
        BelowXij += fmax(Yi[n] * species(n, 6), 1.0e-10);
        // printf("Yi=%lf,Wi=%lf\n", Yi[n], species(n, 6));
    }
    for (size_t m = 0; m < NUM_SPECIES; m++)
    {
        Xij[m] = Yi[m] * species(m, 6) / BelowXij;
        // printf("Xij[m]=%lf", Xij[m]);
    }
    real_t Sum_XjDij = 0.0;
    for (size_t i = 0; i < NUM_SPECIES; i++)
    {
        for (size_t j = 0; j < NUM_SPECIES; j++)
        {
            if (j != i)
            {
                Sum_XjDij += Xij[j] / fmax(CopDij[i * NUM_SPECIES + j], 1.0e-10);
                // printf("Dij=%lf\n", CopDij[i * NUM_SPECIES + j]);
            }
        }
    }
    for (size_t ii = 0; ii < NUM_SPECIES; ii++)
    {
        Dim[ii] = fmax((1 - Xij[ii]) / Sum_XjDij, 1.0e-10);
    }
}

// NOTE:粘性项在不同主方向上的ux和ux是不一样的，如在计算x方向通量的ux时使用网格点u直接计算，但是在y方向通量时要使用周围网格点的ux计算
// NOTE:所以本函数中的x仅表示主方向的意思，在计算不同方向通量时取对应方向上的变量值
KOKKOS_INLINE_FUNCTION
void add_vflux(const HydroState2d &qleft,
               const HydroState2d &q,
               const HydroState2d &qright1,
               const HydroState2d &qright2,
               const HydroParams &params,
               const DataArrayGPUf &uy,
               const DataArrayGPUf &vy,
               HydroState &flux,
               const real_t dx,
               const int i,
               const int j)
{
    HydroSpecies Yi, Yi0, Yi1, Yi2, Yi3, Yix, Dim, CopVflux;
    Yi[0] = 1.0;
    Yi0[0] = 1.0;
    Yi1[0] = 1.0;
    Yi2[0] = 1.0;
    Yi3[0] = 1.0;
    for (size_t n = 1; n < NUM_SPECIES; n++)
    {
        Yi0[n] = qleft[IDY1 + n - 1];                                          // fmax(fmin(qleft[IDY1 + n - 1], 1.0), 1e-10);
        Yi1[n] = q[IDY1 + n - 1];                                              // fmax(fmin(q[IDY1 + n - 1], 1.0), 1e-10);
        Yi2[n] = qright1[IDY1 + n - 1];                                        // TODO:赋值有问题                                        // fmax(fmin(qright1[IDY1 + n - 1], 1.0), 1e-10);
        Yi3[n] = qright2[IDY1 + n - 1];                                        // fmax(fmin(qright2[IDY1 + n - 1], 1.0), 1e-10);
        Yi[n] = 1.0 / 16.0 * (-Yi0[n] + 9.0 * Yi1[n] + 9.0 * Yi2[n] - Yi3[n]); // fmax(fmin(1.0 / 16.0 * (-Yi0[n] + 9.0 * Yi1[n] + 9.0 * Yi2[n] - Yi3[n]), 1.0), 1e-10); // 1.0 / 16.0 * (-Yi0[n] + 9.0 * Yi1[n] + 9.0 * Yi2[n] - Yi3[n]); //
        Yix[n] = 1.0 / 24.0 / dx * (Yi0[n] - 27.0 * Yi1[n] + 27.0 * Yi2[n] - Yi3[n]);
        Yi[0] = Yi[0] - Yi[n];    // fmax(fmin(Yi[0] - Yi[n], 1.0), 1e-10);
        Yi0[0] = Yi0[0] - Yi0[n]; // fmax(fmin(Yi0[0] - Yi0[n], 1.0), 1e-10);
        Yi1[0] = Yi1[0] - Yi1[n]; // fmax(fmin(Yi1[0] - Yi1[n], 1.0), 1e-10);
        Yi2[0] = Yi2[0] - Yi2[n]; // fmax(fmin(Yi2[0] - Yi2[n], 1.0), 1e-10);
        Yi3[0] = Yi3[0] - Yi3[n]; // fmax(fmin(Yi3[0] - Yi3[n], 1.0), 1e-10);
        // printf("Yi0[n]=%lf,Yi1[n]=%lf,Yi2[n]=%lf,Yi3[n]=%lf\n", Yi0[n], Yi1[n], Yi2[n], Yi3[n]);
    }
    Yix[0] = 1.0 / 24.0 / dx * (Yi0[0] - 27.0 * Yi1[0] + 27.0 * Yi2[0] - Yi3[0]);
    real_t P_ = 1.0 / 16.0 * (-qleft[IP] + 9.0 * q[IP] + 9.0 * qright1[IP] - qright2[IP]);
    const real_t T = 298.15; // p / r / Ru / ((1 - y1) / params.species(0, 6) + y1 / params.species(1, 6)); //单位K
    get_CopDim(Dim, params.species, Yi, P_, T);
    real_t Sum_DimYix = 0.0;
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {
        Sum_DimYix += Dim[n] * Yix[n];
    }
    // CopVflux[1] = Dim[1] * Yix[1] - Yi[1] * Sum_DimYix;
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {
        CopVflux[n] = Dim[n] * Yix[n] - Yi[n] * Sum_DimYix;
        // printf(" CopVflux[%d]=%lf \n", n, CopVflux[n]);
    }
    real_t miu_ = get_CopMiu(params.species, Yi, T);
    real_t rho_ = 1.0 / 16.0 * (-qleft[ID] + 9.0 * q[ID] + 9.0 * qright1[ID] - qright2[ID]);
    real_t u_ = 1.0 / 16.0 * (-qleft[IU] + 9.0 * q[IU] + 9.0 * qright1[IU] - qright2[IU]);
    real_t v_ = 1.0 / 16.0 * (-qleft[IV] + 9.0 * q[IV] + 9.0 * qright1[IV] - qright2[IV]);
    real_t ux_ = 1.0 / 24.0 / dx * (qleft[IU] - 27.0 * q[IU] + 27.0 * qright1[IU] - qright2[IU]);
    real_t vx_ = 1.0 / 24.0 / dx * (qleft[IV] - 27.0 * q[IV] + 27.0 * qright1[IV] - qright2[IV]);
    real_t uy_ = 1.0 / 16.0 * (-uy[0] + 9.0 * uy[1] + 9.0 * uy[2] - uy[3]);
    real_t vy_ = 1.0 / 16.0 * (-vy[0] + 9.0 * vy[1] + 9.0 * vy[2] - vy[3]);
    real_t tao11 = 2.0 / 3.0 * miu_ * (2.0 * ux_ - vy_);
    real_t tao12 = miu_ * (uy_ + vx_);

    flux[IP] += -(u_ * tao11 + v_ * tao12);
    flux[IU] += -tao11;
    flux[IV] += -tao12;
    flux[IDY1] += -rho_ * CopVflux[1];
    if (params.if_debug)
    {
        for (size_t n = 0; n < NUM_SPECIES; n++)
        {
            printf("网格点(%d,%d),Yi0[%d]=%lf,Yi1[%d]=%lf,Yi2[%d]=%lf,Yi3[%d]=%lf,Yi[%d]=%lf,Yix[%d]=%lf\n", i, j, n, Yi0[n], n, Yi1[n], n, Yi2[n], n, Yi3[n], n, Yi[n], n, Yix[n]);
        }
    }
}

KOKKOS_INLINE_FUNCTION
void get_Cpi(const HydroParams &params, HydroSpecies &Cpi, real_t T) // ii用来区分是哪个气体的对应属性，后期可以加一个组分名字的字符串数组
{
    HydroSpecies Ri;
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {                                      // NOTE:注意由Hia和Hib计算得到的Cp单位是J/kg/K
        Ri[n] = Ru / params.species(n, 6); // NOTE:注意Ri的单位，与Wi的单位有关
        if (T < (200.0 / Tref))
            T = 200.0 / Tref;
#if Thermo
        if (T >= (1000.0 / Tref) && T < (6000.0 / Tref))
            Cpi[n] = Ri[n] * (params.Hia(n, 0, 1) / T / T + params.Hia(n, 1, 1) / T + params.Hia(n, 2, 1) + params.Hia(n, 3, 1) * T + params.Hia(n, 4, 1) * T * T + params.Hia(n, 5, 1) * std::pow(T, 3) + params.Hia(n, 6, 1) * std::pow(T, 4));
        else if (T < (1000.0 / Tref))
            Cpi[n] = Ri[n] * (params.Hia(n, 0, 0) / T / T + params.Hia(n, 1, 0) / T + params.Hia(n, 2, 0) + params.Hia(n, 3, 0) * T + params.Hia(n, 4, 0) * T * T + params.Hia(n, 5, 0) * std::pow(T, 3) + params.Hia(n, 6, 0) * std::pow(T, 4));
        else if (T >= (6000.0 / Tref) && T < (15000.0 / Tref))
        {
            Cpi[n] = Ri[n] * (params.Hia(n, 0, 2) / T / T + params.Hia(n, 1, 2) / T + params.Hia(n, 2, 2) + params.Hia(n, 3, 2) * T + params.Hia(n, 4, 2) * T * T + params.Hia(n, 5, 2) * std::pow(T, 3) + params.Hia(n, 6, 2) * std::pow(T, 4));
        }
        else
        {
            printf(" T > 15000 K,please check!!!NO Cpi[n) for T>15000 K");
        }
#else
        // Cpi[n)/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
        if (T > (1000.0 / Tref))
            Cpi[n] = Ri * (params.Hia(n, 0, 0) + params.Hia(n, 1, 0) * T + params.Hia(n, 2, 0) * T * T + params.Hia(n, 3, 0) * T * T * T + params.Hia(n, 4, 0) * T * T * T * T);
        else
            Cpi[n] = Ri * (params.Hia(n, 0, 1) + params.Hia(n, 1, 1) * T + params.Hia(n, 2, 1) * T * T + params.Hia(n, 3, 1) * T * T * T + params.Hia(n, 4, 1) * T * T * T * T);
#endif
    }
}

KOKKOS_INLINE_FUNCTION
void get_Hi(const HydroParams &params, HydroSpecies &hi, real_t T)
{
    // NOTE：温度的无量纲，对温度无量纲化之后参数a也要进行无量纲化
    real_t TT = 30000.0 / Tref;
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {
        hi[n] = 0.0;
        real_t Ri = Ru / params.species(n, 6); // NOTE：Non_dim_Ri=Ri*Tref/rs.u_ref/rs.u_ref
#if Thermo
        // NOTE：Non_dim of Hia && Hib, only for h&Cp not for S ATTENTATION
        //  200K~1000K
        params.Hia(n, 0, 0) = params.Hia(n, 0, 0) / std::pow(Tref, 2);
        params.Hia(n, 1, 0) = params.Hia(n, 1, 0) / Tref;

        params.Hia(n, 3, 0) = params.Hia(n, 3, 0) * Tref;
        params.Hia(n, 4, 0) = params.Hia(n, 4, 0) * std::pow(Tref, 2);
        params.Hia(n, 5, 0) = params.Hia(n, 5, 0) * std::pow(Tref, 3);
        params.Hia(n, 6, 0) = params.Hia(n, 6, 0) * std::pow(Tref, 4);
        params.Hib(n, 0, 0) = params.Hib(n, 0, 0) / Tref + params.Hia(n, 1, 0) * std::log(Tref);
        // 1000K~6000K
        params.Hia(n, 0, 1) = params.Hia(n, 0, 1) / std::pow(Tref, 2);
        params.Hia(n, 1, 1) = params.Hia(n, 1, 1) / Tref;

        params.Hia(n, 3, 1) = params.Hia(n, 3, 1) * Tref;
        params.Hia(n, 4, 1) = params.Hia(n, 4, 1) * std::pow(Tref, 2);
        params.Hia(n, 5, 1) = params.Hia(n, 5, 1) * std::pow(Tref, 3);
        params.Hia(n, 6, 1) = params.Hia(n, 6, 1) * std::pow(Tref, 4);
        params.Hib(n, 0, 1) = params.Hib(n, 0, 1) / Tref + params.Hia(n, 1, 1) * std::log(Tref);
        // 6000K~15000K
        params.Hia(n, 0, 2) = params.Hia(n, 0, 2) / std::pow(Tref, 2);
        params.Hia(n, 1, 2) = params.Hia(n, 1, 2) / Tref;

        params.Hia(n, 3, 2) = params.Hia(n, 3, 2) * Tref;
        params.Hia(n, 4, 2) = params.Hia(n, 4, 2) * std::pow(Tref, 2);
        params.Hia(n, 5, 2) = params.Hia(n, 5, 2) * std::pow(Tref, 3);
        params.Hia(n, 6, 2) = params.Hia(n, 6, 2) * std::pow(Tref, 4);
        params.Hib(n, 0, 2) = params.Hib(n, 0, 2) / Tref + params.Hia(n, 1, 2) * std::log(Tref);
#else
        params.Hia(n, 1, 0) = params.Hia(n, 1, 0) * Tref;
        params.Hia(n, 2, 0) = params.Hia(n, 2, 0) * std::pow(Tref, 2);
        params.Hia(n, 3, 0) = params.Hia(n, 3, 0) * std::pow(Tref, 3);
        params.Hia(n, 4, 0) = params.Hia(n, 4, 0) * std::pow(Tref, 4);
        params.Hia(n, 5, 0) = params.Hia(n, 5, 0) / Tref;

        params.Hia(n, 1, 1) = params.Hia(n, 1, 1) * Tref;
        params.Hia(n, 2, 1) = params.Hia(n, 2, 1) * std::pow(Tref, 2);
        params.Hia(n, 3, 1) = params.Hia(n, 3, 1) * std::pow(Tref, 3);
        params.Hia(n, 4, 1) = params.Hia(n, 4, 1) * std::pow(Tref, 4);
        params.Hia(n, 5, 1) = params.Hia(n, 5, 1) / Tref;
#endif
        if (T < 200.0 / Tref)
        {
            TT = T;
            T = 200.0 / Tref;
        }
#if Thermo
        if (T >= (1000.0 / Tref) && T < (6000.0 / Tref))
            hi[n] = Ri * (-params.Hia(n, 0, 1) * 1.0 / T + params.Hia(n, 1, 1) * std::log(T) + params.Hia(n, 2, 1) * T + 0.5 * params.Hia(n, 3, 1) * T * T + params.Hia(n, 4, 1) * std::pow(T, 3) / 3.0 + params.Hia(n, 5, 1) * std::pow(T, 4) / 4.0 + params.Hia(n, 6, 1) * std::pow(T, 5) / 5.0 + params.Hib(n, 0, 1));
        else if (T < (1000.0 / Tref))
            hi[n] = Ri * (-params.Hia(n, 0, 0) * 1.0 / T + params.Hia(n, 1, 0) * std::log(T) + params.Hia(n, 2, 0) * T + 0.5 * params.Hia(n, 3, 0) * T * T + params.Hia(n, 4, 0) * std::pow(T, 3) / 3.0 + params.Hia(n, 5, 0) * std::pow(T, 4) / 4.0 + params.Hia(n, 6, 0) * std::pow(T, 5) / 5.0 + params.Hib(n, 0, 0));
        else if (T >= (6000.0 / Tref) && T < (15000.0 / Tref))
        {
            hi[n] = Ri * (-params.Hia(n, 0, 2) * 1.0 / T + params.Hia(n, 1, 2) * std::log(T) + params.Hia(n, 2, 2) * T + 0.5 * params.Hia(n, 3, 2) * T * T + params.Hia(n, 4, 2) * std::pow(T, 3) / 3.0 + params.Hia(n, 5, 2) * std::pow(T, 4) / 4.0 + params.Hia(n, 6, 2) * std::pow(T, 5) / 5.0 + params.Hib(n, 0, 2));
        }
        else
        {
            printf("T > 15000 K,please check!!!NO h for T>15000 K");
        }
#else
        // H/RT = a1 + a2/2*T + a3/3*T^2 + a4/4*T^3 + a5/5*T^4 + a6/T
        if (T > (1000.0 / Tref))
            hi[n] = Ri * (params.Hia(n, 0, 0) * T + params.Hia(n, 1, 0) * T * T / 2.0 + params.Hia(n, 2, 0) * T * T * T / 3.0 + params.Hia(n, 3, 0) * T * T * T * T / 4.0 + params.Hia(n, 4, 0) * T * T * T * T * T / 5.0 + params.Hia(n, 5, 0));
        else
            hi[n] = Ri * (params.Hia(n, 0, 1) * T + params.Hia(n, 1, 1) * T * T / 2.0 + params.Hia(n, 2, 1) * T * T * T / 3.0 + params.Hia(n, 3, 1) * T * T * T * T / 4.0 + params.Hia(n, 4, 1) * T * T * T * T * T / 5.0 + params.Hia(n, 5, 1));
#endif
        if (TT < 200.0 / Tref)
        {
            HydroSpecies Cpi;
            get_Cpi(params, Cpi, 200.0 / Tref);
            hi[n] += Cpi[n] * (TT - 200.0 / Tref);
        }
    }
}

KOKKOS_INLINE_FUNCTION
void get_CopCp(const HydroParams &params, const HydroSpecies &Yi, const real_t T, real_t *CopCp)
{
    HydroSpecies Cpi;
    get_Cpi(params, Cpi, T);
    *CopCp = 0.0;
    for (size_t i = 0; i < NUM_SPECIES; i++)
    {
        *CopCp += Yi[i] * Cpi[i];
    }
    return;
}

KOKKOS_INLINE_FUNCTION
void get_CopW(const DataArrayScalar &species, const HydroSpecies &Yi, real_t *CopW)
{
    real_t _W = 0.0;
    for (size_t i = 0; i < NUM_SPECIES; i++)
    {
        _W += Yi[i] / species(i, 6);
    }
    *CopW = 1.0 / _W;
    return;
}

KOKKOS_INLINE_FUNCTION
void get_CopGamma(const HydroParams &params, const HydroSpecies &Yi, const real_t T, real_t *CopGamma)
{
    real_t Cp;
    get_CopCp(params, Yi, T, &Cp); // NOTE：混合物的定压比热容Cp
    real_t CopW;
    get_CopW(params.species, Yi, &CopW); // NOTE:混合物的摩尔质量
    *CopGamma = Cp / (Cp - Ru / CopW);
    if (*CopGamma > 1)
        return;
    else
    {
        printf("Gamma计算错误1\n");
        printf("CopGamma=%lf,Yi of qloc =%lf,%lf,Cp=%lf,CopW=%lf\n", *CopGamma, Yi[0], Yi[1], Cp, CopW);
        // abort();
    }
}

KOKKOS_INLINE_FUNCTION
void get_CopC2(const real_t T, const HydroParams &params, HydroState &qLoc, real_t *CopC2)
{
    real_t D = qLoc[ID]; // NOTE:qleft=i,qright=i+1,这里的left和right相对于i+1/2来说。
    real_t u = qLoc[IU];
    real_t v = qLoc[IV];
    real_t w = 0.0;
    real_t q2 = u * u + v * v + w * w;
    real_t p = qLoc[IP];
    HydroSpecies Yi = {1 - qLoc[IDY1], qLoc[IDY1]};
    real_t gamma;
    // printf("Yi of qLoc=%lf,%lf", Yi[0], Yi[1]);
    get_CopGamma(params, Yi, T, &gamma);
    real_t h = (p / (gamma - 1.0) + p) / qLoc[ID];
    real_t H = h + 0.5 * q2;
    real_t e = p / (gamma - 1.0) / qLoc[ID];
    // NOTE:change Yi[],zi[],get c;
    real_t Sum_dpdrhoi = 0.0; // Sum_dpdrhoi:first of c2,存在累加项
    HydroSpecies hi, Ri, _dpdrhoi;
    // enthalpy h_i (unit: J/kg) and mass fraction Y_i of each specie
    get_Hi(params, hi, T);
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {
        Ri[n] = Ru / params.species(n, 6);
        _dpdrhoi[n] = (gamma - 1.0) * (hi[0] - hi[n]) + gamma * (Ri[n] - Ri[0]) * T;
        if (0 != n)
            Sum_dpdrhoi += Yi[n] * _dpdrhoi[n];
    }
    *CopC2 = Sum_dpdrhoi + (gamma - 1) * (p / D + e - hi[0]) + gamma * Ri[0] * T;
    return;
    //  derivatives
}

KOKKOS_INLINE_FUNCTION
void get_CopGamma0(const real_t T, const HydroParams &params, const HydroState &qleft, const HydroState &qright, real_t *CopGamma)
{
    real_t D = sqrt(qright[ID] / qleft[ID]); // NOTE:qleft=i,qright=i+1,这里的left和right相对于i+1/2来说。
    real_t D1 = 1.0 / (D + 1.0);
    HydroSpecies Yi = {1 - (qleft[IDY1] + D * qright[IDY1]) * D1, (qleft[IDY1] + D * qright[IDY1]) * D1};
    real_t Cp;
    get_CopCp(params, Yi, T, &Cp); // NOTE：混合物的定压比热容Cp
    real_t CopW;
    get_CopW(params.species, Yi, &CopW); // NOTE:混合物的摩尔质量
    *CopGamma = Cp / (Cp - Ru / CopW);
    if (*CopGamma > 1)
        return;
    else
    {
        // printf("Yi=%lf,%lf,D=%lf,D1=%lf,Cp=%lf,CopW=%lf,Gamma计算错误2", Yi[0], Yi[1], D, D1, Cp，CopW);
        //  abort();
    }
}

KOKKOS_INLINE_FUNCTION
void get_CopC2(const real_t T, const HydroParams &params, HydroState &qleft, HydroState &qright, real_t *CopC2)
{
    real_t D = sqrt(qright[ID] / qleft[ID]); // NOTE:qleft=i,qright=i+1,这里的left和right相对于i+1/2来说。
    real_t D1 = 1.0 / (D + 1.0);
    real_t u = (qleft[IU] + D * qright[IU]) * D1;
    real_t v = (qleft[IU] + D * qright[IU]) * D1;
    real_t w = 0.0;
    real_t q2 = u * u + v * v + w * w;
    real_t pr = qright[IP];
    real_t pl = qleft[IP];
    real_t p = (pl + D * pr) * D1;
    HydroSpecies Yi = {1 - (qleft[IDY1] + D * qright[IDY1]) * D1, (qleft[IDY1] + D * qright[IDY1]) * D1};
    real_t gamma;
    get_CopGamma(params, Yi, T, &gamma);
    real_t hr = (pr / (gamma - 1.0) + pr) / qright[ID];
    real_t hl = (pl / (gamma - 1.0) + pl) / qleft[ID];
    real_t er = (pr / (gamma - 1.0)) / qright[ID];
    real_t el = (pl / (gamma - 1.0)) / qleft[ID];
    real_t Hr = hr + 0.5 * (qright[IU] * qright[IU] + qright[IV] * qright[IV]);
    real_t Hl = hl + 0.5 * (qleft[IU] * qleft[IU] + qleft[IV] * qleft[IV]);
    real_t H = (Hl + D * Hr) * D1;
    real_t h = (hl + D * hr) * D1;
    real_t e = (el + D * er) * D1;
    // NOTE:change Yi[],zi[],get c;
    real_t Sum_dpdrhoi = 0.0; // Sum_dpdrhoi:first of c2,存在累加项
    HydroSpecies hi, hi_L, hi_R, Ri, _dpdrhoi;
    // enthalpy h_i (unit: J/kg) and mass fraction Y_i of each specie
    get_Hi(params, hi_L, T);
    get_Hi(params, hi_R, T);
    for (size_t n = 0; n < NUM_SPECIES; n++)
    {
        hi[n] = (hi_L[n] + D * hi_R[n]) * D1;
        Ri[n] = Ru / params.species(n, 6);
        _dpdrhoi[n] = (gamma - 1.0) * (hi[0] - hi[n]) + gamma * (Ri[n] - Ri[0]) * T;
        if (0 != n)
            Sum_dpdrhoi += Yi[n] * _dpdrhoi[n];
    }
    *CopC2 = Sum_dpdrhoi + (gamma - 1) * (p / D + e - hi[0]) + gamma * Ri[0] * T;
    return;
    // NOTE:不能把原始变量的均值都算出来以后进行统一计算进阶的部分
    //  derivatives
}

// KOKKOS_INLINE_FUNCTION void
// calculate_dqa(real_t *dq, real_t qPlus, real_t qMinus, real_t qPlus2, real_t qMinus2, real_t dx)
// {
//     *dq = (8.0 * (qPlus - qMinus) - (qPlus2 - qMinus2)) / 12.0 / dx;
//     //*dq = (qPlus - *dq) / dx;
// }

// KOKKOS_INLINE_FUNCTION
// void GetQuadraticInterCoeff(BField T_s, BField Omega_t, BField aa)
// {
//     aa[2] = ((Omega_t[0] - Omega_t[1]) / (T_s[0] - T_s[1]) - (Omega_t[1] - Omega_t[2]) / (T_s[1] - T_s[2])) / (T_s[0] - T_s[2]);
//     aa[1] = (Omega_t[0] - Omega_t[1]) / (T_s[0] - T_s[1]) - aa[2] * (T_s[0] + T_s[1]);
//     aa[0] = Omega_t[0] - aa[1] * T_s[0] - aa[2] * T_s[0] * T_s[0];
// }

// KOKKOS_INLINE_FUNCTION
// real_t Omega_interpolated(const DataArray2d &Omega_table, const DataArrayOne &T_star, const DataArrayOne &delta_star,
//                           const real_t Tstar, const real_t deltastar, int index)
// { // NOTE：Omega form LYX
//     int ti1, ti2, ti3;
//     if (Tstar > T_star(0) && Tstar < T_star(36))
//     {
//         int ii = 1;
//         {
//             while (Tstar > T_star(ii))
//                 ii = ii + 1;
//         }
//         ti1 = ii - 1;
//         ti2 = ii;
//         ti3 = ii + 1;
//     }
//     else if (Tstar <= T_star(0))
//     {
//         ti1 = 0;
//         ti2 = 1;
//         ti3 = 2;
//     }
//     else if (Tstar >= T_star(36))
//     {
//         ti1 = 34;
//         ti2 = 35;
//         ti3 = 36;
//     }
//     int tj1, tj2, tj3;
//     if (deltastar > delta_star(0) && deltastar < delta_star(7))
//     {
//         int jj = 1;
//         {
//             while (deltastar > delta_star(jj))
//                 jj = jj + 1;
//         }
//         tj1 = jj - 1;
//         tj2 = jj;
//         tj3 = jj + 1;
//     }
//     else if (deltastar <= delta_star(0))
//     {
//         tj1 = 0;
//         tj2 = 1;
//         tj3 = 2;
//     }
//     else if (deltastar >= delta_star(7))
//     {
//         tj1 = 5;
//         tj2 = 6;
//         tj3 = 7;
//     }
//     // real_t aa(3);
//     BField T_s, delta_s, Omega_t, tmp, aa;
//     T_s[0] = T_star(ti1);
//     T_s[1] = T_star(ti2);
//     T_s[2] = T_star(ti3);
//     delta_s[0] = delta_star(tj1);
//     delta_s[1] = delta_star(tj2);
//     delta_s[2] = delta_star(tj3);
//     Omega_t[0] = Omega_table(index, ti1, tj1);
//     Omega_t[1] = Omega_table(index, ti2, tj1);
//     Omega_t[2] = Omega_table(index, ti3, tj1);

//     GetQuadraticInterCoeff(T_s, Omega_t, aa);
//     tmp[0] = aa[0] + aa[1] * Tstar + aa[2] * Tstar * Tstar;

//     Omega_t[0] = Omega_table(index, ti1, tj2);
//     Omega_t[1] = Omega_table(index, ti2, tj2);
//     Omega_t[2] = Omega_table(index, ti3, tj2);

//     GetQuadraticInterCoeff(T_s, Omega_t, aa);
//     tmp[1] = aa[0] + aa[1] * Tstar + aa[2] * Tstar * Tstar;

//     Omega_t[0] = Omega_table(index, ti1, tj3);
//     Omega_t[1] = Omega_table(index, ti2, tj3);
//     Omega_t[2] = Omega_table(index, ti3, tj3);

//     GetQuadraticInterCoeff(T_s, Omega_t, aa);
//     tmp[2] = aa[0] + aa[1] * Tstar + aa[2] * Tstar * Tstar;

//     GetQuadraticInterCoeff(delta_s, tmp, aa);

//     return aa[0] + aa[1] * deltastar + aa[2] * deltastar * deltastar;
// }

/**
 * @brief get binary(specie j&specie k) diffusion coefficient at T temperature per pressure
 * @para T temperature
 * @para specie
 * unit: 1 pa=1 kg/(m.s2)=10 g/(cm.s2)
   [Wi]=g/mol;   [T]=K;  [PP]=pa;   [Djk]=cm2/s
 */
// KOKKOS_INLINE_FUNCTION
// real_t Dkj(const real_t T, real_t PP,
//            DataArrayspecies specie_k, DataArrayspecies specie_j,
//            const DataArray2d &Omega_table, const DataArrayOne &T_star, const DataArrayOne &delta_star) // PP:pressure,unit:Pa
// {                                                                                                      // NOTE：Omega是查表所得
//     real_t epsilon_jk_kB, d_jk, mue_jk_sqr;
//     // either both nonpolar or both polar
//     if ((specie_j[3] > 0 && specie_k[3] > 0) || (specie_j[3] == 0 && specie_k[3] == 0))
//     {
//         epsilon_jk_kB = std::sqrt(specie_j[1] * specie_k[1]); // unit:K,equation5-6
//         d_jk = (specie_j[2] + specie_k[2]) / 2.0;
//         mue_jk_sqr = specie_j[3] * specie_k[3];
//     }
//     // polar molecule interacting with a nonpolar molecule
//     else
//     {
//         real_t epsilon_n_kB, epsilon_p_kB, alpha_n, mue_p, d_n, d_p; // equation 5-9~~5-14
//         if (specie_k[3] > 0 && specie_j[3] == 0)
//         {
//             epsilon_n_kB = specie_j[1];
//             epsilon_p_kB = specie_k[1];
//             alpha_n = specie_j[4];
//             d_n = specie_j[2];
//             d_p = specie_k[2];
//             mue_p = specie_k[3];
//         }
//         if (specie_j[3] > 0 && specie_k[3] == 0)
//         {
//             epsilon_n_kB = specie_k[1];
//             epsilon_p_kB = specie_j[1];
//             alpha_n = specie_k[4];
//             d_n = specie_k[2];
//             d_p = specie_j[2];
//             mue_p = specie_j[3];
//         }
//         real_t alpha_n_star = alpha_n / std::pow(d_n, 3);                                           // equation5-13
//         real_t mue_p_star = mue_p / std::pow(epsilon_p_kB * kb, 0.5) / std::pow(d_p, 1.5) * 1.0e-6; // equation5-14
//         real_t ksi = 1 + 0.25 * alpha_n_star * mue_p_star * std::sqrt(epsilon_p_kB / epsilon_n_kB); // equation5-12

//         epsilon_jk_kB = ksi * ksi * std::sqrt(epsilon_n_kB * epsilon_p_kB); // equation5-9
//         d_jk = std::pow(ksi, -1.0 / 6.0) * (specie_j[2] + specie_k[2]) / 2.0;
//         mue_jk_sqr = 0.0;
//     }
//     real_t T_jk_star = T / epsilon_jk_kB;                                                        // equation5-15
//     real_t delta_jk_star = 0.5 * mue_jk_sqr / d_jk / d_jk / d_jk / epsilon_jk_kB / kb * 1.0e-12; // equation5-16
//     real_t W_jk = specie_k[6] * specie_j[6] / (specie_k[6] + specie_j[6]) / NA;                  // unit,g;equation5-5
//     // real_t W_jk=specie_k(Wi)*specie_j(Wi)/(specie_k(Wi)+specie_j(Wi));
//     real_t Omega1 = Omega_interpolated(Omega_table, T_star, delta_star, T_jk_star, delta_jk_star, 1);
//     PP = PP * 10.0;                                                                                                    // pa==>g/(cm.s2)
//     real_t temp = 3 * std::sqrt(2 * pi * std::pow(T * kb, 3) / W_jk) / (16 * PP * pi * d_jk * d_jk * Omega1) * 1.0e12; // equation5-4 //unit:m*m/s
//     return temp;                                                                                                       //
// }

// KOKKOS_INLINE_FUNCTION
// real_t PHIkj(DataArrayspecies specie_k, DataArrayspecies specie_j, const real_t T) // equation 5-50
// {                                                                            // NOTE:得到混合物粘性系数中间变量Phi
//     real_t phi = 0.0;
//     phi = std::pow(specie_j[6] / specie_k[6], 0.25) * std::pow(specie_k[7] / specie_j[7], 0.5);
//     phi = (phi + 1.0) * (phi + 1.0) * 0.5 / std::sqrt(2.0);
//     phi = phi * std::pow(1.0 + specie_k[6] / specie_j[6], -0.5);
//     return phi; //国际制单位
// }

// KOKKOS_INLINE_FUNCTION
// void get_Miui(HydroSpecies &Miu, const DataArray2d &Omega_table, const DataArrayOne &T_star, const DataArrayOne &delta_star,
//               const DataArrayspecies &species, const real_t T)
// {
//     for (size_t n = 0; n < NUM_SPECIES; n++)
//     {
//         Wi[n] = species(n, 6);                                                                                                                     // kg/mol
//         sigmai[n] = species(n, 2);                                                                                                                 // 1e-10m
//         Tstar[n] = T * 1.0 / species(n, 1);                                                                                                        // 1
//         deltastar[n] = 0.5 * species(n, 3) * species(n, 3) / species(n, 1) / (kb * 1e5) / (species(n, 2) * species(n, 2) * species(n, 2)) * 1e-12; // NOTE：从文件中读到的值我都没有修改，kb与源程序单位不一样
//         real_t Omega2[n] = Omega_interpolated(Omega_table, T_star, delta_star, Tstar[n], deltastar[n], 0);
//         Miu[n] = 0.5 / 16.0 * sqrt(pi * (Wi[n] * 1e3) * (kb * 1e5) * T / NA) / (pi * sigmai[n] * sigmai[n] * Omega2[n]); // unit: Pa.s=kg/(m.s)
//     }
// }

// KOKKOS_INLINE_FUNCTION
// real_t get_CopMiu(const DataArrayScalar &species, DataArrayspecies &species0, DataArrayspecies &species1,
//                   DataArrayGPUt &X, real_t y1, const real_t T)
// {
//     HydroSpecies miu;
//     DataArrayGPUSpecies2 phi;
//     get_Miu(miu, Omega_table, T_star, delta_star, species, T);
//     phi[0] = PHIkj(species0, species0, T);
//     phi[1] = PHIkj(species0, species1, T);
//     phi[2] = PHIkj(species1, species0, T);
//     phi[3] = PHIkj(species1, species1, T);
//     real_t tempmiu = 0;
//     if (int(y1 * 100) * 1.0 / 100 < 1e-10)
//     {
//         tempmiu = miu[0];
//     }
//     else if (int(y1 * 100) * 1.0 / 100 - 1.0 < 1e-10)
//     {
//         tempmiu = miu[1];
//     }
//     else
//     {
//         for (size_t k = 0; k < 2; k++)
//         {
//             real_t tempXPhi = 0;
//             for (size_t j = 0; j < 2; j++)
//             {
//                 tempXPhi += X[j] * phi[k * 2 + j]; // NOTE:先按一维方式来
//             }
//             tempmiu += X[k] * miu[k] / tempXPhi;
//         }
//     }
//     return tempmiu * 1e-6; //转换为国际制单位
// }

// KOKKOS_INLINE_FUNCTION
// // static void getDkm(real_t *Dkm, const Specie *specie, real_t y1, const real_t T, real_t p, Omega getOmega)
// void getDkm(const DataArrayScalar &species, const DataArray2d &Omega_table, const DataArrayOne &T_star, const DataArrayOne &delta_star,
//             real_t y1, const real_t T, real_t p,
//             DataArrayGPUt &Dkm, DataArrayGPUt &X, DataArrayGPUt &y, DataArrayGPUf &Djk, DataArrayspecies &species0, DataArrayspecies &species1,
//             BField &aa, BField &T_s, BField &Omega_t, BField &tmp, BField &delta_s)
// {
//     // DataArrayOne X("X", 2), y("y", 2);
//     // DataArrayScalar Djk("Djk", 2, 2);
//     // DataArrayOne species0("species0", 8), species1("species1", 8);
//     Djk[0] = Dkj(T, p, species0, species0, Omega_table, T_star, delta_star, aa, T_s, Omega_t, tmp, delta_s);
//     Djk[1] = Dkj(T, p, species0, species1, Omega_table, T_star, delta_star, aa, T_s, Omega_t, tmp, delta_s);
//     Djk[3] = Dkj(T, p, species1, species0, Omega_table, T_star, delta_star, aa, T_s, Omega_t, tmp, delta_s);
//     Djk[4] = Dkj(T, p, species1, species1, Omega_table, T_star, delta_star, aa, T_s, Omega_t, tmp, delta_s);
//     for (size_t k = 0; k < 2; k++)
//     {
//         real_t tempbelow = 0;
//         for (size_t j = 0; j < 2; j++)
//         {
//             if (j != k)
//             {
//                 tempbelow += X[j] / Djk[j * 2 + k];
//             }
//         }
//         Dkm[k] = (1 - y[k]) / tempbelow * 1e8; //国际制单位
//     }
// }

// KOKKOS_INLINE_FUNCTION
// void get_MatrixEigen(const HydroState &qleft, const HydroState &qright, const HydroParams &params,
//                      HydroEigen2d eigen_l, HydroEigen2d eigen_r, HydroState2d eigen_value, int A)
// {
//     real_t D = sqrt(qright[ID] / qleft[ID]); // NOTE:qleft=i,qright=i+1,这里的left和right相对于i+1/2来说。
//     real_t D1 = 1.0 / (D + 1.0);
//     real_t u = (qleft[IU] + D * qright[IU]) * D1;
//     real_t v = (qleft[IU] + D * qright[IU]) * D1;
//     real_t w = 0.0;
//     int B = 1 - A; // flux_x when A==1,flux_y when A==0;
//     BField _uvw = {u, v, w};
//     {
//         real_t q2 = u * u + v * v + w * w;
//         real_t u_ = A * u + B * v;
//         real_t v_ = A * v - B * u;
//         real_t pr = qright[IP];
//         real_t pl = qleft[IP];
//         real_t hr = (pr / (gamma - 1.0) + pr) / qright[ID];
//         real_t hl = (pl / (gamma - 1.0) + pl) / qleft[ID];
//         real_t Hr = hr + 0.5 * (qright[IU] * qright[IU] + qright[IV] * qright[IV]);
//         real_t Hl = hl + 0.5 * (qleft[IU] * qleft[IU] + qleft[IV] * qleft[IV]);
//         real_t H = (Hl + D * Hr) * D1;
//         real_t h = (hl + D * hr) * D1;
// #if Mixture
//         real_t c2 = get_SoundSpeedMultiSpecies();
// #else
//         real_t c2 = (gamma - 1.0) * (H - 0.5 * q2); // sound speed for H
// #endif
//         // NOTE:这里的gamma是混合物的gamma值，需要进行计算
//         real_t cl = sqrt(gamma * pl / qleft[ID]);
//         real_t cr = sqrt(gamma * pr / qright[ID]);
//         c = (cl + D * cr) * D1;
//         if (c2 < 0)
//         {
//             c = (cl + D * cr) * D1;
//             std::cout << "c2 < 0,please modify _c" << std::endl;
//         }
//         real_t b1 = (gamma - 1.0) / c2;
//         real_t b2 = 1.0 + b1 * q2 - b1 * H;
//         HydroSpecies zi, Yil, Yir, Yi, hil, hir, hi; // TODO:z
//         Yil[NUM_SPECIES - 1] = 1;
//         Yir[NUM_SPECIES - 1] = 1;
//         for (size_t yi = 0; yi < NUM_SPECIES - 1; yi++)
//         {
//             Yil[yi] = 1 - qleft[IDY1];
//             Yir[yi] = 1 - qright[IDY1];
//             Yil[NUM_SPECIES - 1] -= Yil[yi];
//             Yir[NUM_SPECIES - 1] -= Yir[yi];
//             Yi[yi] = (Yil[yi] + D * Yir[yi]) * D1;
//         }
//         Yi[NUM_SPECIES - 1] = (Yil[NUM_SPECIES - 1] + D * Yir[NUM_SPECIES - 1]) * D1;

//         real_t b3 = 0;
//     }
//     real_t _c1 = 1.0 / _c;
//     using HYDRO_NBVAR = HYDRO_2D_NBVAR;

//     // left eigen vectors
//     eigen_l[0 * HYDRO_NBVAR + ID] = 0.5 * (b2 + M * _c1 + b3);
//     eigen_l[0 * HYDRO_NBVAR + IU] = -0.5 * (b1 * _u + n[0] * _c1);
//     eigen_l[0 * HYDRO_NBVAR + IV] = -0.5 * (b1 * _v + n[1] * _c1);
//     eigen_l[0 * HYDRO_NBVAR + IP] = 0.5 * b1;

//     eigen_l[1 * HYDRO_NBVAR + ID] = (1.0 - b2 - b3); //-q2 + _H;
//     eigen_l[1 * HYDRO_NBVAR + IU] = _u * b1;
//     eigen_l[1 * HYDRO_NBVAR + IV] = _v * b1;
//     eigen_l[1 * HYDRO_NBVAR + IP] = -1.0 * b1;

//     eigen_l[2 * HYDRO_NBVAR + ID] = (M * n[1] - _v) / (n[0] + n[1] + n[2]);
//     eigen_l[2 * HYDRO_NBVAR + IU] = -1.0 * n[1];
//     eigen_l[2 * HYDRO_NBVAR + IV] = (1 - n[1] * n[1]) / (n[0] + n[1] + n[2]);
//     eigen_l[2 * HYDRO_NBVAR + IP] = 0.0;

//     eigen_l[3 * HYDRO_NBVAR + ID] = 0.5 * (b2 - M * _c1 + b3);
//     eigen_l[3 * HYDRO_NBVAR + IU] = 0.5 * (-b1 * _u + n[0] * _c1);
//     eigen_l[3 * HYDRO_NBVAR + IV] = 0.5 * (-b1 * _v + n[1] * _c1);
//     eigen_l[3 * HYDRO_NBVAR + IP] = 0.5 * b1;

//     // right eigen vectors
//     eigen_r[ID * HYDRO_NBVAR + 0] = 1.0;
//     eigen_r[IP * HYDRO_NBVAR + 0] = _H - _c * M;
//     eigen_r[IU * HYDRO_NBVAR + 0] = _u - _c * n[0];
//     eigen_r[IV * HYDRO_NBVAR + 0] = -v - _c * n[1];

//     eigen_r[ID * HYDRO_NBVAR + 1] = 1.0;
//     eigen_r[IP * HYDRO_NBVAR + 1] = _H - 1.0 / b1;
//     eigen_r[IU * HYDRO_NBVAR + 1] = _u;
//     eigen_r[IV * HYDRO_NBVAR + 1] = _v;

//     eigen_r[ID * HYDRO_NBVAR + 2] = 0.0;
//     eigen_r[IP * HYDRO_NBVAR + 2] = _v * n[0] - _u * n[1];
//     eigen_r[IU * HYDRO_NBVAR + 2] = -1.0 * n[1];
//     eigen_r[IV * HYDRO_NBVAR + 2] = n[0];

//     eigen_r[ID * HYDRO_NBVAR + 3] = 1.0;
//     eigen_r[IP * HYDRO_NBVAR + 3] = _H + _c * M;
//     eigen_r[IU * HYDRO_NBVAR + 3] = _u + _c * n[0];
//     eigen_r[IV * HYDRO_NBVAR + 3] = _v + c * n[1];

//     eigen_value[ID] = fabs(_u - _c);
//     eigen_value[IP] = fabs(_u + _c);
//     eigen_value[IU] = fabs(_u);
//     eigen_value[IV] = eigen_value[IU];

// #if Mixture
//     for (size_t i0 = HYDRO_NBVAR_EXCOP; i0 < HYDRO_NBVAR; i0++)
//     {
//         eigen_l[0 * HYDRO_NBVAR + i0] = -0.5 * b1 * z[i0 - HYDRO_NBVAR_EXCOP];
//     }
//     for (size_t i1 = HYDRO_NBVAR_EXCOP; i1 < HYDRO_NBVAR; i1++)
//     {
//         eigen_l[1 * HYDRO_NBVAR + i1] = b1 * z[i1 - HYDRO_NBVAR_EXCOP];
//     }
//     for (size_t i2 = HYDRO_NBVAR_EXCOP; i2 < HYDRO_NBVAR; i2++)
//     {
//         eigen_l[2 * HYDRO_NBVAR + i2] = 0.0;
//     }
//     for (size_t i3 = HYDRO_NBVAR_EXCOP; i3 < HYDRO_NBVAR; i3++)
//     {
//         eigen_l[3 * HYDRO_NBVAR + i3] = -0.5 * b1 * z[i3 - HYDRO_NBVAR_EXCOP];
//     }
//     for (size_t ir = 0; ir < NUM_SPECIES; ir++)
//     {
//         eigen_l[(ir + HYDRO_NBVAR_EXCOP) * HYDRO_NBVAR] = -1.0 * y[ir];
//         for (size_t ico1 = 1; ico1 < HYDRO_NBVAR_EXCOP; ico1++)
//         {
//             eigen_l[(ir + HYDRO_NBVAR_EXCOP) * HYDRO_NBVAR + ico1] = 0.0;
//         }
//         for (size_t ico2 = 0; ico2 < NUM_SPECIES; ico2++)
//         {
//             if (ir == ico2)
//                 eigen_l[(ir + HYDRO_NBVAR_EXCOP) * HYDRO_NBVAR + HYDRO_NBVAR_EXCOP + ico2] = 1.0;
//             else
//                 eigen_l[(ir + HYDRO_NBVAR_EXCOP) * HYDRO_NBVAR + HYDRO_NBVAR_EXCOP + ico2] = 0.0;
//         }
//     }

//     for (size_t j0 = HYDRO_NBVAR_EXCOP; j0 < HYDRO_NBVAR; j0++)
//     {
//         eigen_r[j0 * HYDRO_NBVAR] = y[j0 - HYDRO_NBVAR_EXCOP];
//     }
//     for (size_t j1 = HYDRO_NBVAR_EXCOP; j1 < HYDRO_NBVAR; j1++)
//     {
//         eigen_r[j1 * HYDRO_NBVAR + 1] = y[j1 - HYDRO_NBVAR_EXCOP];
//     }
//     for (size_t j2 = HYDRO_NBVAR_EXCOP; j2 < HYDRO_NBVAR; j2++)
//     {
//         eigen_r[j2 * HYDRO_NBVAR + 2] = 0.0;
//     }
//     for (size_t j2 = HYDRO_NBVAR_EXCOP; j2 < HYDRO_NBVAR; j2++)
//     {
//         eigen_r[j2 * HYDRO_NBVAR + 3] = y[j2 - HYDRO_NBVAR_EXCOP];
//     }
//     for (size_t jco = 0; jco < NUM_SPECIES; jco++)
//     {
//         eigen_r[ID * HYDRO_NBVAR + HYDRO_NBVAR_EXCOP + jco] = 0.0;
//         eigen_r[IU * HYDRO_NBVAR + HYDRO_NBVAR_EXCOP + jco] = 0.0;
//         eigen_r[IV * HYDRO_NBVAR + HYDRO_NBVAR_EXCOP + jco] = 0.0;
//         eigen_r[IP * HYDRO_NBVAR + HYDRO_NBVAR_EXCOP + jco] = z[jco - HYDRO_NBVAR_EXCOP];
//         for (size_t jro = 0; jro < NUM_SPECIES; jro++)
//         {
//             if (jro == jco)
//                 eigen_r[(jro + HYDRO_NBVAR_EXCOP) * HYDRO_NBVAR + HYDRO_NBVAR_EXCOP + jco] = 1.0;
//             else
//                 eigen_r[(jro + HYDRO_NBVAR_EXCOP) * HYDRO_NBVAR + HYDRO_NBVAR_EXCOP + jco] = 0.0;
//         }
//     }

//     for (size_t iev = HYDRO_NBVAR_EXCOP; iev < HYDRO_NBVAR; iev++)
//     {
//         eigen_l[iev] = eigen_value[IU];
//     }
// #endif
// }

//-------------------------------------------------------------------------------------------------
//					Non-dimensionalize
//-------------------------------------------------------------------------------------------------
// pressure
class DMS
{
public:
    real_t p, rho, t, l, u, v, w, miu;
    real_t _v, _rho, _length;
    DMS(/* args */);
    ~DMS();
    real_t non_dms_p(double p)
    {
        return p / _v / _v / _rho;
    }
    // density
    real_t non_dms_rho(double rho)
    {
        return rho / _rho;
    }
    // time
    real_t non_dms_time(double time)
    {
        return time * _v / _length;
    }
    // length
    real_t non_dms_length(double length)
    {
        return length / _length;
    }
    // velocity
    real_t non_dms_velocity(double velocity)
    {
        return velocity / _v;
    }
    // viscosity
    real_t non_dms_viscosity(double mu)
    {
        return mu / _v / _rho / _length;
    }
    //-------------------------------------------------------------------------------------------------
    //									Dimensionalize
    //-------------------------------------------------------------------------------------------------
    // pressure
    real_t dms_p(double p_non)
    {
        return p_non * _v * _v * _rho;
    }
    // density
    real_t dms_rho(double rho_non)
    {
        return rho_non * _rho;
    }
    // time
    real_t dms_time(double time_non)
    {
        return time_non / _v * _length;
    }
    // length
    real_t dms_length(double length_non)
    {
        return length_non * _length;
    }
    // velocity
    real_t dms_velocity(double velocity_non)
    {
        return velocity_non * _v;
    }
};
