#pragma once
#include "shared/kokkos_shared.h"
//#include "shared/HydroParams.h"
#include "shared/HydroState.h"
#include <string.h>
#include <iostream>
#include <fstream>

void readOmega(DataArray2dHost &Omega_tableHost, DataArrayOneHost &T_starHost, DataArrayOneHost &delta_starHost)
{
    std::ifstream fin("./collision_integral.dat");
    for (int n = 0; n < 8; n++)
    {
        fin >> delta_starHost(n); // reduced dipole moment;
        // std::cout << delta_starHost(n) << std::endl;
    }
    for (int i = 0; i < 37; i++)
    {
        fin >> T_starHost(i); // reduced temperature;
        // std::cout << T_starHost(i) << std::endl;
        for (int j = 0; j < 8; j++)
        {
            fin >> Omega_tableHost(1, i, j); // collision integral for binary diffusion coefficient;
            // std::cout << Omega_tableHost(1, i, j) << std::endl;
        }
    }
    for (int p = 0; p < 37; p++)
        for (int q = 0; q < 8; q++)
        {
            fin >> Omega_tableHost(0, p, q); // collision integral for viscosity and thermal conductivity;
            // std::cout << Omega_tableHost(0, p, q) << std::endl;
        }
    fin.close();
    // std::cout << Omega_tableHost(1, 2, 4) << std::endl;
};
//};

void readSpecie(std::string *species_name, DataArrayScalarHost &CharaHost) //: Species_name(species_name), CharaHost(6)(wi), CharaHost(7)(Miu) // Wi(wi), miu(Miu)
{
    for (size_t i = 0; i < NUM_SPECIES; i++)
    {
        char Key_word[128];
        std::fstream fin("./transport_data.dat");
        // check the name of the species "n"
        std::string my_string = species_name[i];
        char *species_name_n = new char[my_string.size()];
        strcpy(species_name_n, my_string.c_str());
        // reset file point location
        fin.seekg(0);
        while (!fin.eof())
        {
            fin >> Key_word;
            if (!strcmp(Key_word, "*END"))
                break;
            if (!strcmp(Key_word, species_name_n))
            {
                fin >> CharaHost(i, 0); // geo
                fin >> CharaHost(i, 1); // epsilon_kB;
                fin >> CharaHost(i, 2); // d;
                fin >> CharaHost(i, 3); // mue;
                fin >> CharaHost(i, 4); // alpha;
                fin >> CharaHost(1, 5); // Zrot_298;
                break;
            }
        }
        fin.close();
        CharaHost(i, 7) = 0; // Miu[i]; // miu}
                             //  CharaHost(i, 6) = wi[i];  // Wi,在readHiArgus内读取
                             //  double epsilon_kB;//epsilon: Lennard-Jones potential well depth;unit:K//势井深度
                             //  double mue;//dipole moment,unit:Debye(m);//偶极距
                             //  double d;//Lennard-Jones collision diameter ,unit: angstroms,10e-10m//碰撞直径in 4-3
                             //  int geo;//0:monoatom,1:nonpolar(linear) molecule,2:polar molecule//极性
                             //  double alpha;//polarizability;unit:cubic angstrom//极化率
                             //  double Zrot_298;//rotational relaxation collision Zrot at 298K;
                             //  NOTE:气体性质转换为国际制单位,麻烦
        // CharaHost(i, 0) *= 1.0;   // geo
        // CharaHost(i, 1) *= 1.0;   // epsilon_kB;
        // CharaHost(i, 2) *= 1.0;   // d;
        // CharaHost(i, 3) *= 1.0;   // mue;
        // CharaHost(i, 4) *= 1.0;   // alpha;//FAQ?(0~1)?
        // CharaHost(1, 5) *= 1.0;   // Zrot_298;//FAQ?
    }
};

void readHiArgus(std::string *species_name, DataArray2dHost &a, DataArray2dHost &b, DataArrayScalarHost &CharaHost)
{
#define Thermo 1 // 1 for NASA and 0 for JANAF
    char Key_word[128];
#if Thermo
    std::fstream fin("./thermal_dynamics.dat");
#else
    std::fstream fin("./thermal_dynamics_janaf.dat");
#endif
    for (int n = 0; n < NUM_SPECIES; n++)
    {
        // check the name of the species "n"
        std::string my_string = "*" + species_name[n];
        char *species_name_n = new char[my_string.size() + 1];
        std::strcpy(species_name_n, my_string.c_str());
        // reset file point location
        fin.seekg(0);
        while (!fin.eof())
        {
            fin >> Key_word;
            if (!std::strcmp(Key_word, "*END"))
                break;
            if (!std::strcmp(Key_word, species_name_n))
            {
#if Thermo
                // low temperature parameters, 200K<T<1000K
                for (int m = 0; m < 7; m++)
                    fin >> a(n, m, 0); // a1-a7
                for (int m = 0; m < 2; m++)
                    fin >> b(n, m, 0); // b1,b2
                // high temperature parameters, 1000K<T<6000K
                for (int m = 0; m < 7; m++)
                    fin >> a(n, m, 1); // a1-a7
                for (int m = 0; m < 2; m++)
                    fin >> b(n, m, 1); // b1,b2
                // high temperature parameters, 6000K<T<15000K
                for (int m = 0; m < 7; m++)
                    fin >> a(n, m, 2); // a1-a7
                for (int m = 0; m < 2; m++)
                    fin >> b(n, m, 2); // b1,b2
#else
                // high temperature parameters
                for (int m = 0; m < 7; m++)
                    fin >> species[n].a[m][0]; // a1-a7
                // low temperature parameters
                for (int m = 0; m < 7; m++)
                    fin >> species[n].a[m][1]; // a1-a7
#endif
                // molar mass, unit: g/mol
                fin >> CharaHost(n, 6); // species[n].Wi;
                CharaHost(n, 6) *= 1e-3;
                break;
            }
        }
    }
    fin.close();
}