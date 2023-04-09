 cmake -DUSE_MPI=ON -DCMAKE_INSTALL_PREFIX=./kokkos/ -DCMAKE_CXX_COMPILER=/public/software/compiler/rocm/dtk-22.04.2/bin/hipcc -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_VEGA906=ON ..
