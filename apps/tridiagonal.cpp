#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <iostream>


struct data_type {
    using View2D = Kokkos::View<double**, Kokkos::LayoutRight>;
    using View1D = Kokkos::View<double*>;
    using ViewInt1D = Kokkos::View<int*>;

    View2D a;  // sub-diagonal
    View2D b;  // diagonal
    View2D c;  // super-diagonal
    View2D r;
    View2D x;  // solution
    View2D gam;
    ViewInt1D jtop;

    // Constructor
    data_type(int ncol, int nz, const std::string &name_prefix)
    : a(name_prefix + "_a", ncol, nz),
      b(name_prefix + "_b", ncol, nz),
      c(name_prefix + "_c", ncol, nz),
      r(name_prefix + "_r", ncol, nz),
      x(name_prefix + "_x", ncol, nz),
      gam(name_prefix + "_gam", ncol, nz),
      jtop(name_prefix + "_jtop", ncol) {}

    void initialize(int ncol, int nz) {
        auto a_h = Kokkos::create_mirror_view(a);
        auto b_h = Kokkos::create_mirror_view(b);
        auto c_h = Kokkos::create_mirror_view(c);
        auto r_h = Kokkos::create_mirror_view(r);
        auto x_h = Kokkos::create_mirror_view(x);
        auto gam_h = Kokkos::create_mirror_view(gam);
        auto j_h = Kokkos::create_mirror_view(jtop);
        //     1 1 0      6
        // A = 2 7 8  r = 9   
        //     0 3 5      6
        for (int ci = 0; ci < ncol; ++ci) {
          for (int j = 0; j < nz; ++j) {
            a_h(ci,j) = 0.0;
            b_h(ci,j) = 0.0;
            c_h(ci,j) = 0.0;
            r_h(ci,j) = 0.0;
            x_h(ci,j) = 0.0;
            gam_h(ci,j) = 0.0;
          }
          j_h(ci) = 0;
        }

        Kokkos::deep_copy(a, a_h);
        Kokkos::deep_copy(b, b_h);
        Kokkos::deep_copy(c, c_h);
        Kokkos::deep_copy(r, r_h);
        Kokkos::deep_copy(x, x_h);
        Kokkos::deep_copy(gam, gam_h);
        Kokkos::deep_copy(jtop,j_h);
    }

    KOKKOS_INLINE_FUNCTION
    void fill(int i, int j, double aij, double bij, double cij, double rij) const {
        a(i,j) = aij;
        b(i,j) = bij;
        c(i,j) = cij;
        r(i,j) = rij;
    }

    void print_column(int ci, int nz) {
        auto x_h = Kokkos::create_mirror_view(x);
        Kokkos::deep_copy(x_h, x);
        std::cout << "Column " << ci << ": ";
        for (int j = 0; j < nz; ++j)
          std::cout << x_h(ci,j) << " ";
        std::cout << std::endl;
      }

};

struct TridiagLUFunctor {
    // Thomas Algorithm for solving Tridiagonal system equations
    // Ax = r 
    // A = LU, L the Lower-triangular matrix, and U the Upper-triangular matrix
    // LUx = r -> Ly = r, Ux = y
    // 1. Forward elimination: solve Ly = r
    // 2. Back substitution:   solve Ux = y
    using View1D = Kokkos::View<int*>;
    using View1DBool = Kokkos::View<bool*>;

    int nz;
    data_type data;   // store a, b, c, x, jtop

    // Constructor
    TridiagLUFunctor(int nz_,
                   data_type data_)
    : nz(nz_), data(data_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int ci) const {

        const int jt = data.jtop(ci);
        double bet = data.b(ci,jt);
        
        // Forward LU factorization
        for (int j = jt; j < nz; ++j) {
            if (j == jt) {
                data.x(ci,j) = data.r(ci,j) / bet;
            } else {
                data.gam(ci,j) = data.c(ci,j-1) / bet;
                bet = data.b(ci,j) - data.a(ci,j) * data.gam(ci,j);
                data.x(ci,j)   = (data.r(ci,j) - data.a(ci,j)*data.x(ci,j-1)) / bet;
            }
        }

        for (int j = nz - 2; j >= jt; --j) {
            data.x(ci,j) = data.x(ci,j) - data.gam(ci,j+1) * data.x(ci,j+1);
        }
    }
};

int main(int argc, char* argv[]) {

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int mpi_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    Kokkos::initialize(argc, argv);
    {
        // Safest: host-only print guard
        if (mpi_rank == 0) {
        // Ensure only master thread prints
        #pragma omp master
        {
            std::cout << "=== Kokkos Configuration (rank 0 only) ===" << std::endl;
            Kokkos::print_configuration(std::cout);
        }

        using ExecSpace = Kokkos::DefaultExecutionSpace;
        using Layout    = Kokkos::LayoutRight;

        const int ncol = 1;   // total columns
        const int nz   = 3;   // vertical levels

        // ------------------
        // Allocate views
        // ------------------
        data_type data(ncol, nz, "data");
        data.initialize(ncol, nz);
        // Fill A 
        Kokkos::parallel_for(
            "Fill_A",
            Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{ncol,nz}),
            KOKKOS_LAMBDA(int i, int j) {
                if (j == 0) {
                    data.fill(i,j,0.0,1.0,1.0,6.0);
                } else if (j == 1) {
                    data.fill(i,j,2.0,7.0,8.0,9.0);
                } else if (j == 2) {
                    data.fill(i,j,3.0,5.0,0.0,6.0);
                }
            }
        );
        Kokkos::fence(); 
        // Launch
        Kokkos::parallel_for(
            "TridiagLU",
            Kokkos::RangePolicy<>(0, ncol), // Kokkos will use the default execution space (Kokkos::DefaultExecutionSpace).
            TridiagLUFunctor(nz, data)
        );
        Kokkos::fence();
        // Print results for column 0
        data.print_column(0, nz);
        }
    }

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
