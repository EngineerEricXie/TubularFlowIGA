// Wrap existing tests without changing their equations or weakening the caller gate.
#include <mpi.h>
#include <iostream>
#define main fsi_wrapped_test_main
#include IGA_WRAPPED_TEST
#undef main
int main(int argc, char** argv)
{
    int provided = MPI_THREAD_SINGLE;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS)
        return 2;
    if (provided < MPI_THREAD_FUNNELED) {
        std::cerr << "MPI_THREAD_FUNNELED unavailable\n";
        MPI_Finalize(); return 2;
    }
    const int status = fsi_wrapped_test_main(argc, argv);
    MPI_Finalize(); return status;
}
