#include <fstream>
#include "test_functions_for_file_output.h"


int main() {
    std::ofstream fout("file.dat");
    std::uint32_t Ni[] = {30000, 100000, 300000, 1000000};
    for(short type = 0; type<6; ++type) {
        if (type == STD_SORT) {
            fout << "std::sort\n";
        } else if (type == STD_STABLE_SORT) {
            fout << "std::stable_sort\n";
        } else if (type == STD_PARTIAL_SORT) {
            fout << "std::partial_sort\n";
        } else if (type == MSD_RADIX_SORT) {
            fout << "msd_radix_sort\n";
        } else if (type == INTRO_SORT) {
            fout << "intro_sort\n";
        } else if (type == SHELL_SORT) {
            fout << "shell_sort\n";
        }
        for (auto &N : Ni) {
            fout << "N = " << N << std::endl;
            fstreamVersion::uniform_0to2_in_pow_31_std_sort(N, type, fout);
            fstreamVersion::uniform_0to2_in_pow_15_std_sort(N, type, fout);
            fstreamVersion::uniform_0to2_N_1_std_sort(N, type, fout);
            fstreamVersion::normal_0to2_in_pow_31_std_sort(N, type, fout);
        }
    }

    return 0;
}
