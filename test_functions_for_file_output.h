#ifndef TESTSORTS_TEST_FUNCTIONS_FOR_FILE_OUTPUT_H
#define TESTSORTS_TEST_FUNCTIONS_FOR_FILE_OUTPUT_H

#include <chrono>
#include <random>
#include <fstream>
#include <algorithm>
#include <vector>
#include "sorts.h"

const short STD_SORT = 0;
const short STD_STABLE_SORT = 1;
const short STD_PARTIAL_SORT = 2;
const short MSD_RADIX_SORT = 3;
const short INTRO_SORT = 4;
const short SHELL_SORT = 5;

namespace fstreamVersion {
    void uniform_0to2_in_pow_31_std_sort(const std::uint32_t N, const short SORT_TYPE, std::ostream& fout) {
        std::uint32_t array[N];
        std::vector<std::uint32_t> vec(N);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<std::uint32_t> dis(0u, (1u << 31)); // 2^31
        long testDataArray[5];
        long testDataVector[5];
        auto now = std::chrono::high_resolution_clock::now;
        fout << "vec[i] = array[i] ~ U(0, 2^31)" << std::endl;
        for (auto i = 0u; i < 5; ++i) {
            for (auto j = 0u; j<N; ++j) {
                vec[j] = array[j] = dis(gen);
            }

            auto startOfSort = now();
            if (SORT_TYPE == STD_SORT)
                std::sort(array, array + N);
            else if (SORT_TYPE == STD_STABLE_SORT)
                std::stable_sort(array, array + N);
            else if (SORT_TYPE == STD_PARTIAL_SORT)
                std::partial_sort(array, array + (N / 2), array + N);
            else if(SORT_TYPE == MSD_RADIX_SORT) {
                MsdRadix(array, array + N);
            }
            else if(SORT_TYPE == INTRO_SORT) {
                introsort(array, array + N);
            }
            else if(SORT_TYPE == SHELL_SORT) {
                shell_sort(array, array + N);
            }
            auto endOfSort = now();
            testDataArray[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSort - startOfSort).count();

            startOfSort = now();
            if (SORT_TYPE == STD_SORT) {
                std::sort(vec.begin(), vec.end());
            }
            else if (SORT_TYPE == STD_STABLE_SORT) {
                std::stable_sort(vec.begin(), vec.end());
            }
            else if (SORT_TYPE == STD_PARTIAL_SORT) {
                std::partial_sort(vec.begin(), vec.begin()+vec.size()/2,vec.end());
            }
            else if(SORT_TYPE == MSD_RADIX_SORT) {
                MsdRadix(vec.begin(), vec.end());
            }
            else if(SORT_TYPE == INTRO_SORT) {
                introsort(vec.begin(), vec.end());
            }
            else if(SORT_TYPE == SHELL_SORT) {
                shell_sort(vec.begin(), vec.end());
            }

            endOfSort = now();
            testDataVector[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSort - startOfSort).count();

        }
        fout << "Average (array): " << std::accumulate(testDataArray, testDataArray + 5, 0l) / 5. << "ms"  << std::endl;
        fout << "Average (vector): " << std::accumulate(testDataVector, testDataVector + 5, 0l) / 5. << "ms" << std::endl << std::endl;
    }

    void uniform_0to2_in_pow_15_std_sort(const std::uint32_t N, const short SORT_TYPE, std::ostream& fout) {
        std::uint32_t array[N];
        std::vector<std::uint32_t> vec(N);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<std::uint32_t> dis(0u, (1u << 15)); // 2^15
        long testDataArray[5];
        long testDataVector[5];
        auto now = std::chrono::high_resolution_clock::now;
        fout << "vec[i] = array[i] ~ U(0, 2^31)" << std::endl;
        for (auto i = 0u; i < 5; ++i) {
            for (auto j = 0u; j<N; ++j) {
                vec[j] = array[j] = dis(gen);
            }

            auto startOfSort = now();
            if (SORT_TYPE == STD_SORT)
                std::sort(array, array + N);
            else if (SORT_TYPE == STD_STABLE_SORT)
                std::stable_sort(array, array + N);
            else if (SORT_TYPE == STD_PARTIAL_SORT)
                std::partial_sort(array, array + (N / 2), array + N);
            else if(SORT_TYPE == MSD_RADIX_SORT) {
                MsdRadix(array, array + N);
            }
            else if(SORT_TYPE == INTRO_SORT) {
                introsort(array, array + N);
            }
            else if(SORT_TYPE == SHELL_SORT) {
                shell_sort(array, array + N);
            }
            auto endOfSort = now();
            testDataArray[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSort - startOfSort).count();

            startOfSort = now();
            if (SORT_TYPE == STD_SORT) {
                std::sort(vec.begin(), vec.end());
            }
            else if (SORT_TYPE == STD_STABLE_SORT) {
                std::stable_sort(vec.begin(), vec.end());
            }
            else if (SORT_TYPE == STD_PARTIAL_SORT) {
                std::partial_sort(vec.begin(), vec.begin()+vec.size()/2,vec.end());
            }
            else if(SORT_TYPE == MSD_RADIX_SORT) {
                MsdRadix(vec.begin(), vec.end());
            }
            else if(SORT_TYPE == INTRO_SORT) {
                introsort(vec.begin(), vec.end());
            }
            else if(SORT_TYPE == SHELL_SORT) {
                shell_sort(vec.begin(), vec.end());
            }
            endOfSort = now();
            testDataVector[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSort - startOfSort).count();

        }
        fout << "Average (array): " << std::accumulate(testDataArray, testDataArray + 5, 0l) / 5. << "ms"  << std::endl;
        fout << "Average (vector): " << std::accumulate(testDataVector, testDataVector + 5, 0l) / 5. << "ms" << std::endl << std::endl;
    }

    void uniform_0to2_N_1_std_sort(const std::uint32_t N, const short SORT_TYPE, std::ostream& fout) {
        std::uint32_t array[N];
        std::vector<std::uint32_t> vec(N);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<std::uint32_t> dis(0u, N - 1);
        long testDataArray[5];
        long testDataVector[5];
        auto now = std::chrono::high_resolution_clock::now;
        fout << "vec[i] = array[i] ~ U(0, " << N - 1 << ")" << std::endl;
        for (auto i = 0u; i < 5; ++i) {
            for (auto j = 0u; j<N; ++j) {
                vec[j] = array[j] = dis(gen);
            }

            auto startOfSort = now();
            if (SORT_TYPE == STD_SORT)
                std::sort(array, array + N);
            else if (SORT_TYPE == STD_STABLE_SORT)
                std::stable_sort(array, array + N);
            else if (SORT_TYPE == STD_PARTIAL_SORT)
                std::partial_sort(array, array + (N / 2), array + N);
            else if(SORT_TYPE == MSD_RADIX_SORT) {
                MsdRadix(array, array + N);
            }
            else if(SORT_TYPE == INTRO_SORT) {
                introsort(array, array + N);
            }
            else if(SORT_TYPE == SHELL_SORT) {
                shell_sort(array, array + N);
            }
            auto endOfSort = now();
            testDataArray[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSort - startOfSort).count();

            startOfSort = now();
            if (SORT_TYPE == STD_SORT) {
                std::sort(vec.begin(), vec.end());
            }
            else if (SORT_TYPE == STD_STABLE_SORT) {
                std::stable_sort(vec.begin(), vec.end());
            }
            else if (SORT_TYPE == STD_PARTIAL_SORT) {
                std::partial_sort(vec.begin(), vec.begin()+vec.size()/2,vec.end());
            }
            else if(SORT_TYPE == MSD_RADIX_SORT) {
                MsdRadix(vec.begin(), vec.end());
            }
            else if(SORT_TYPE == INTRO_SORT) {
                introsort(vec.begin(), vec.end());
            }
            else if(SORT_TYPE == SHELL_SORT) {
                shell_sort(vec.begin(), vec.end());
            }
            endOfSort = now();
            testDataVector[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSort - startOfSort).count();

        }
        fout << "Average (array): " << std::accumulate(testDataArray, testDataArray + 5, 0l) / 5. << "ms"  << std::endl;
        fout << "Average (vector): " << std::accumulate(testDataVector, testDataVector + 5, 0l) / 5. << "ms" << std::endl << std::endl;
    }

    void normal_0to2_in_pow_31_std_sort(const std::uint32_t N, const short SORT_TYPE, std::ostream& fout) {
        std::uint32_t array[N];
        std::vector<std::uint32_t> vec(N);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dis{(1u << 31) / 2, (1u << 31)}; // 2^31
        long testDataArray[5];
        long testDataVector[5];
        auto now = std::chrono::high_resolution_clock::now;
        fout << "vec[i] = array[i] ~ N(2^31/2, 2^31)" << std::endl;
        for (auto i = 0u; i < 5; ++i) {
            for (auto j = 0u; j<N; ++j) {
                vec[j] = array[j] = dis(gen);
            }
            auto startOfSort = now();
            if (SORT_TYPE == STD_SORT)
                std::sort(array, array + N);
            else if (SORT_TYPE == STD_STABLE_SORT)
                std::stable_sort(array, array + N);
            else if (SORT_TYPE == STD_PARTIAL_SORT)
                std::partial_sort(array, array + (N / 2), array + N);
            else if(SORT_TYPE == MSD_RADIX_SORT) {
                MsdRadix(array, array + N);
            }
            else if(SORT_TYPE == INTRO_SORT) {
                introsort(array, array + N);
            }
            else if(SORT_TYPE == SHELL_SORT) {
                shell_sort(array, array + N);
            }
            auto endOfSort = now();
            testDataArray[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSort - startOfSort).count();

            startOfSort = now();
            if (SORT_TYPE == STD_SORT) {
                std::sort(vec.begin(), vec.end());
            }
            else if (SORT_TYPE == STD_STABLE_SORT) {
                std::stable_sort(vec.begin(), vec.end());
            }
            else if (SORT_TYPE == STD_PARTIAL_SORT) {
                std::partial_sort(vec.begin(), vec.begin()+vec.size()/2,vec.end());
            }
            else if(SORT_TYPE == MSD_RADIX_SORT) {
                MsdRadix(vec.begin(), vec.end());
            }
            else if(SORT_TYPE == INTRO_SORT) {
                introsort(vec.begin(), vec.end());
            }
            else if(SORT_TYPE == SHELL_SORT) {
                shell_sort(vec.begin(), vec.end());
            }
            endOfSort = now();
            testDataVector[i] = std::chrono::duration_cast<std::chrono::milliseconds>(endOfSort - startOfSort).count();

        }
        fout << "Average (array): " << std::accumulate(testDataArray, testDataArray + 5, 0l) / 5. << "ms"  << std::endl;
        fout << "Average (vector): " << std::accumulate(testDataVector, testDataVector + 5, 0l) / 5. << "ms" << std::endl << std::endl;
    }
};

#endif //TESTSORTS_TEST_FUNCTIONS_FOR_FILE_OUTPUT_H
