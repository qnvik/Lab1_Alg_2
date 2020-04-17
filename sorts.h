#ifndef VICTORIATESTSORTS_SORTS_H
#define VICTORIATESTSORTS_SORTS_H

#include <array>
#include <cmath>
#include <numeric>
#include <random>
#include <cassert>

template <typename Iterator>
void introsort(Iterator begin, Iterator end);

namespace introsort_implementation {

    template<typename RandomIterator, typename Comparator>
    bool insertion_sort(RandomIterator begin, RandomIterator end, Comparator comp,
                        const int MAX_DIST = std::numeric_limits<int>::max()) {

        for (auto i = begin + 1; i < end; ++i){
            auto x = std::move(*i);
            decltype(i) j;
            int d = 0;

            for (j = i; j > begin && comp(x, *(j - 1)); j--){
                *j = std::move(*(j - 1));

                if (++d > MAX_DIST)
                    return false;
            }

            *j = x;
        }

        return true;
    };

#define sort_net2(x, y) if (comp(y, x)) std::swap(x, y)

    template<typename T, typename Comparator>
    void sort2(T &a, T &b, Comparator comp) {
        sort_net2(a, b);
    };

    template<typename T, typename Comparator>
    void sort3(T &a, T &b, T &c, Comparator comp) {
        sort_net2(b, c);
        sort_net2(a, c);
        sort_net2(a, b);
    };

    template<typename T, typename Comparator>
    void sort4(T &a, T &b, T &c, T &d, Comparator comp) {
        sort_net2(a, b);
        sort_net2(c, d);
        sort_net2(a, c);
        sort_net2(b, d);
        sort_net2(b, c);
    };

    template<typename T, typename Comparator>
    void sort5(T &a, T &b, T &c, T &d, T &e, Comparator comp) {
        sort_net2(a, b);
        sort_net2(d, e);
        sort_net2(c, e);
        sort_net2(c, d);
        sort_net2(b, e);
        sort_net2(a, d);
        sort_net2(a, c);
        sort_net2(b, d);
        sort_net2(b, c);
    };

    template<typename T, typename Comparator>
    void sort6(T &a, T &b, T &c, T &d, T &e, T &f, Comparator comp) {
        sort_net2(b, c);
        sort_net2(e, f);
        sort_net2(a, c);
        sort_net2(d, f);
        sort_net2(a, b);
        sort_net2(d, e);
        sort_net2(c, f);
        sort_net2(a, d);
        sort_net2(b, e);
        sort_net2(c, e);
        sort_net2(b, d);
        sort_net2(c, d);
    };

    constexpr int CUT_OFF = 25;
    constexpr int MAX_DIST = 15;

    template<typename RandomIterator, typename Comparator>
    typename std::iterator_traits<RandomIterator>::value_type
    choose_pivot(RandomIterator low, RandomIterator high, Comparator comp) {
        static std::random_device rd;
        static std::mt19937 mt(rd());
        static std::uniform_int_distribution<> distribution(0, std::numeric_limits<int>::max());

        const auto length = high - low + 1;

        typename std::iterator_traits<RandomIterator>::value_type x = *(low + distribution(mt) % length);
        typename std::iterator_traits<RandomIterator>::value_type y = *(low + distribution(mt) % length);
        typename std::iterator_traits<RandomIterator>::value_type z = *(low + distribution(mt) % length);

        sort3(x, y, z, comp);
        return y;
    }

    template<typename RandomIterator, typename Comparator>
    std::pair<RandomIterator, RandomIterator>
    three_way_partition(RandomIterator low, RandomIterator high, Comparator comp) {
        assert(high - low + 1 >= 7);

        RandomIterator lt = low;
        RandomIterator gt = high;
        auto pivot = choose_pivot(low, high, comp);

        auto i = low;
        while (i <= gt) {
            if (comp(*i, pivot))
                std::iter_swap(lt++, i++);
            else if (comp(pivot, *i))
                std::iter_swap(i, gt--);
            else
                i++;
        }

        return std::make_pair(lt, gt);
    }

    template<typename RandomIterator, typename Comparator>
    void impl_intro_sort(RandomIterator low, RandomIterator high, unsigned depth, Comparator comp) {
        while (low <= high) {
            const auto length = high - low + 1;

            /* Trivial case */
            if (length <= 1)
                return;

            /* For very small ranges (< 7 elements) we use a sorting network */
            switch (length) {
                case 2:
                    sort2(*low, *(low + 1), comp);
                    return;
                case 3:
                    sort3(*low, *(low + 1), *(low + 2), comp);
                    return;
                case 4:
                    sort4(*low, *(low + 1), *(low + 2), *(low + 3), comp);
                    return;
                case 5:
                    sort5(*low, *(low + 1), *(low + 2), *(low + 3), *(low + 4), comp);
                    return;
                case 6:
                    sort6(*low, *(low + 1), *(low + 2), *(low + 3), *(low + 4), *(low + 5), comp);
                    return;
                default:
                    break;
            }

            /* We have less than CUT_OFF elements so we must use insertion sort */
            if (length <= CUT_OFF) {
                insertion_sort(low, high + 1, comp);
                return;
            }

            /* Quicksort is not behaving properly so we switch to heapsort */
            if (depth == 0) {
                std::make_heap(low, high + 1, comp);
                std::sort_heap(low, high + 1, comp);
                return;
            }

            // try gambling ? (optimization)
            if (insertion_sort(low, high + 1, comp, MAX_DIST))
                return;

            /*
             * We use 3-way partition in order to get 3 ranges:
             *  a[low..pit.first) ++ a[pit.first..pit.second] ++ a(pit.second..high]
             *      < pivot                  == pivot                  > pivot
             */
            auto pit = three_way_partition(low, high, comp);

            /*
             * If len(a[low..pit.first)) < len(a(pit.second..high]) then
             * transform sort(a(pit.second..high]) into tail recursion
             * Otherwise transform sort(a[low..pit.first)) into tail recursion
             */
            if (pit.first - low < high - pit.second){
                impl_intro_sort(low, pit.first - 1, depth - 1, comp);
                low = pit.second + 1;
            }
            else{
                impl_intro_sort(pit.second + 1, high, depth - 1, comp);
                high = pit.first - 1;
            }

            depth -= 1;
        }
    };

    template<typename RandomIterator>
    unsigned find_maximal_depth(RandomIterator low, RandomIterator high) {
        return (sizeof(int) * __CHAR_BIT__  - 1 - __builtin_clz(high - low)) * 2;
    }

    template <typename RandomIterator>
    void introsort(RandomIterator begin, RandomIterator end, std::random_access_iterator_tag){
        if (begin != end)
            impl_intro_sort(begin, --end, find_maximal_depth(begin, end),
                            std::less<typename std::iterator_traits<RandomIterator>::value_type>());
    };

    template <typename BidirectionalIterator>
    void introsort(BidirectionalIterator begin, BidirectionalIterator end, std::bidirectional_iterator_tag){
        std::vector<typename std::iterator_traits<BidirectionalIterator>::value_type> container(begin, end);
        introsort(container.begin(), container.end(), std::random_access_iterator_tag());
        std::move(container.begin(), container.end(), begin);
    };

    template <typename ForwardIterator>
    void introsort(ForwardIterator begin, ForwardIterator end, std::forward_iterator_tag){
        std::vector<typename std::iterator_traits<ForwardIterator>::value_type> container(begin, end);
        introsort(container.begin(), container.end(), std::random_access_iterator_tag());
        std::move(container.begin(), container.end(), begin);
    };

    template <typename Iterator>
    void introsort(Iterator begin, Iterator end){
        introsort_implementation::introsort(begin, end, typename std::iterator_traits<Iterator>::iterator_category());
    };
}

template <typename Iterator>
void introsort(Iterator begin, Iterator end){
    introsort_implementation::introsort(begin, end);
};

namespace details
{
    const int kNoDigit = -1;
    int ExtractDigit(int i, int pos) {
        const int digitsCount = log10(i) + 1;
        if (pos > digitsCount) return kNoDigit;
        return (int)(i / pow(10, digitsCount - pos)) % 10;
    }

    // For pos equals to 2 and {10 20 1} -> {0, 1, 3, 3, ...},
    // 1 ends with empty digit in the second digit and 2 ends with 0
    template <class It>
    auto CountingSort(It begin, It end, int pos) {
        std::array<int, 12> bins;
        std::fill(bins.begin(), bins.end(), 0);
        for (auto it = begin; it < end; ++it) {
            const int digit = ExtractDigit(*it, pos);
            ++bins[digit + 1];
        }
        std::partial_sum(bins.cbegin(), bins.cend(), bins.begin());
        std::move(bins.cbegin(), bins.cend() - 1, bins.begin() + 1);
        return bins;
    }

    template <class It>
    void MsdRadixInternal(It begin, It end, int pos) {
        const auto bins = CountingSort(begin, end, pos);
        // We finish when i is 1, because the last part ends up sorted anyway
        for (int i = 10; i > 0; --i) {
            const int digit = i - 1;
            const auto local_begin = begin + bins[i];
            const auto local_end = begin + bins[i + 1];
            if (local_begin == begin) break;
            if (std::distance(local_begin, local_end) > 0) {
                auto crsrForward = begin;
                auto crsrBackward = local_end - 1;
                while (crsrForward < crsrBackward) {
                    while (ExtractDigit(*crsrBackward, pos) == digit) --crsrBackward;
                    while (ExtractDigit(*crsrForward, pos) != digit) ++crsrForward;
                    if (crsrForward < local_begin) {
                        std::swap(*crsrBackward, *crsrForward);
                    }
                    ++crsrForward;
                }
            }
        }
        // Start from 1 as we don't want to sort numbers wich are out of digits in pos already
        for (int i = 1; i < 11; ++i) {
            if (bins[i + 1] - bins[i] > 1)
                MsdRadixInternal(begin + bins[i], begin + bins[i + 1], pos + 1);
        }
    }
}

template <class It>
void MsdRadix(It begin, It end) {
    details::MsdRadixInternal(begin, end, 1);
}


    static void shell_sort( auto first, auto last)
    {
        auto size = last - first;
        auto h = 1;
        while( h < size / 3 ) h = 3 * h + 1;

        while( h >= 1 )
        {
            for( auto i = first + h; i != last; i++ )
            {
                for( auto j = i; (j - first) >= h && (*j < *(j-h)); j -= h )
                    std::swap( *j, *(j-h) );
            }

            h /= 3;
        }
    }


#endif //VICTORIATESTSORTS_SORTS_H
