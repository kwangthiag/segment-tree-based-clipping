#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>


#include "segment_tree/data_types.h"
#include "segment_tree/edge.h"


// Define a custom comparison functor to compare by the 'start' attribute
struct compareByStart {
    __host__ __device__
    bool operator()(const Edge& a, const Edge& b) const {
        return a.start < b.start;
    }
};

/*device methods*/
void thrustDeviceStableSortByKey(REAL *arr, int *arrIdx, int size){
    thrust::stable_sort_by_key(thrust::device, arr, arr+size, arrIdx);
}
void thrustDeviceStableSortByKey(int *arr, int *arrIdx, int size){
    thrust::stable_sort_by_key(thrust::device, arr, arr+size, arrIdx);
}
void thrustDeviceStableSortByKey(LONG *arr, int *arrIdx, int size){
    thrust::stable_sort_by_key(thrust::device, arr, arr+size, arrIdx);
}
void thrustDeviceExclusiveScan(int *arr, int size){
    thrust::exclusive_scan(thrust::device, arr, arr + size, arr);
}
void thrustDeviceExclusiveScan(long *arr, int size){
    thrust::exclusive_scan(thrust::device, arr, arr + size, arr);
}
// void thrustDeviceSort(REAL *arr, int size){
//     thrust::sort(thrust::device, arr, arr+size, compareByStart());
// }

/*Host methods*/

void thrustHostStableSortByKey(REAL *arr, int *arrIdx, int size){
    thrust::stable_sort_by_key(thrust::host, arr, arr+size, arrIdx);
}
// void thrustHostSort(REAL *arr, int size){
//     thrust::sort(thrust::host, arr, arr+size, compareByStart());
// }
