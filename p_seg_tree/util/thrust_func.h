#ifndef THURST_FUNC_H_
#define THURST_FUNC_H_

void thrustDeviceStableSortByKey(REAL *arr, int *arrIdx, int size);
void thrustDeviceStableSortByKey(int *arr, int *arrIdx, int size);
void thrustDeviceStableSortByKey(LONG *arr, int *arrIdx, int size);
void thrustDeviceExclusiveScan(int *arr, int size);
void thrustDeviceExclusiveScan(long *arr, int size);
void thrustDeviceSort(REAL *arr, int size);

void thrustHostSort(REAL *arr, int size);
void thrustHostStableSortByKey(REAL *arr, int *arrIdx, int size);

#endif /* THURST_FUNC_H_ */