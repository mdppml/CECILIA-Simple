//
// Created by Mete Akgun on 18.07.20.
//

#ifndef PML_LLIB_H
#define PML_LLIB_H

#include <stdlib.h>
#include <stdint.h>

#define PARTCOUNT 100
#define MAXSAMPLESIZE 50000

void swap(uint64_t* a, uint64_t* b)
{
    uint64_t t = *a;
    *a = *b;
    *b = t;
}

int partition (uint64_t real_arr[], int low, int high)
{
    uint64_t pivot = real_arr[high]; // pivot
    int i = (low - 1); // Index of smaller element

    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (real_arr[j] > pivot)
        {
            i++; // increment index of smaller element
            swap(&real_arr[i], &real_arr[j]);
        }
    }
    swap(&real_arr[i + 1], &real_arr[high]);
    return (i + 1);
}
void quickSort(uint64_t real_arr[], int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
         at right place */
        int pi = partition(real_arr, low, high);

        // Separately sort elements before
        // partition and after partition
        quickSort(real_arr, low, pi - 1);
        quickSort(real_arr, pi + 1, high);
    }
}

void sort_values(uint64_t real_arr[], int n){
    quickSort(real_arr, 0, n);
}


#endif //PML_LLIB_H
