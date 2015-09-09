/* quicksort */

#include <stdio.h>
#include <stdlib.h>
#include <mags3d.h>

//typedef int T;          /* type of item to be sorted */
//typedef int int;   /* type of subscript */

#define compGT(a,b) (a > b)

void insertSort(TYPE *a, int lb, int ub) {
    TYPE t;
    int i, j;

   /**************************
    *  sort array a[lb..ub]  *
    **************************/
    for (i = lb + 1; i <= ub; i++) {
        t = a[i];

        /* Shift elements down until */
        /* insertion point found.    */
        for (j = i-1; j >= lb && compGT(a[j], t); j--)
            a[j+1] = a[j];

        /* insert */
        a[j+1] = t;
    }
}

int partition(TYPE *a, int lb, int ub) {
    TYPE t, pivot;
    int i, j, p;

   /*******************************
    *  partition array a[lb..ub]  *
    *******************************/

    /* select pivot and exchange with 1st element */
    p = lb + ((ub - lb)>>1);
    pivot = a[p];
    a[p] = a[lb];

    /* sort lb+1..ub based on pivot */
    i = lb+1;
    j = ub;
    while (1) {
        while (i < j && compGT(pivot, a[i])) i++;
        while (j >= i && compGT(a[j], pivot)) j--;
        if (i >= j) break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
        j--; i++;
    }

    /* pivot belongs in a[j] */
    a[lb] = a[j];
    a[j] = pivot;

    return j;
}

void quickSort(TYPE *a, int lb, int ub) {
    int m;

   /**************************
    *  sort array a[lb..ub]  *
    **************************/

    while (lb < ub) {

        /* quickly sort short lists */
        if (ub - lb <= 12) {
            insertSort(a, lb, ub);
            return;
        }

        /* partition into two segments */
        m = partition (a, lb, ub);

        /* sort the smallest partition    */
        /* to minimize stack requirements */
        if (m - lb <= ub - m) {
            quickSort(a, lb, m - 1);
            lb = m + 1;
        } else {
            quickSort(a, m + 1, ub);
            ub = m - 1;
        }
    }
}


int existValue(TYPE *array, int sizearray, TYPE value){
	int i;
	for(i=0;i<sizearray;i++){
		if(array[i]==value)
			return i;
	}
	return -1;
}



//void fill(T *a, int lb, int ub) {
//    int i;
//    srand(1);
//    for (i = lb; i <= ub; i++) a[i] = rand();
//}
//
//int main(int argc, char *argv[]) {
//    int maxnum, lb, ub;
//    T *a;
//
//    /* command-line:
//     *
//     *   qui maxnum
//     *
//     *   qui 2000
//     *       sorts 2000 records
//     *
//     */
//
//    maxnum = atoi(argv[1]);
//    lb = 0; ub = maxnum - 1;
//    if ((a = malloc(maxnum * sizeof(T))) == 0) {
//        fprintf (stderr, "insufficient memory (a)\n");
//        exit(1);
//    }
//
//    fill(a, lb, ub);
//    quickSort(a, lb, ub);
//
//    return 0;
//}
