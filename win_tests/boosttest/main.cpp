#include <iostream>
#include <gsl/gsl_cblas.h>

int main() {
	double dirA[3] = {1,2,3};
	double dirC[3] = {1,2,3};
	cblas_dcopy(3, &dirA[0], 1, &dirC[0], 1);
	// cblas_dcopy(3, &c_xyz[0], 1, &dirC[0], 1);
	// cblas_daxpy(3, -1.0, &b_xyz[0], 1, &dirA[0], 1);
	// cblas_daxpy(3, -1.0, &b_xyz[0], 1, &dirC[0], 1);
	return 0;
}