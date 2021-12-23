#include <cmath>
#include <iostream>
#include <iomanip>
#include <omp.h>

using namespace std;

long double f(long double x) {
	return 1 / (1 + x * x);
}

long double f2(long double x) {
	return exp(x);
}

const long double A = 0, B = 1;

int main() {
	long double h, integral = 0, itime, ftime, exec_time;
	int n = 1500000000;
	h = (B - A) / double(n);

	itime = omp_get_wtime();
	#pragma omp parallel reduction(+:integral)
	{
		#pragma omp single
		{int k = omp_get_num_threads(); cout << k << endl;}
		#pragma omp for schedule(static)
		for (int i = 0; i <= n - 1; ++i)
			integral += f2(h * i);
	}
	integral *= h;
	ftime = omp_get_wtime();
	exec_time = ftime - itime;
	cout << fixed << setprecision(16) << integral << " " << exec_time << endl;
	return 0;
}
