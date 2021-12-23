#include <iostream>
#include <cmath>
#include <omp.h>
#include <vector>
#include <random>
#include <iomanip>
#include <cstdio>
#include <string>
#include <algorithm>
#include <numeric>

using namespace std;
using ll = long long;
using vi = vector<int>;
using vii = vector<vi>;
using vll = vector<ll>;

#define sqr(x) ((x) * (x))
#define sz(c) ll((c).size())

const int N = 2e7, M = 100;

void print(const vi& vec) {
        for(auto x: vec) {
                cout << x << " ";
        }
        cout << endl;
}

void print(const vii& matrix) {
	for(auto vec: matrix) {
		for(auto x: vec)
			cout << x << " ";
		cout << endl;
	}
	cout << endl;
}

void check_reduction(vi& arr) {
	double ftime, itime, exec_time;
	long long sum;
	cout << "Reduction: " << endl;
	//itime = omp_get_wtime();
	 for(int p = 1; p <= 16; ++p) {
                omp_set_num_threads(p);
                sum = 0;
                itime = omp_get_wtime();
        	#pragma omp parallel for shared(arr) schedule(static) reduction(+: sum)
                	for (int i = 0; i < N;  ++i)
                        	sum += arr[i];
		ftime = omp_get_wtime();
        	exec_time = ftime - itime;
        	cout <<  exec_time << endl;
	}
}

void check_atomic(vi& arr) {
        double ftime, itime, exec_time;
	long long sum, sum_cur;
	cout << "Atomic: " << endl;
        //itime = omp_get_wtime();
	 for(int p = 1; p <= 16; ++p) {
                omp_set_num_threads(p);
                sum = 0;
                itime = omp_get_wtime();
                #pragma omp parallel shared(arr)
                {
                        #pragma omp for schedule(static)
                        for (int i = 0; i < N;  ++i) {
				#pragma omp atomic
                        		sum += arr[i];
			}
                }
		ftime = omp_get_wtime();
        	exec_time = ftime - itime;
        	cout <<  exec_time << endl;
        }
}

void check_locks(vi& arr) {
        double ftime, itime, exec_time;
	omp_lock_t lock;
	long long sum, sum_cur;
	cout << "Clocks: " << endl;
        //itime = omp_get_wtime();
	 for(int p = 1; p <= 16; ++p) {
                omp_set_num_threads(p);
                sum = 0;
                itime = omp_get_wtime();
		//omp_lock_t lock;
		omp_init_lock(&lock);
                #pragma omp parallel shared(arr)
		{
			//omp_set_lock(&lock);
			#pragma omp for schedule(static)
                        for (int i = 0; i < N;  ++i) {
				omp_set_lock(&lock);
                        	sum += arr[i];
                        	omp_unset_lock(&lock);
			}
			//omp_unset_lock(&lock);
		}
		omp_destroy_lock(&lock);
		ftime = omp_get_wtime();
        	exec_time = ftime - itime;
        	cout <<  exec_time << endl;
        }
}

void check_critical(vi& arr) {
	double ftime, itime, exec_time;
	long long sum, sum_cur;
	cout << "Critical: " << endl;
        //itime = omp_get_wtime();
        for(int p = 1; p <= 16; ++p) {
		omp_set_num_threads(p);
                sum = 0;
		itime = omp_get_wtime();
                #pragma omp parallel shared(arr)
                {
                        #pragma omp for schedule(static)
                        for (int i = 0; i < N;  ++i) {
				#pragma omp critical
                        		sum += arr[i];
			}
                }
		 ftime = omp_get_wtime();
        	 exec_time = ftime - itime;
       		 cout << exec_time << endl;
        }
}

int main() {
        vi arr(N);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(-M, M);
        auto generate_vector = [&distrib, &gen](auto& vec) {
                for_each(vec.begin(), vec.end(), [&distrib, &gen](int &n) {
                        n = distrib(gen);
                });
        };

	generate_vector(arr);

	#pragma omp parallel
	{
		#pragma omp single
                {  cout << omp_get_num_threads() << endl;}
	}

	check_reduction(arr);
	check_atomic(arr);
	check_locks(arr);
	check_critical(arr);
	//check_locks(arr);

        return 0;
}
