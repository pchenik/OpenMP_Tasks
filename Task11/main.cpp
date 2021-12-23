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
using vvi = vector<vi>;
using vll = vector<ll>;

#define sqr(x) ((x) * (x))
#define sz(c) ll((c).size())

const int N = 1100000, M = 100;

void print(const vi& vec) {
        for(auto x: vec) {
                cout << x << " ";
        }
        cout << endl;
}

void print(const vvi& matrix) {
	for(auto vec: matrix) {
		for(auto x: vec)
			cout << x << " ";
		cout << endl;
	}
	cout << endl;
}

template<typename DistribType>
void init_matrix(mt19937& gen, DistribType& distrib, vvi& matrix) {
        for(int i = 0; i < sz(matrix); ++i)
                for(int j = 0; j < sz(matrix[i]); ++j)
                        matrix[i][j] = distrib(gen);
}


template<typename DistribType>
void check_dynamic(const vi& a,  DistribType& distrib, mt19937& gen) {
	double ftime, itime, exec_time;
	int k, sum;
	//ll sum;
	cout << "check_dynamic" << endl;
	for (int p = 1; p <= 16; ++p) {
	vi kk(N);
	for(auto& x: kk) {
		x = distrib(gen) % 15 + 1;
		//cout << x << " ";
	}
	cout << endl << "here " << kk[0] << " " << " " << kk[1] << " " << kk[2] << " " << kk[3] << endl;
	omp_set_num_threads(p);
	//sum = 0;
	itime = omp_get_wtime();
	#pragma omp parallel for shared(a, kk) private(k, sum) schedule(static)
                for(int i = 0; i < N; ++i) {
                        k = 1ll << kk[i];
			sum = 0;
			for (int j = 0; j < k; ++j)
				sum += a[j];
		}
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout  << exec_time << endl;

	//sum = 0;
        itime = omp_get_wtime();
	#pragma omp parallel for shared(a, kk) private(k, sum) schedule(dynamic)
                for(int i = 0; i < N; ++i) {
			k = 1 << kk[i];
                        sum = 0;
                        for (int j = 0; j < k; ++j)
                                sum += a[j];
                }
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout  << exec_time << endl;

	//sum = 0;
        itime = omp_get_wtime();
	#pragma omp parallel for shared(a, kk) private(k, sum) schedule(guided)
                for(int i = 0; i < N; ++i) {
			k = 1 << kk[i];
                        sum = 0;
                        for (int j = 0; j < k; ++j)
                                sum += a[j];
                }
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout  << exec_time << endl;
	cout << endl;
	}
}

void check_static(const vi& a) {
        double ftime, itime, exec_time;
        ll sum;
	cout << "check_static" << endl;
	for (int p = 1; p <= 16; ++p) {
        omp_set_num_threads(p);
        sum = 0;
	itime = omp_get_wtime();
        #pragma omp parallel for schedule(static) reduction(+:sum)
                for(int i = 0; i < N; ++i)
                        sum += a[i];
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout <<  exec_time << endl;

        itime = omp_get_wtime();
        sum = 0;
        #pragma omp parallel for schedule(dynamic) reduction(+:sum)
                for(int i = 0; i < N; ++i)
                        sum += a[i];
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout <<  exec_time << endl;

        itime = omp_get_wtime();
        sum = 0;
        #pragma omp parallel for schedule(guided) reduction(+:sum)
                for(int i = 0; i < N; ++i)
                        sum += a[i];
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout << exec_time << endl << endl;
	}
}
void check_guided(const vi& a) {
        double ftime, itime, exec_time;
       	const ll max = 32, factor = 10000000l;
	cout << "check_quided" << endl;
	for (int p = 1; p <= 16; ++p) {
        omp_set_num_threads(p);
        itime = omp_get_wtime();
    	#pragma omp parallel for schedule(static)
    	for (int i = 0; i < max; i++) {
		ll sum = 0;
    		for (int w = 0; w < (max - i) * factor; w++) sum += w;
    	}
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout << exec_time << endl;

        itime = omp_get_wtime();
	#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < max; i++) {
                ll sum = 0;
                for (int w = 0; w < (max - i) * factor; w++) sum += w;
        }
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout << exec_time << endl;

        itime = omp_get_wtime();
	#pragma omp parallel for schedule(guided, 1)
        for (int i = 0; i < max; i++) {
                ll sum = 0;
                for (int w = 0; w < (max - i) * factor; w++) sum += w;
        }
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout << exec_time << endl << endl;
	}
}

void check_guided2(const vi& a) {
        double ftime, itime, exec_time;
        const ll max = 32, factor = 10000000l;
        cout << "check_quided" << endl;
        for (int p = 1; p <= 16; ++p) {
        omp_set_num_threads(p);
        itime = omp_get_wtime();
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < max; i++) {
                ll sum = 0;
                for (int w = 0; w < (max - i) * factor; w++) sum += w;
        }
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout << exec_time << endl;

        itime = omp_get_wtime();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < max; i++) {
                ll sum = 0;
                for (int w = 0; w < (max - i) * factor; w++) sum += w;
        }
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout << exec_time << endl;

        itime = omp_get_wtime();
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i < max; i++) {
                ll sum = 0;
                for (int w = 0; w < (max - i) * factor; w++) sum += w;
        }
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        cout << exec_time << endl << endl;
        }
}

int main() {
        double ftime, itime, exec_time;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(0, M);
        auto generate_vector = [&distrib, &gen](auto& vec) {
                for_each(vec.begin(), vec.end(), [&distrib, &gen](int &n) {
                        n = distrib(gen);
                });
        };

	auto make_matrix = [=](auto& matrix) {
                for_each(matrix.begin(), matrix.end(), [](auto& matrix_row) {
                        matrix_row.resize(N);
                });
        };
	#pragma omp parallel
	{
		#pragma omp single
		cout << omp_get_num_threads() << endl;
	}
	vi a(N);
	//vi a_static(N), a_guided(N), a_dynamic(N);
	generate_vector(a);
	//generate_vector(a_static);
	//generate_vector(a_guided);
 	//generate_vector(a_dynamic);

	//Проверка static, dynamic, guided
	check_static(a);
	geometric_distribution<> distrib2(0.2);
	check_dynamic<geometric_distribution<> >(a, distrib2, gen);
	check_guided(a);
        return 0;
}
