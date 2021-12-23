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
#define sz(c) int((c).size())

const int N = 20000, M = 1e9;

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

template<typename DistribType>
void init_triangle_matrix(mt19937& gen, DistribType& distrib, vii& matrix) {
	for(int i = 0; i < sz(matrix); ++i)
                for(int j = i; j < sz(matrix); ++j)
        		matrix[i][j] = distrib(gen);
}

template<typename DistribType>
void init_band_matrix(mt19937& gen, DistribType& distrib, vii& matrix, int k) {
	for(int i = 0; i < sz(matrix); ++i) {
		for(int j = 0; j < sz(matrix); ++j)
			matrix[i][j] = (abs(i - j) > k) ? 0 : distrib(gen);
	}
}

void compute_band(vii& matrix, int k) {
	double ftime, itime, exec_time;
	int min_per_row, maxi;
	cout << "band" << endl;
        for(int p = 1; p <= 16; ++p) {
                //cout << "Количество потоков: " << p << endl;
		omp_set_num_threads(p);
		maxi = -M;
                itime = omp_get_wtime();
                #pragma omp parallel for shared(matrix) private(min_per_row) schedule(runtime) \
			reduction(max: maxi)
                for (int i = 0; i < N;  ++i) {
                        min_per_row = 0;
                        for(int j = max(i - k, 0); j <= min(N - 1, i + k); ++j)
                                min_per_row = min(min_per_row, matrix[i][j]);
			maxi = max(maxi, min_per_row);
                }
                ftime = omp_get_wtime();
                exec_time = ftime - itime;
                cout  << exec_time << endl;
        }
}

void compute_triangle(vii& matrix) {
	double ftime, itime, exec_time;
	int min_per_row, maxi;
	cout << "triangle" << endl;
	for(int p = 1; p <= 16; ++p) {
		//cout << "Количество потоков: " << p << endl;
		omp_set_num_threads(p);
		maxi = -M;
		itime = omp_get_wtime();
		#pragma omp parallel for shared(matrix) private(min_per_row) schedule(runtime) \
			reduction(max: maxi)
                for (int i = 0; i < N;  ++i) {
                        min_per_row = 0;
                        for(int j = i; j < N; ++j)
                                min_per_row = min(min_per_row, matrix[i][j]);
			maxi = max(maxi, min_per_row);
                        /*if (maxi < min_per_row) {
                                #pragma critical
                                {
                                        if (maxi < min_per_row)
                                                maxi = min_per_row;
                                }
                        }*/
                }
        	ftime = omp_get_wtime();
        	exec_time = ftime - itime;
        	cout << exec_time << endl;
	}
}

int main() {
        vii triangle_matrix(N), band_matrix(N);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(-M, 0);

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

	make_matrix(triangle_matrix);
	make_matrix(band_matrix);

	init_triangle_matrix<uniform_int_distribution<> >(gen, distrib, triangle_matrix);
	init_band_matrix<uniform_int_distribution<> >(gen, distrib, band_matrix, 1);

        //print(triangle_matrix);
	//print(band_matrix);

	compute_triangle(triangle_matrix);
	compute_band(band_matrix, 1);

	/*#pragma omp parallel
	{
		#pragma omp single
                {  cout << omp_get_num_threads() << endl;}
	}*/

        return 0;
}
