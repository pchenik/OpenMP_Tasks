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

//const int N = 10000, K = 10000, M = 1e9;
const int K = 100, N = 1000000, M = 1e9;
//const int N = 100, K = 1000000, M = 1e9;

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

int main() {
        double ftime, itime, exec_time;
        vii matrix(N);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(-M, 0);
        auto generate_vector = [&distrib, &gen](auto& vec) {
                for_each(vec.begin(), vec.end(), [&distrib, &gen](int &n) {
                        n = distrib(gen);
                });
        };


	for(auto& row: matrix) {
		row.resize(K);
		generate_vector(row);
	}

	//print(matrix);

 	int min_per_row1, maxi1 = -M;

	for (int i = 0; i < N;  ++i) {
			min_per_row1 = M;
                        for(int j = 0; j < K; ++j)
                                if (min_per_row1 > matrix[i][j])
                                    min_per_row1 =  matrix[i][j];
                        maxi1 = std::max(min_per_row1, maxi1);
        }
	cout << maxi1 << " actual answer " << endl;


	for (int p = 1; p <= 16; ++p) {
	vector<int> min_per_row(N, M);
	int maxi;
	omp_set_num_threads(p);
	itime = omp_get_wtime();
	//for(int t = 0; t < 10; ++t) {
        #pragma omp parallel for shared(matrix, min_per_row) schedule(static) collapse(2)
        for (int i = 0; i < N;  ++i)
		for(int j = 0; j < K; ++j) {
			if (matrix[i][j] < min_per_row[i]) {
                                #pragma critical
                                {
                                        if (matrix[i][j] < min_per_row[i])
                                                min_per_row[i] = matrix[i][j];
                                }
                        }
			//min_per_row[i] = min(min_per_row[i], matrix[i][j]);
		}
	maxi = -M;
	#pragma omp parallel for schedule(static) reduction(max:maxi)
	for(int i = 0; i < N; ++i)
		maxi = max(maxi, min_per_row[i]);
	//}
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        //cout << p << " " << maxi << " " << exec_time << endl;
	cout << " " << exec_time << endl;
	}
	//cout << " " << exec_time << endl;
        return 0;
}
