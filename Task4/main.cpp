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

const int N = 10000, M = 1e9;

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
		row.resize(N);
		generate_vector(row);
	}

        //print(matrix);

	int min_per_row, maxi = -M;
	/*#pragma omp parallel
	{
		#pragma omp single
                {  cout << omp_get_num_threads() << endl;}
	}*/

	itime = omp_get_wtime();
	for(int t = 0; t < 1; ++ t) {
	maxi = -M; 
        #pragma omp parallel for shared(matrix, maxi) private(min_per_row) schedule(static)
                for (int i = 0; i < N;  ++i) {
			min_per_row = M;
			for(int j = 0; j < N; ++j)
				min_per_row = min(min_per_row, matrix[i][j]);
			if (maxi < min_per_row) {
				#pragma critical
				{
					if (maxi < min_per_row)
						maxi = min_per_row;
				}
			}
		}
	}
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        //cout << maxi << " " << exec_time << endl;
	cout << maxi << " " << exec_time << endl;
        return 0;
}
