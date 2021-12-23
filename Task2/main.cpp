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
using vll = vector<ll>;

#define sqr(x) ((x) * (x))
#define sz(c) ll((c).size())

const int N = 1ll << 24, M = 100;

void print(const vi& vec) {
        for(auto x: vec) {
                cout << x << " ";
        }
        cout << endl;
}

int main() {
        double ftime, itime, exec_time;
        ll total_scal = 0ll;
        vi vec1(N), vec2(N), temp(N);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(-M, M);
	auto generate_vector = [&distrib, &gen](auto& vec) {
		for_each(vec.begin(), vec.end(), [&distrib, &gen](int &n) {
			n = distrib(gen);
		});
	};

	generate_vector(vec1);
	generate_vector(vec2);

	for(int i = 0; i < N; ++i)
		temp[i] = vec1[i] * vec2[i];
	cout << "Верное скалярное: " << accumulate(temp.begin(), temp.end(), 0ll) << endl;

	itime = omp_get_wtime();
	for (int t = 0; t < 100; ++t) {
        total_scal = 0;
        #pragma omp parallel shared(vec1, vec2) reduction(+:total_scal)
        {
                /*#pragma omp single
                {  cout << omp_get_num_threads() << endl;}*/
                #pragma omp for
                for (int i = 0; i < N;  ++i)
                        total_scal = total_scal + vec1[i] * vec2[i];
        }
	}
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
        //cout << total_scal << " " << exec_time << endl;
	cout << total_scal << " " << exec_time << endl;
        return 0;
}
