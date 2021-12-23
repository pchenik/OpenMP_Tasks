#include <iostream>
#include <cmath>
#include <omp.h>
#include <vector>
#include <random>
#include <iomanip>
#include <cstdio>
#include <string>

using namespace std;
using ll = long long;
using vi = vector<int>;
using vll = vector<ll>;

#define sqr(x) ((x) * (x))
#define sz(c) ll((c).size())

const int N = 1 << 25;

void print(const vi& vec) {
	for(auto x: vec) {
		cout << x << " ";
	}
	cout << "\n";
}

int main() {
	double ftime, itime, exec_time;
	int max_per_thread, total_max = 0;
	vi arr(N);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distrib(-2 * N, 2 * N);
        for(auto& x: arr)
                x = distrib(gen);

	itime = omp_get_wtime();
	for(int t = 0; t < 30; ++t) {
	total_max = 0;
	#pragma omp parallel shared(arr, total_max) private(max_per_thread)
	{
		/*#pragma omp single
                {  cout << omp_get_num_threads() << endl;}*/
		max_per_thread = 0ll;
		#pragma omp for schedule(static)
		for (int i = 0; i < sz(arr);  ++i)
			if (arr[i] > max_per_thread) max_per_thread = arr[i];
		if (total_max < max_per_thread) {
			#pragma omp critical
			if (total_max < max_per_thread)
				total_max = max_per_thread;
		}
	}
	}
	ftime = omp_get_wtime();
	exec_time = ftime - itime;
	cout << " " << exec_time << endl;
        return 0;
}
