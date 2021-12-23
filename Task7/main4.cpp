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
#include <fstream>

using namespace std;
using ll = long long;
using vi = vector<int>;
using vvi = vector<vi>;
using vvl = vector<ll>;

#define sqr(x) ((x) * (x))
#define sz(c) int((c).size())

//const int N = 100, M = 100000, K = 10;
//const int N = 1e7, M = 10, K = 10;
//const int N = 6, M = 10, K = 10;
//const int N = 6, M = 20, K = 10;

const int K = 10;

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

int main() {
	int N, M;
	ios_base::sync_with_stdio(0);
	cin.tie(0);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(0, K);

	freopen("input5.txt", "r", stdin);
	ofstream myfile("output2.txt");
        cin >> M >> N;

	vvi vs(M);
        auto generate_vector = [&distrib, &gen](auto& vec) {
                for_each(vec.begin(), vec.end(), [&distrib, &gen](int &n) {
                        n = distrib(gen);
                });
        };

        auto make_matrix = [=](auto& matrix) {
                for_each(matrix.begin(), matrix.end(), [=](auto& matrix_row) {
                        matrix_row.resize(N);
                });
        };

	#pragma omp parallel
	{
		#pragma omp single
                {  cout << omp_get_num_threads() << " threads for consequtive" << endl;}
	}

	ll sum = 0, cur_sum;
	make_matrix(vs);
	//omp_set_nested(1);

	double ftime, itime, exec_time;
        itime = omp_get_wtime();
	//#pragma omp parallel for schedule(static) shared(vs) reduction(+:sum) private(cur_sum)
	for(int i = 0; i < M; ++i) {
		cur_sum = 0;
		for(int j = 0; j < N; ++j)
			cin >> vs[i][j];
		for(int j = 0; j < N / 2; ++j)
			cur_sum += vs[i][j] * vs[i][j + N / 2];
		sum = max(cur_sum, sum);
		myfile << sum << " i = " << i + 1 << endl;
		//sum += cur_sum;
	}
	ftime = omp_get_wtime();
        exec_time = ftime - itime;
	myfile.close();
        cout << sum << " " << exec_time << endl;
        return 0;
}
