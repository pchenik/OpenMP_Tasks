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
using vvl = vector<ll>;

#define sqr(x) ((x) * (x))
#define sz(c) int((c).size())

//const int N = 100, M = 100000, K = 10;
const int N = 100, M = 100, K = 10;
//const int N = 1e7, M = 10, K = 10;

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
	vvi vs(M);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(0, K);

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

	make_matrix(vs);
	init_matrix<uniform_int_distribution<> >(gen, distrib, vs);
	cout << M << " " << N << endl;
        print(vs);
        return 0;
}
