// dotproduct.cpp
// Compile with:
// g++ -std=c++11 -pthread dotproduct.cpp -o dotproduct

#include <thread>
#include <vector>
#include <iostream>
#include <atomic>
using Vector = std::vector<int>;

struct DotProduct
{
	DotProduct(const Vector &a, const Vector &b) : a(a), b(b)
	{
		if (a.size() != b.size())
			throw "The vectors are of unequal length";
	}

	double operator()()
	{
		size_t nr_threads = 2;
		size_t length = a.size();

		int delta = length / nr_threads;
		int remainder = length % nr_threads;
		int L = 0, R = 0;

		result = 0;

		// Create threads
		for (int i = 0; i < nr_threads; ++i) {
			R = L + delta;
			if (i == nr_threads - 1)
				R += remainder;
			workers.push_back(std::thread(&DotProduct::partial_dot_product, this, L, R));
			L = R;
		}
		// Join threads
		for (auto &thread : workers) {
			thread.join();
		}

		return result;
	}

private:
	const Vector &a;
	const Vector &b;

	std::vector<std::thread> workers;
	std::atomic <int> result;

	void partial_dot_product(int L, int R)
	{
		for (int i = L; i < R; ++i) {
			result += a[i] * b[i];
		}
	}

};

int main()
{
	int nr_elements = 100000;

	// Fill two vectors with some values 
	Vector v1(nr_elements, 1), v2(nr_elements, 2);

	// Create Functor object
	DotProduct dp(v1, v2);

	// Print the result
	std::cout << dp() << std::endl;

	return 0;
}

