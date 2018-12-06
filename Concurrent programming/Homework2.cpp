// mandelbrot.cpp
// compile with: g++ -std=c++11 mandelbrot.cpp -o mandelbrot
// view output with: eog mandelbrot.ppm

#include <fstream>
#include <complex> // if you make use of complex number facilities in C++
#include <iostream>
#include <thread>
#include <mutex>
#include <vector>

using namespace std;

template <class T> struct RGB { T r, g, b; };

template <class T>
class Matrix 
{
public:
	Matrix(const size_t rows, const size_t cols) : _rows(rows), _cols(cols) 
	{
		_matrix = new T*[rows];
		for (size_t i = 0; i < rows; ++i) 
		{
			_matrix[i] = new T[cols];
		}
	}
	Matrix(const Matrix &m) : _rows(m._rows), _cols(m._cols) 
	{
		_matrix = new T*[m._rows];
		for (size_t i = 0; i < m._rows; ++i) 
		{
			_matrix[i] = new T[m._cols];
			for (size_t j = 0; j < m._cols; ++j) 
			{
				_matrix[i][j] = m._matrix[i][j];
			}
		}
	}
	~Matrix() {
		for (size_t i = 0; i < _rows; ++i) 
		{
			delete[] _matrix[i];
		}
		delete[] _matrix;
	}
	T *operator[] (const size_t nIndex)
	{
		return _matrix[nIndex];
	}
	size_t width() const { return _cols; }
	size_t height() const { return _rows; }
protected:
	size_t _rows, _cols;
	T **_matrix;
};

// Portable PixMap image
class PPMImage : public Matrix<RGB<unsigned char> >
{
public:
	unsigned int size;

	PPMImage(const size_t height, const size_t width) : Matrix(height, width) { }
	void save(const string &filename)
	{
		ofstream out(filename, ios_base::binary);
		out << "P6" << endl << _cols << " " << _rows << endl << 255 << endl;
		for (size_t y = 0; y < _rows; y++)
			for (size_t x = 0; x < _cols; x++)
				out << _matrix[y][x].r << _matrix[y][x].g << _matrix[y][x].b;
	}
};

void draw_Mandelbrot(PPMImage & image, const unsigned width, const unsigned height, double cxmin, double cxmax, double cymin, double cymax, unsigned int max_iterations)
{
	for (size_t ix = 0; ix < width; ++ix)
		for (size_t iy = 0; iy < height; ++iy)
		{
			complex<double> c(cxmin + ix / (width - 1.0)*(cxmax - cxmin), cymin + iy / (height - 1.0)*(cymax - cymin));
			complex<double> z = 0;
			unsigned int iterations;

			for (iterations = 0; iterations < max_iterations && abs(z) < 2.0; ++iterations)
				z = z * z + c;

			image[iy][ix].r = image[iy][ix].g = image[iy][ix].b = iterations;

		}
}

int main()
{
	const unsigned width = 512;
	const unsigned height = 512;

	PPMImage image(height, width);


	int parts = 4;

	vector<int>bnd(parts, image.size);

	thread *tt = new thread[parts - 1];

	time_t start, end;
	time(&start);
	//Lauch parts-1 threads
	for (int i = 0; i < parts - 1; ++i)
	{
		tt[i] = thread(draw_Mandelbrot, ref(image), width, height, -2.0, 0.5, -1.0, 1.0, 30);
	}

	//Use the main thread to do part of the work !!!
	for (int i = parts - 1; i < parts; ++i)
	{
		draw_Mandelbrot(ref(image), width, height, -2.0, 0.5, -1.0, 1.0, 30);
	}

	//Join parts-1 threads
	for (int i = 0; i < parts - 1; ++i)
		tt[i].join();

	time(&end);
	cout << difftime(end, start) << " seconds" << endl;


	image.save("mandelbrot.ppm");

	delete[] tt;

	return 0;
}