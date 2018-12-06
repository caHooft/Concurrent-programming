// mandelbrot.cpp
// compile with: g++ -std=c++11 mandelbrot.cpp -o mandelbrot
// view output with: eog mandelbrot.ppm
/*
#include <fstream>
#include <complex> // if you make use of complex number facilities in C++
#include <iostream>

using namespace std;

template <class T> struct RGB { T r, g, b; };

template <class T>
class Matrix {
public:
	Matrix(const size_t rows, const size_t cols) : _rows(rows), _cols(cols) {
		_matrix = new T*[rows];
		for (size_t i = 0; i < rows; ++i) {
			_matrix[i] = new T[cols];
		}
	}
	Matrix(const Matrix &m) : _rows(m._rows), _cols(m._cols) {
		_matrix = new T*[m._rows];
		for (size_t i = 0; i < m._rows; ++i) {
			_matrix[i] = new T[m._cols];
			for (size_t j = 0; j < m._cols; ++j) {
				_matrix[i][j] = m._matrix[i][j];
			}
		}
	}
	~Matrix() {
		for (size_t i = 0; i < _rows; ++i) {
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
	PPMImage(const size_t height, const size_t width) : Matrix(height, width) { }
	void save(const std::string &filename)
	{
		std::ofstream out(filename, std::ios_base::binary);
		out << "P6" << std::endl << _cols << " " << _rows << std::endl << 255 << std::endl;
		for (size_t y = 0; y < _rows; y++)
			for (size_t x = 0; x < _cols; x++)
				out << _matrix[y][x].r << _matrix[y][x].g << _matrix[y][x].b;
	}
};

void draw_Mandelbrot(PPMImage & image, const unsigned width, const unsigned height, double cxmin, double cxmax, double cymin, double cymax, unsigned int max_iterations)
{
	for (std::size_t ix = 0; ix < width; ++ix)
		for (std::size_t iy = 0; iy < height; ++iy)
		{
			std::complex<double> c(cxmin + ix / (width - 1.0)*(cxmax - cxmin), cymin + iy / (height - 1.0)*(cymax - cymin));
			std::complex<double> z = 0;
			unsigned int iterations;

			for (iterations = 0; iterations < max_iterations && std::abs(z) < 2.0; ++iterations)
				z = z * z + c;

			image[iy][ix].r = image[iy][ix].g = image[iy][ix].b = iterations;

		}
}

int main()
{
	const unsigned width = 1600;
	const unsigned height = 1600;

	PPMImage image(height, width);

	draw_Mandelbrot(image, width, height, -2.0, 0.5, -1.0, 1.0, 10);

	image.save("mandelbrot.ppm");
	return 0;
}
*/