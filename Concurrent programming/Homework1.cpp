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
int findMandelbrot(double cr, double ci, int max_iterations)
{
	int i = 0;
	double zr = 0.0, zi = 0.0;
	while (i < max_iterations && zr * zr + zi * zi < 4.0)
	{
		double temp = zr * zr - zi * zi + cr;
		zi = 2.0 * zr * zi + ci;
		zr = temp;
		i++;
	}

	return i;
}

double mapToReal(int x, int imageWidth, double minR, double maxR)
{
	double range = maxR - minR;
	return x * (range / imageWidth) + minR;
}

double mapToImaginary(int y, int imageHeight, double minI, double maxI)
{
	double range = maxI - minI;
	return y * (range / imageHeight) + minI;
}

int main()
{
	const unsigned width = 512;
	const unsigned height = 512;
	int maxN = 255;
	double minR = -1.5, maxR = 0.7, minI = -1, maxI = 1;
	PPMImage image(height, width);	

	// Open the output file, write the PPM header...
	ofstream fout("output_image.ppm");
	fout << "P3" << endl; // "Magic Number" - PPM file
	fout << width << " " << height << endl; // Dimensions
	fout << "255" << endl; // Maximum value of a pixel R,G,B value...

	// For every pixel...
	for (int y = 0; y < height; y++) // Rows...
	{
		for (int x = 0; x < width; x++) // Pixels in row (columns)...
		{
			// ... Find the real and imaginary values for c, corresponding to that
			//     x, y pixel in the image.
			double cr = mapToReal(x, width, minR, maxR);
			double ci = mapToImaginary(y, height, minI, maxI);

			// ... Find the number of iterations in the Mandelbrot formula
			//     using said c.
			int n = findMandelbrot(cr, ci, maxN);

			// ... Map the resulting number to an RGP value
			int r = (n % 256);
			int g = (n % 256);
			int b = (n % 256);

			// ... Output it to an image
			fout << r << " " << g << " " << b << " ";
		}
		fout << endl;
	}
	fout.close();

	image.save("mandelbrot.ppm");
	return 0;
}