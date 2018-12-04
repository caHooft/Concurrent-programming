#include <fstream>
#include <iostream>
#include <complex>

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
		//int max_iteration = _rows * _cols;
		std::ofstream out(filename, std::ios_base::binary);
		out << "P6" << std::endl << _cols << " " << _rows << std::endl << 255 << std::endl;
		for (size_t y = 0; y < _rows; y++)
			for (size_t x = 0; x < _cols; x++)
				out << _matrix[y][x].r << _matrix[y][x].g << _matrix[y][x].b;
	}
};

int value(int x, int y, int width, int height) {
	complex<float> point((float)x / width - 1.5, (float)y / height - 0.5);
	complex<float> z(0, 0);
	int nb_iter = 0;
	while (abs(z) < 2 && nb_iter <= 20) {
		z = z * z + point;
		nb_iter++;
	}
	if (nb_iter < 20)
		return (255 * nb_iter) / 5;
	else
		return 0;
}

int main()
{
	const unsigned width = 400;
	const unsigned height = 400;

	ofstream my_Image("mandelbrot.ppm");

	if (my_Image.is_open()) {
		my_Image << "P3\n" << width << " " << height << " 255\n";
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++)
			{
				int val = value(i, j, width, height);
				my_Image << val << ' ' << 0 << ' ' << 0 << "\n";
			}
		}
		my_Image.close();
	}
	else
		cout << "Could not open the file";

	return 0;

}