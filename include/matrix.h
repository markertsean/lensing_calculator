/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * Should you need to contact me, the author, you can do so by
 * e-mail - mail your message to Viktor T. Toth at <vttoth@vttoth.com>.
 */


#include <iostream>
#include <complex>

template <typename T> class matrix;

template<typename T> inline ::std::ostream &operator<<(::std::ostream &os, matrix<T> const &op)
{
	for (int r = 0; r < op.rows; r++)
	{
		os << '|';
		for (int c = 0; c < op.cols; c++)
		{
			os << ((c > 0) ? "," : "") << op.rdata[r][c];
		}
		os << '|' << endl;
	}
	return os;
}

template <typename T> class matrix
{
private:
	class row
	{
	public:
		T *data;
		row() {data = NULL;};
		row(int n) {data = new T[n];};
		~row() {delete[] data;};
		inline T &operator[](int c) {return data[c];};
	} *rdata;
	int rows, cols;
public:
	matrix(int r, int c)
	{
		rows = r;
		cols = c;
		rdata = new row[r];
		for (r = 0; r < rows; r++)
		{
			rdata[r].data = new T[c];
		};
	};
	matrix(const matrix &op)
	{
		rows = op.rows;
		cols = op.cols;
		rdata = new row[rows];
		for (int r = 0; r < rows; r++)
		{
			rdata[r].data = new T[cols];
			for (int c = 0; c < cols; c++)
			{
				rdata[r][c] = op.rdata[r][c];
			}
		}
	};
	~matrix() {delete[] rdata;};
	inline matrix &operator=(const matrix &op)
	{
		if (rows != op.rows || cols != op.cols) throw 0;
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				rdata[r][c] = op.rdata[r][c];
			}
		}
	};
	inline matrix sub(int r, int c)
	{
		matrix t(r, c);
		for (int i = 0; i < (r < rows ? r : rows); i++)
		{
			for (int j = 0; j < (c < cols ? c : cols); j++)
			{
				t[i][j] = rdata[i][j];
			}
		}
		return t;
	};
	inline row &operator[](int r) {return rdata[r];};
	inline matrix &operator*=(const T &op)
	{
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				rdata[r][c] *= op;
			}
		}
		return *this;
	};
	inline matrix operator*(const T &op)
	{
		matrix r(*this);
		r *= op;
		return r;
	};
	inline matrix &operator/=(const T &op)
	{
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				rdata[r][c] /= op;
			}
		}
		return *this;
	};
	inline matrix operator/(const T &op)
	{
		matrix r(*this);
		r /= op;
		return r;
	};
	inline matrix operator*(const matrix &op)
	{
		if (cols != op.rows) throw 0;
		matrix<T> t(rows, op.cols);
		for (int r = 0; r < t.rows; r++)
		{
			for (int c = 0; c < t.cols; c++)
			{
				t[r][c] = T(0);
				for (int i = 0; i < cols; i++)
				{
					t[r][c] += rdata[r][i] * op.rdata[i][c];
				}
			}
		}
		return t;
	}
	inline operator matrix< complex<T> > () const
	{
		matrix< complex<T> > t(rows, cols);
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				t[r][c] = rdata[r][c];
			}
		}
		return t;
	};
   	friend ::std::ostream &operator<< <T>(::std::ostream &os, matrix<T> const &op);
};
