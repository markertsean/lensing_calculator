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


#include <complex>
#include "mpdouble.h"

using namespace std;

typedef complex<mpdouble> mpcomplex;

inline mpcomplex operator+(mpcomplex const &op1, int op2)
{
	return op1 + mpdouble(op2);
}

mpcomplex operator+(mpcomplex const &op1, double op2)
{
	return op1 + mpdouble(op2);
}

inline mpcomplex &operator+=(mpcomplex &op1, int op2)
{
	return op1 += mpcomplex(op2, 0);
}

inline mpcomplex operator-(mpcomplex const &op1, int op2)
{
	return op1 - mpdouble(op2);
}

inline mpcomplex operator-(mpcomplex const &op1, double op2)
{
	return op1 - mpdouble(op2);
}

inline mpcomplex &operator-=(mpcomplex &op1, int op2)
{
	return op1 -= mpcomplex(op2, 0);
}

inline mpcomplex operator*(mpcomplex const &op1, int op2)
{
	return op1 * mpdouble(op2);
}

inline mpcomplex operator*(mpcomplex const &op1, double op2)
{
	return op1 * mpdouble(op2);
}

inline mpcomplex &__doadv(mpcomplex *op1, mpcomplex const &op2)
{
	mpdouble D = op2.real() * op2.real() + op2.imag() * op2.imag();
	mpdouble Re = (op1->real() * op2.real() + op1->imag() * op2.imag()) / D;
	mpdouble Im = (op1->imag() * op2.real() - op1->real() * op2.imag()) / D;
	*op1 = mpcomplex(Re, Im);
	return *op1;
}

inline mpcomplex operator/(mpcomplex const &op1, mpcomplex const &op2)
{
	mpcomplex t(op1);
	t /= op2;
	return t;
}

inline mpcomplex operator/(mpcomplex const &op1, int op2)
{
	return op1 / mpdouble(op2);
}

inline mpcomplex operator/(mpcomplex const &op1, double op2)
{
	return op1 / mpdouble(op2);
}

inline mpcomplex operator/(int op1, mpcomplex const &op2)
{
	return mpcomplex(op1, 0) / op2;
}

inline mpcomplex operator/(double op1, mpcomplex const &op2)
{
	return mpcomplex(op1,0) / op2;
}

inline ::std::ostream &operator<<(::std::ostream &os, mpcomplex const &op)
{
	os << '(' << op.real() << ',' << op.imag() << ')';
	return os;
}

inline mpcomplex log(mpcomplex op)
{
	return mpcomplex(0.5 * log(op.real()*op.real()+op.imag()*op.imag()),
					 atan2(op.imag(), op.real()));
}

inline mpcomplex exp(mpcomplex op)
{
	return mpcomplex(cos(op.imag()), sin(op.imag())) * exp(op.real());
}

inline mpcomplex sin(mpcomplex op)
{
	mpdouble ex = exp(op.imag());
	mpdouble sx = sin(op.real());
	mpdouble cx = cos(op.real());
	return mpcomplex(sx*(ex+1/ex)/2, cx*(ex-1/ex)/2);
}
