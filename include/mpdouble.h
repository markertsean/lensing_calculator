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


#include <ctype.h>
#include <gmp.h>

class mpdouble
{
private:
	mpf_t data;
	static mpf_t _pi;
	static mpf_t _e;
	static int _base;
	static void compute_pi(unsigned long prec);
	static void compute_e(unsigned long prec);
	mpdouble(mpf_t op)
	{
		mpf_init_set(data, op);
	};
public:
	static void set_default_prec(unsigned long prec)
	{
		mpf_set_default_prec(prec);
		compute_pi(prec);
		compute_e(prec);
	};
	static void set_default_base(int base) {_base = base;};
	static inline mpdouble pi() {return mpdouble(_pi);};
	static inline mpdouble e() {return mpdouble(_e);};
	mpdouble()
	{
		mpf_init(data);
	};
	mpdouble(const mpdouble &op)
	{
		mpf_init_set(data, op.data);
	};
	mpdouble(const char *op)
	{
		mpf_init_set_str(data, op, _base);
	};
	mpdouble(unsigned long op)
	{
		mpf_init_set_ui(data, op);
	};
	mpdouble(long op)
	{
		mpf_init_set_si(data, op);
	};
	mpdouble(int op)
	{
		mpf_init_set_si(data, op);
	};
	mpdouble(double op)
	{
		mpf_init_set_d(data, op);
	};
	~mpdouble()
	{
		mpf_clear(data);
	};
	inline mpdouble &operator=(mpdouble const &op)
	{
		mpf_set(data, op.data);
		return *this;
	};
	inline mpdouble &operator=(unsigned long op)
	{
		mpf_set_ui(data, op);
		return *this;
	};
	inline mpdouble &operator=(long op)
	{
		mpf_set_si(data, op);
		return *this;
	};
	inline mpdouble &operator=(int op)
	{
		mpf_set_si(data, op);
		return *this;
	};
	inline mpdouble &operator=(double op)
	{
		mpf_set_d(data, op);
		return *this;
	};
	inline mpdouble &operator=(const char *op)
	{
		mpf_set_str(data, op, _base);
		return *this;
	};
	inline mpdouble &operator+() { return *this; };
	inline mpdouble operator-() const
	{
		mpdouble r;
		mpf_neg(r.data, data);
		return r;
	};
	inline mpdouble operator+(mpdouble const &op) const
	{
		mpdouble r;
		mpf_add(r.data, data, op.data);
		return r;
	};
	inline mpdouble operator+(unsigned long op) const
	{
		mpdouble r;
		mpf_add_ui(r.data, data, op);
		return r;
	};
	inline mpdouble operator+(long op) const
	{
		mpdouble r;
		mpf_add(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble operator+(int op) const
	{
		mpdouble r;
		mpf_add(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble operator+(double op) const
	{
		mpdouble r;
		mpf_add(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble& operator+=(mpdouble const &op)
	{
		return *this = *this + op;
	}
	inline mpdouble& operator+=(unsigned long op)
	{
		return *this = *this + op;
	}
	inline mpdouble& operator+=(int op)
	{
		return *this = *this + op;
	}
	inline mpdouble& operator+=(double op)
	{
		return *this = *this + op;
	}
	friend mpdouble operator+(unsigned long op1, mpdouble const &op2);
	friend mpdouble operator+(int op1, mpdouble const &op2);
	friend mpdouble operator+(double, mpdouble const &op2);
	inline mpdouble operator-(mpdouble const &op) const
	{
		mpdouble r;
		mpf_sub(r.data, data, op.data);
		return r;
	};
	inline mpdouble operator-(unsigned long op) const
	{
		mpdouble r;
		mpf_sub_ui(r.data, data, op);
		return r;
	};
	inline mpdouble operator-(long op) const
	{
		mpdouble r;
		mpf_sub(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble operator-(int op) const
	{
		mpdouble r;
		mpf_sub(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble operator-(double op) const
	{
		mpdouble r;
		mpf_sub(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble& operator-=(mpdouble const &op)
	{
		return *this = *this - op;
	}
	inline mpdouble& operator-=(unsigned long op)
	{
		return *this = *this - op;
	}
	inline mpdouble& operator-=(int op)
	{
		return *this = *this - op;
	}
	inline mpdouble& operator-=(double op)
	{
		return *this = *this - op;
	}
	friend mpdouble operator-(unsigned long op1, mpdouble const &op2);
	friend mpdouble operator-(int op1, mpdouble const &op2);
	friend mpdouble operator-(double, mpdouble const &op2);
	inline mpdouble operator*(mpdouble const &op) const
	{
		mpdouble r;
		mpf_mul(r.data, data, op.data);
		return r;
	};
	inline mpdouble operator*(unsigned long op) const
	{
		mpdouble r;
		mpf_mul_ui(r.data, data, op);
		return r;
	};
	inline mpdouble operator*(long op) const
	{
		mpdouble r;
		mpf_mul(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble operator*(int op) const
	{
		mpdouble r;
		mpf_mul(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble operator*(double op) const
	{
		mpdouble r;
		mpf_mul(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble& operator*=(mpdouble const &op)
	{
		return *this = *this * op;
	}
	inline mpdouble& operator*=(unsigned long op)
	{
		return *this = *this * op;
	}
	friend mpdouble operator*(unsigned long op1, mpdouble const &op2);
	friend mpdouble operator*(int op1, mpdouble const &op2);
	friend mpdouble operator*(double op1, mpdouble const &op2);
	inline mpdouble operator/(mpdouble const &op) const
	{
		mpdouble r;
		mpf_div(r.data, data, op.data);
		return r;
	};
	inline mpdouble operator/(unsigned long op) const
	{
		mpdouble r;
		mpf_div_ui(r.data, data, op);
		return r;
	};
	inline mpdouble operator/(long op) const
	{
		mpdouble r;
		mpf_div(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble operator/(int op) const
	{
		mpdouble r;
		mpf_div(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble operator/(double op) const
	{
		mpdouble r;
		mpf_div(r.data, data, mpdouble(op).data);
		return r;
	};
	inline mpdouble& operator/=(mpdouble const &op)
	{
		return *this = *this / op;
	}
	inline mpdouble& operator/=(unsigned long op)
	{
		return *this = *this / op;
	}
	friend mpdouble operator/(unsigned long op1, mpdouble const &op2);
	friend mpdouble operator/(int op1, mpdouble const &op2);
	friend mpdouble operator/(double op1, mpdouble const &op2);
	friend mpdouble sqrt(mpdouble const &op);
	friend mpdouble pow(mpdouble const &op1, unsigned long op2);
	friend mpdouble abs(mpdouble const &op);
	friend mpdouble floor(mpdouble const &op);
	inline bool operator<(mpdouble const &op) const
	{
		return mpf_cmp(data, op.data) < 0;
	};
	inline bool operator>(mpdouble const &op) const
	{
		return mpf_cmp(data, op.data) > 0;
	};
	inline bool operator<=(mpdouble const &op) const
	{
		return mpf_cmp(data, op.data) <= 0;
	};
	inline bool operator>=(mpdouble const &op) const
	{
		return mpf_cmp(data, op.data) >= 0;
	};
	inline bool operator==(mpdouble const &op) const
	{
		return mpf_cmp(data, op.data) == 0;
	};
	inline bool operator!=(mpdouble const &op) const
	{
		return mpf_cmp(data, op.data) != 0;
	};
	inline bool operator<(unsigned long op) const
	{
		return mpf_cmp_ui(data, op) < 0;
	};
	inline bool operator>(unsigned long op) const
	{
		return mpf_cmp_ui(data, op) > 0;
	};
	inline bool operator<=(unsigned long op) const
	{
		return mpf_cmp_ui(data, op) <= 0;
	};
	inline bool operator>=(unsigned long op) const
	{
		return mpf_cmp_ui(data, op) >= 0;
	};
	inline bool operator==(unsigned long op) const
	{
		return mpf_cmp_ui(data, op) == 0;
	};
	inline bool operator!=(unsigned long op) const
	{
		return mpf_cmp_ui(data, op) != 0;
	};
	inline bool operator<(long op) const
	{
		return mpf_cmp_si(data, op) < 0;
	};
	inline bool operator>(long op) const
	{
		return mpf_cmp_si(data, op) > 0;
	};
	inline bool operator<=(long op) const
	{
		return mpf_cmp_si(data, op) <= 0;
	};
	inline bool operator>=(long op) const
	{
		return mpf_cmp_si(data, op) >= 0;
	};
	inline bool operator==(long op) const
	{
		return mpf_cmp_si(data, op) == 0;
	};
	inline bool operator!=(long op) const
	{
		return mpf_cmp_si(data, op) != 0;
	};
	inline bool operator<(int op) const
	{
		return mpf_cmp_si(data, op) < 0;
	};
	inline bool operator>(int op) const
	{
		return mpf_cmp_si(data, op) > 0;
	};
	inline bool operator<=(int op) const
	{
		return mpf_cmp_si(data, op) <= 0;
	};
	inline bool operator>=(int op) const
	{
		return mpf_cmp_si(data, op) >= 0;
	};
	inline bool operator==(int op) const
	{
		return mpf_cmp_si(data, op) == 0;
	};
	inline bool operator!=(int op) const
	{
		return mpf_cmp_si(data, op) != 0;
	};
	inline bool operator<(double op) const
	{
		return mpf_cmp(data, mpdouble(op).data) < 0;
	};
	inline bool operator>(double op) const
	{
		return mpf_cmp(data, mpdouble(op).data) > 0;
	};
	inline bool operator<=(double op) const
	{
		return mpf_cmp(data, mpdouble(op).data) <= 0;
	};
	inline bool operator>=(double op) const
	{
		return mpf_cmp(data, mpdouble(op).data) >= 0;
	};
	inline bool operator==(double op) const
	{
		return mpf_cmp(data, mpdouble(op).data) == 0;
	};
	inline bool operator!=(double op) const
	{
		return mpf_cmp(data, mpdouble(op).data) != 0;
	};
	friend ::std::ostream &operator<<(::std::ostream &os, mpdouble const &op);
	friend ::std::istream &operator>>(::std::istream &is, mpdouble &op);
	inline double getdouble() { return mpf_get_d(data); };	/* No operator, to avoid unwanted implicit conversions */
};

::std::ostream &operator<<(::std::ostream &os, const mpdouble &op)
{
	mp_exp_t e;
	char *p = mpf_get_str(NULL, &e, mpdouble::_base, 0, op.data);
	int l = strlen(p);
	if (l == 0) os << '0';
	else if (e >= l)
	{
		if (p[0] == '-') e++;
		os.write(p, l);
		while (e-- > l) os << '0';
	}
	else if (e <= 0)
	{
		if (p[0] == '-') os << p[0];
		os << '.';
		while (e++ < 0) os << '0';
		os << (p[0] == '-' ? p + 1 : p);
	}
	else
	{
		if (p[0] == '-') e++;
		os.write(p, e);
		os << '.';
		os << p + e;
	}
	free(p);
	return os;
}

::std::istream &operator>>(::std::istream &is, mpdouble &op)
{
	char c;
	int n = 1000;
	char *b = (char *)malloc(n);
	char *p = b;
	do
	{
		if (is.eof()) goto DONE;
		is.read(&c, 1);
	} while (isspace(c));

	if (c == '+' || c == '-')
	{
		*p++ = c;
		is.read(&c, 1);
	}
	while (isdigit(c))
	{
		if (p-b >= n-1)
		{
			n += 1000;
			b = (char *)realloc(b, n);
		}
		*p++ = c;
		is.read(&c, 1);
		if (is.eof()) goto DONE;
	}
	if (c == '.')
	{
		if (p-b >= n-1)
		{
			n += 1000;
			b = (char *)realloc(b, n);
		}
		*p++ = c;
		is.read(&c, 1);
		if (is.eof()) goto DONE;
	}
	while (isdigit(c))
	{
		if (p-b >= n-1)
		{
			n += 1000;
			b = (char *)realloc(b, n);
		}
		*p++ = c;
		is.read(&c, 1);
		if (is.eof()) goto DONE;
	}
DONE:
	*p = '\0';
	op = b;
	free(b);
}

inline mpdouble operator+(unsigned long op1, mpdouble const &op2)
{
	return op2 + op1;
}

inline mpdouble operator+(double op1, mpdouble const &op2)
{
	return op2 + op1;
}

inline mpdouble operator+(int op1, mpdouble const &op2)
{
	return op2 + op1;
}

inline mpdouble operator-(unsigned long op1, mpdouble const &op2)
{
	mpdouble r;
	mpf_ui_sub(r.data, op1, op2.data);
	return r;
}

inline mpdouble operator-(double op1, mpdouble const &op2)
{
	mpdouble r;
	mpf_sub(r.data, mpdouble(op1).data, op2.data);
	return r;
}

inline mpdouble operator-(int op1, mpdouble const &op2)
{
	mpdouble r;
	mpf_sub(r.data, mpdouble(op1).data, op2.data);
	return r;
}

inline mpdouble operator*(unsigned long op1, mpdouble const &op2)
{
	return op2 * op1;
}

inline mpdouble operator*(double op1, mpdouble const &op2)
{
	return op2 * op1;
}

inline mpdouble operator*(int op1, mpdouble const &op2)
{
	return op2 * op1;
}

inline mpdouble operator/(unsigned long op1, mpdouble const &op2)
{
	mpdouble r;
	mpf_ui_div(r.data, op1, op2.data);
	return r;
}

inline mpdouble operator/(double op1, mpdouble const &op2)
{
	mpdouble r;
	mpf_div(r.data, mpdouble(op1).data, op2.data);
	return r;
}

inline mpdouble operator/(int op1, mpdouble const &op2)
{
	mpdouble r;
	mpf_div(r.data, mpdouble(op1).data, op2.data);
	return r;
}

inline mpdouble sqrt(mpdouble const &op)
{
	mpdouble r;
	mpf_sqrt(r.data, op.data);
	return r;
}

inline mpdouble pow(mpdouble const &op1, unsigned long op2)
{
	mpdouble r;
	mpf_pow_ui(r.data, op1.data, op2);
	return r;
}

inline mpdouble abs(mpdouble const &op)
{
	mpdouble r;
	mpf_abs(r.data, op.data);
	return r;
}

inline mpdouble floor(mpdouble const &op)
{
	mpdouble r;
	mpf_floor(r.data, op.data);
	return r;
}

mpf_t mpdouble::_pi;
mpf_t mpdouble::_e;
int mpdouble::_base = 10;

class mpdouble_init
{
public:
	mpdouble_init() {mpdouble::set_default_prec(1024);};
} mpdouble_initializer;

void mpdouble::compute_pi(unsigned long prec)
{
	mpdouble x = 6 - 4 * sqrt(mpdouble(2));
	mpdouble y = sqrt(mpdouble(2)) - 1;
	mpdouble two = 2;

	while (true)
	{
		two *= 4;
		mpdouble t = sqrt(sqrt(1 - y*y*y*y));
		mpdouble y_ = (1 - t) / (1 + t);
		t = 1 + y_;
		mpdouble x_ = (t*t*t*t * x) - two * y_ * (1 + y_ + y_*y_);
		if (x_ == x) break;
		x = x_;
		y = y_;
	}
	mpf_set_prec(_pi, prec);
	mpf_set(_pi, (1/x).data);
}

void mpdouble::compute_e(unsigned long prec)
{
	mpdouble s = 1;
	mpdouble e = 1;
	for (int i = 1; ; i++)
	{
		s /= i;
		mpdouble t = e + s;
		if (t == e) break;
		e = t;
	}
	mpf_set_prec(_e, prec);
	mpf_set(_e, e.data);
}

inline mpdouble log(mpdouble op)
{
    if (op < 0) return mpdouble(0);

    mpdouble r;
    while (op > 1.5)
    {
        r += 1;
        op /= mpdouble::e();
    }
    while (op < .5)
    {
        r -= 1;
        op *= mpdouble::e();
    }

    op = op - 1;
    mpdouble t = r;
    mpdouble o = 1;
    int s = -1;
    for (int i = 1; ; i++)
    {
        o *= op;
        s = -s;
        t += o / (s*i);
        if (t == r) break;
        r = t;
    }
    return r;
}

mpdouble sin(mpdouble x)
{
	mpdouble t = floor(x / (2*mpdouble::pi()));
	x -= t * 2 * mpdouble::pi();
	if (x > mpdouble::pi()) x -= 2 * mpdouble::pi();

	t = x;
	mpdouble p = x;
	mpdouble r = x;
	mpdouble xx = x*x;
	int i = 1;
	while (true)
	{
		// p *= -xx / ++i / ++i; // FAILS in GCC 4.x
		p *= -xx / ++i;
		p /= ++i;
		t = r + p;
		if (t == r) break;
		r = t;
	}
	return r;
}

inline mpdouble cos(mpdouble x)
{
	return sin(mpdouble::pi() / 2 - x);
}

mpdouble exp(mpdouble x)
{
	if (x == 0) return mpdouble(1);
	if (x < 0) return 1/exp(-x);
	mpdouble t = floor(x);
	x -= t;
	mpdouble p = pow(mpdouble::e(), (unsigned long)t.getdouble());
	mpdouble s = 1;
	mpdouble r = 1;
	for (int i = 1; ; i++)
	{
		s *= x / i;
		t = r + s;
		if (t == r) break;
		r = t;
	}
	return r * p;
}

mpdouble atan2(mpdouble y, mpdouble x)
{
	mpdouble r;
	if (abs(y) > abs(x)) r = mpdouble::pi() / 2 - atan2(x, y);
	else
	{
		if (abs(y) == abs(x)) r = mpdouble::pi() / 4;
		else if (abs(y) == 0) r = 0;
		else
		{
			r = abs(y/x);
			mpdouble yx = r*r;
			mpdouble s = r;
			for (int i = 3; ; i += 2)
			{
				s = -s * yx;
				mpdouble t = r + s/i;
				if (t == r) break;
				r = t;
			}
		}
		if (y < 0) r = -r;
		if (x < 0) r = mpdouble::pi() - r;
	}
	if (r < -mpdouble::pi()) r += 2*mpdouble::pi();
	if (r >= mpdouble::pi()) r -= 2*mpdouble::pi();
	return r;
}
