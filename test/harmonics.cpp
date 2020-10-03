#include <memory>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
using namespace std;

unsigned N;

constexpr static inline unsigned ilog2(unsigned i, unsigned bits)
{	return !bits ? 0 : i >> bits ? bits + ilog2(i >> bits, bits >> 1) : ilog2(i, bits >> 1);
}
constexpr static inline unsigned ilog2(unsigned i)
{	return ilog2(i, numeric_limits<unsigned>::digits / 2);
}

class iteratablebitset
{	typedef unsigned long basetype;
	constexpr static const unsigned SHIFT = ilog2(numeric_limits<basetype>::digits);
	constexpr static const unsigned MASK = numeric_limits<basetype>::digits -1;
	unique_ptr<basetype[]> buffer;
 public:
	class const_iterator
	{	const basetype* root;
		unsigned cur;
	 public:
		const_iterator(const iteratablebitset& parent, unsigned cur) : root(parent.buffer.get()), cur(cur) {}
		unsigned operator*() const
		{	return cur;
		}
		void operator++();
		friend bool operator==(const const_iterator& l, const const_iterator& r)
		{	return l.cur == r.cur;
		}
		friend bool operator!=(const const_iterator& l, const const_iterator& r)
		{	return l.cur != r.cur;
		}
	};
	iteratablebitset() {}
	void reset(unsigned N) { buffer.reset(new basetype[(N+MASK)>>SHIFT]); memset(buffer.get(), 0, (N+MASK)>>SHIFT); }
	bool operator[](unsigned i) const
	{	return (bool)((buffer[i>>SHIFT] >> (i&MASK)) & 1U);
	}
	void set(unsigned i)
	{	buffer[i>>SHIFT] |= (basetype)1 << (i&MASK);
	}
	const_iterator begin() const
	{	const_iterator ret(*this, 0U);
		if (!(*this)[0U])
			++ret;
		return ret;
	}
	const_iterator end() const
	{	return const_iterator(*this, N);
	}
};
void iteratablebitset::const_iterator::operator++()
{	++cur;
	if (cur == N)
		return;
	const basetype* ptr = root + (cur >> SHIFT);
	unsigned bit = cur & MASK;
	if (bit == 0)
	{next:
		while (!*ptr) // fast forward
		{	cur += MASK + 1;
			if (cur >= N)
				goto end;
			++ptr;
	}	}
	{	basetype val = *ptr;
		do
		{	if (val & (basetype)1U << bit)
			{	cur = (cur & ~MASK) | bit;
				return;
			}
			++bit;
		} while (bit <= MASK);
	}
	cur = (cur | MASK) + 1;
	if (cur < N)
	{	++ptr;
		bit = 0;
		goto next;
	}
 end:
	cur = N;
}


iteratablebitset order1;
iteratablebitset order2;
iteratablebitset order3;

static inline bool check(const iteratablebitset& target, unsigned val)
{	return val < N && target[val];
}
static inline bool check(const iteratablebitset& target, unsigned num1, unsigned num2)
{	return check(target, num1 + num2)
		|| check(target, num1 - num2);
}
static bool check(unsigned num)
{	/*if (check(order3, 3U*num))
	{	//printf("%u: collision at %u in order3\n", num, 3U*num);
		return true;
	}
	for (unsigned i : order2)
		if (check(order3, num, i))
		{	//printf("%u: collision at %u+/-%u in order3\n", num, num, i);
			return true;
		}*/
	if (check(order2, 2U*num))
	{	printf("%u: collision at %u in order2\n", num, 2U*num);
		return true;
	}
	for (unsigned i : order1)
	{	if (check(order2, num, i))
		{	printf("%u: collision at %u-%u in order2\n", num, num, i);
			return true;
		}
		/*if (check(order3, 2U*num, i))
		{	//printf("%u: collision at %u+/-%u in order3\n", num, 2U*num, i);
			return true;
		}*/
	}
	/* implicitely false
	if (order1[num])
	{	printf("%u: collision in order1\n", num);
		return true;
	}*/
	return false;
}

static inline void apply(iteratablebitset& target, unsigned val)
{	if (val < N)
		target.set(val);
}
static inline void apply(iteratablebitset& target, unsigned num1, unsigned num2)
{	apply(target, num1 + num2);
	apply(target, num1 - num2);
}
static void apply(unsigned num)
{	/*apply(order3, 3U*num);
	for (auto i : order2)
		apply(order3, num, i);*/
	apply(order2, 2U*num);
	for (auto i : order1)
	{	apply(order2, num, i);
		//apply(order3, 2U*num, i);
	}
	order1.set(num);
}

int main(int argc, char* argv[])
{	unsigned i = 1;
	double finc = 1.;
	double fadd = 1.;
	switch (argc)
	{default:
		puts("usage: %s Nmax [Nmin [finc [fadd]]]");
		return 1;
	 case 5:
		fadd = strtod(argv[4], NULL);
	 case 4:
		finc = strtod(argv[3], NULL);
	 case 3:
		i = (unsigned)strtoul(argv[2], NULL, 0);
	 case 2:;
	}
	N = (unsigned)strtoul(argv[1], NULL, 0);
	order1.reset(N);
	order2.reset(N);
	order3.reset(N);
	
	unsigned count = 0;
	while (i < N)
	{	if (!check(i))
		{	printf("HIT%u: %u\n", ++count, i);
			apply(i);
			i = (unsigned)(finc * i + fadd);
		} else
			++i;
	}
	
	return 0;
}
