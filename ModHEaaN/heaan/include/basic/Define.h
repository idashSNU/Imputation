/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_BASIC_DEFINE_H_
#define HEAAN_BASIC_DEFINE_H_

#include <vector>
#include <complex.h>
#include <stdint.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <utility>
#include <set>
#include <map>
#include <stdexcept>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

// Check uint128 //
#ifndef __SIZEOF_INT128__
#error "HEAAN library only works with __uint128_t"
#endif

typedef uint32_t word32; 

typedef long long unsigned int word64;

typedef __uint128_t word128;

typedef std::vector<word64> PrimePoly;

typedef std::vector<vector<word64>> PrimePolyVec;

typedef std::vector<complex<long double>> Message;


#include <immintrin.h>
#include <x86intrin.h>

inline void __heaan_MUL_64_64_high(word64 op1, word64 op2, word64 &high)
{
	const word128 _prod = static_cast<word128>(op1) * op2;
	high = (_prod >> 64) & 0xFFFFFFFFFFFFFFFF;
}

inline void __heaan_MUL_64_64_high(word64 op1, word64 op2, word64 *high)
{
	const word128 _prod = static_cast<word128>(op1) * op2;
	*high = (_prod >> 64) & 0xFFFFFFFFFFFFFFFF;
}

inline void __heaan_MUL_64_64_128(const word64 &op1, const word64 &op2, word64 &high,
								  word64 &low)
{
	const word128 _prod = static_cast<word128>(op1) * op2;
	low = (_prod) & 0xFFFFFFFFFFFFFFFF;
	high = (_prod >> 64) & 0xFFFFFFFFFFFFFFFF;
}

inline void __heaan_MUL_64_64_128(const word64 &op1, const word64 &op2, word64 *high,
								  word64 *low)
{
	const word128 _prod = static_cast<word128>(op1) * op2;
	*low = (_prod) & 0xFFFFFFFFFFFFFFFF;
	*high = (_prod >> 64) & 0xFFFFFFFFFFFFFFFF;
}


inline void __heaan_ADD_64_64_CARRY(word64 op1, word64 op2, word64 &res,
									unsigned char &carry)
{
	carry = _addcarry_u64(carry, op1, op2, &res);
}


inline void __heaan_DIV_128_64(word64 high, word64 low, word64 op,
							   word64 &res)
{
	word128 tmp = static_cast<word128>(high) << 64;
	tmp = tmp + static_cast<word128>(low);
	tmp = tmp / op;
	res = static_cast<word64>(tmp);
}

#endif
