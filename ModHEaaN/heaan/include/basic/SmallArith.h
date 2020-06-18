/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_BASIC_SMALLARITH_H_
#define HEAAN_BASIC_SMALLARITH_H_

#include "Define.h"

using namespace std;

static void __heaan_barrett_reduction(const word64 &high, const word64 &low,
									  word64 &res, const word64 &prime, const word64 &ratio, long k)
{
	unsigned char carry = 0;
	word64 high1;
	word64 high2, low2;
	__heaan_MUL_64_64_high(low, ratio, high1);
	__heaan_MUL_64_64_128(high, ratio, high2, low2);
	__heaan_ADD_64_64_CARRY(high1, low2, high1, carry);
	__heaan_ADD_64_64_CARRY(high2, 0, high2, carry);
	high1 >>= k - 64;
	high2 <<= 128 - k;
	high1 = high1 + high2;
	high1 = high1 * prime;
	res = low - high1;
	if (res >= prime)
	{
		res -= prime;
	}
}

static void __heaan_barrett_reduction(const word64& op, word64& res,
    const word64& prime, const word64& ratio, long k)
{
    word64 high;
    __heaan_MUL_64_64_high(op, ratio, &high);
    high >>= k - 64;
    res = op - (high * prime);
    if(res >= prime)
    {
        res -= prime;
    }
}

inline static void __heaan_mul_and_reduce_shoup(const word64 &op1, const word64 &op2, const word64 &op2_PsInv, const word64 &prime, word64 &res)
{
	word64 high;
	__heaan_MUL_64_64_high(op2_PsInv, op1, high);
	res = op1 * op2 - high * prime;
}

/// modulus reductionfor q = k * 2^m + 1 (the result is a * b * k mod q)
inline static void __heaan_mul_and_reduce_k(const word64 &op1, const word64 &op2_kInv, long k, long m, word64 &res)
{
	word64 high, low;
	__heaan_MUL_64_64_128(op1, op2_kInv, high, low);
	word64 c0 = low & ((1 << m) - 1);
	word64 c1 = (low >> m) + (high << (64 - m));
	res = k * c0 + c1;
}

static void __heaan_inner_product(const vector<word64> &op1, const vector<word64> &op2, word64 &high, word64 &low)
{
#ifdef DEBUG
	if (op1.size() != op2.size())
	{
		cerr << "__inner_product : The size of two input vector should be same" << endl;
	}
#endif
	unsigned char carry = 0;
	word64 tmp_high, tmp_low;
	for (size_t i = 0; i < op1.size(); ++i)
	{
		__heaan_MUL_64_64_128(op1[i], op2[i], tmp_high, tmp_low);
		carry = 0;
		__heaan_ADD_64_64_CARRY(tmp_low, low, low, carry);
		__heaan_ADD_64_64_CARRY(tmp_high, high, high, carry);
	}
}

inline static void __heaan_butterfly(word64 &op1, word64 &op2, const word64 &W, const word64 &W_PsInv, const word64 &prime, const word64 &two_prime)
{
	if (op1 >= two_prime)
		op1 -= two_prime;
	word64 U;
	__heaan_mul_and_reduce_shoup(op2, W, W_PsInv, prime, U);
	op2 = op1 + two_prime - U; // op2 <- op1 - op2 * W
	op1 += U;				   // op1 <- op1 + op2 * W
}

inline static void __heaan_butterfly_inverse(word64 &op1, word64 &op2, const word64 &W, const word64 &W_PsInv, const word64 &prime, const word64 &two_prime)
{
	uint64_t T = op1 + two_prime - op2;
	op1 += op2;
	if (op1 >= two_prime)
		op1 -= two_prime;
	__heaan_mul_and_reduce_shoup(T, W, W_PsInv, prime, op2);
}

void static __heaan_mult_mod(const word64 &op1, const word64 &op2, const word64 &mod, word64 &res)
{
	word128 tmp = static_cast<word128>(op1) * op2;
	res = tmp % mod;
}

void static __heaan_power_mod(const word64 &op, const word64 &mod, word64 e, word64 &res)
{
	word64 tmp = op;
	res = 1;
	while (e > 0)
	{
		if (e & 1)
		{
			__heaan_mult_mod(tmp, res, mod, res);
		}
		e = e >> 1;
		__heaan_mult_mod(tmp, tmp, mod, tmp);
	}
}

void static __heaan_inverse_mod(const word64 &op, const word64 &mod, word64 &res)
{
	word64 tmp = op > mod ? (op % mod) : op;
	__heaan_power_mod(tmp, mod, mod - 2, res);
}

bool static __heaan_is_prime(const word64 p)
{
	if (p < 2)
		return false;
	if (p != 2 && p % 2 == 0)
		return false;
	word64 s = p - 1;
	while (s % 2 == 0)
	{
		s /= 2;
	}
	for (long i = 0; i < 200; i++)
	{
		word64 temp1 = rand();
		temp1 = (temp1 << 32) | rand();
		temp1 = temp1 % (p - 1) + 1;
		word64 temp2 = s;
		word64 mod;
		__heaan_power_mod(temp1, p, temp2, mod);
		while (temp2 != p - 1 && mod != 1 && mod != p - 1)
		{
			__heaan_mult_mod(mod, mod, p, mod);
			temp2 *= 2;
		}
		if (mod != p - 1 && temp2 % 2 == 0)
			return false;
	}
	return true;
}

#endif
