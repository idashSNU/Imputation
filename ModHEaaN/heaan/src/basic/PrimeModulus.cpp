/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#include "basic/PrimeModulus.h"

namespace heaan
{
	
namespace basic
{

PrimeModulus::PrimeModulus(const word64 &prime) : prime__(prime)
{
	two_times_prime__ = prime__ << 1;
	computeBarrettParameter__();
	primitive_root__ = findPrimitiveRoot__();
}

void PrimeModulus::mult(const word64 &op1, const word64 &op2, const word64 &op2_PsInv, word64 &res)
{
	__heaan_mul_and_reduce_shoup(op1, op2, op2_PsInv, prime__, res);
}

void PrimeModulus::power(const word64 &op, word64 e, word64 &res)
{
	word64 tmp = op;
	res = 1;
	while (e > 0)
	{
		if (e & 1)
		{
			mult(tmp, res, res);
		}
		e = e >> 1;
		mult(tmp, tmp, tmp);
	}
}

void PrimeModulus::inverse(const word64 &op, word64 &res)
{
	word64 tmp = op > prime__ ? (op % prime__) : op;
	power(tmp, prime__ - 2, res);
}

void PrimeModulus::pseudoInverse(const word64 &op, word64 &res)
{
	word128 temp = static_cast<word128>(op) << 64;
	res = static_cast<word64>(temp / prime__);
}

void PrimeModulus::butterfly(word64 &op1, word64 &op2, const word64 &W, const word64 &W_PsInv)
{
	__heaan_butterfly(op1, op2, W, W_PsInv, prime__, two_times_prime__);
}

void PrimeModulus::butterflyInv(word64 &op1, word64 &op2, const word64 &W, const word64 &W_PsInv)
{
	__heaan_butterfly_inverse(op1, op2, W, W_PsInv, prime__, two_times_prime__);
}

void PrimeModulus::innerProduct(const vector<word64> &op1,
								const vector<word64> &op2, word64 &res)
{
	word64 high, low;
	__heaan_inner_product(op1, op2, high, low);
	__heaan_barrett_reduction(high, low, res, prime__, barrett_ratio__, barrett_k__);
}

void PrimeModulus::computeBarrettParameter__()
{
	barrett_k__ = 2 * floor(log2(prime__)) + 1;
	word64 high = static_cast<word64>(1) << (barrett_k__ - 64);
	__heaan_DIV_128_64(high, 0, prime__, barrett_ratio__);
}

word64 PrimeModulus::findPrimitiveRoot__()
{
	vector<word64> s;
	word64 phi = prime__ - 1;
	__heaan_find_prime_factors(s, phi);
	for (word64 r = 2; r <= phi; r++)
	{
		bool flag = false;
		for (auto it = s.begin(); it != s.end(); it++)
		{
			word64 tmp;
			power(r, phi / (*it), tmp);
			if (tmp % prime__ == 1)
			{
				flag = true;
				break;
			}
		}
		if (!flag)
		{
			return r;
		}
	}
	return -1;
}

} // namespace basic
} // namespace heaan
