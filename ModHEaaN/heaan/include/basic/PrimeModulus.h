/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_BASIC_PRIMEMODULUS_H_
#define HEAAN_BASIC_PRIMEMODULUS_H_

#include "Tools.h"
#include "SmallArith.h"

namespace heaan
{
namespace basic
{

class PrimeModulus
{

	friend class PrimeRing;

	friend class Ring;

public:

	PrimeModulus() = default;

	PrimeModulus(const word64 &prime);

	word64 getPrimitiveRoot() const
	{
		return primitive_root__;
	}

	word64 getPrime() const
	{
		return prime__;
	}

	word64 getBarrettRatio() const
	{
		return barrett_ratio__;
	}

	long getBarrettK() const
	{
		return barrett_k__;
	}

	inline void add(const word64 &op1, const word64 &op2, word64 &res)
	{
		res = op1 + op2;
		if (res > prime__)
		{
			res -= prime__;
		}
	}

	inline void sub(const word64 &op1, const word64 &op2, word64 &res)
	{
		if (op1 >= op2)
		{
			res = op1 - op2;
			return;
		}
		else
		{
			res = op2 - op1;
			res = prime__ - res;
			return;
		}
	}

	inline void mult(const word64 &op1, const word64 &op2, word64 &res)
	{
		word64 high, low;
		__heaan_MUL_64_64_128(op1, op2, high, low);
		__heaan_barrett_reduction(high, low, res, prime__, barrett_ratio__, barrett_k__);
	}

	void mult(const word64 &op1, const word64 &op2, const word64 &op2_PsInv, word64 &res);

	void power(const word64 &op, word64 e, word64 &res);

	void inverse(const word64 &op, word64 &res);

	void pseudoInverse(const word64 &op, word64 &res);

	void butterfly(word64 &op1, word64 &op2, const word64 &W, const word64 &W_PsInv);

	void butterflyInv(word64 &op1, word64 &op2, const word64 &W, const word64 &W_PsInv);

	void innerProduct(const vector<word64> &op1, const vector<word64> &op2, word64 &res);

private:

	word64 prime__ = 0; ///< a prime integer p for this integer ring

	word64 two_times_prime__ = 0; ///< = 2 * prime_ which is used in butterfly algorithm

	long barrett_k__ = 0; ///< an integer larger than 2 * log2(prime_) + 1

	word64 barrett_ratio__ = 0; ///< ratio = floor(2^k_ / prime_) for barrett reduction

	word64 primitive_root__ = 0; ///< one of the primitive roots in integer ring

	void computeBarrettParameter__();

	word64 findPrimitiveRoot__();

};

} // namespace basic
} // namespace heaan

#endif
