/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#include "basic/PrimeRing.h"

namespace heaan
{
namespace basic
{

PrimeRing::PrimeRing(const word64 &prime, long degree) : prime__(prime), degree__(degree)
{
	prime_modulus__ = make_shared<PrimeModulus>(prime);
	computeNTTParameters__();
}

PrimeRing::PrimeRing(const word64 &prime, long degree, long logq) : prime__(prime), degree__(degree), logq__(logq)
{
	prime_modulus__ = make_shared<PrimeModulus>(prime);
	computeNTTParameters__();
}

void PrimeRing::negate(const PrimePoly &op, PrimePoly &res) const
{
	negate(op.data(), res.data());
}

void PrimeRing::negateInplace(PrimePoly &op) const
{
	negateInplace(op.data());
}

void PrimeRing::negate(const word64 *op, word64 *res) const
{
	size_t n = degree__;
	for (; n--; op++, res++)
	{
		if (*op != 0)
		{
			*res = prime__ - *op;
		}
		else
		{
			*res = static_cast<word64>(0);
		}
	}
}

void PrimeRing::negateInplace(word64 *op) const
{
	size_t n = degree__;
	for (; n--; op++)
	{
		if (*op != 0)
		{
			*op = prime__ - *op;
		}
	}
}

void PrimeRing::add(const PrimePoly &op1, const PrimePoly &op2, PrimePoly &res) const
{
	add(op1.data(), op2.data(), res.data());
}

void PrimeRing::addInplace(PrimePoly &op1, const PrimePoly &op2) const
{
	addInplace(op1.data(), op2.data());
}

void PrimeRing::add(const word64 *op1, const word64 *op2, word64 *res) const
{
	size_t n = degree__;
	for (; n--; op1++, op2++, res++)
	{
		*res = *op1 + *op2;
		if (((prime__ - *res) >> 63))
		{
			*res -= prime__;
		}
	}
}

void PrimeRing::addInplace(word64 *op1, const word64 *op2) const
{
	size_t n = degree__;
	for (; n--; op1++, op2++)
	{
		*op1 += *op2;
		if (((prime__ - *op1) >> 63))
		{
			*op1 -= prime__;
		}
	}
}

void PrimeRing::sub(const PrimePoly &op1, const PrimePoly &op2, PrimePoly &res) const
{
	sub(op1.data(), op2.data(), res.data());
}

void PrimeRing::subInplace(PrimePoly &op1, const PrimePoly &op2) const
{
	subInplace(op1.data(), op2.data());
}

void PrimeRing::sub(const word64 *op1, const word64 *op2, word64 *res) const
{
	size_t n = degree__;
	for (; n--; op1++, op2++, res++)
	{
		if (*op1 >= *op2)
		{
			*res = *op1 - *op2;
		}
		else
		{
			*res = prime__ - (*op2 - *op1);
		}
	}
}

void PrimeRing::subInplace(word64 *op1, const word64 *op2) const
{
	size_t n = degree__;
	for (; n--; op1++, op2++)
	{
		if (*op1 >= *op2)
		{
			*op1 -= *op2;
		}
		else
		{
			*op1 = prime__ - (*op2 - *op1);
		}
	}
}

void PrimeRing::mod(const PrimePoly &op, PrimePoly &res) const
{
	mod(op.data(), res.data());
}

void PrimeRing::mod(const word64 *op, word64 *res) const
{
	word64 ratio = prime_modulus__->getBarrettRatio();
	long k = prime_modulus__->getBarrettK();
	transform(op, op + degree__, res, [&](auto coeff) {
		word64 _res;
		__heaan_barrett_reduction(coeff, _res, prime__, ratio, k);
		return _res;
	});
}

void PrimeRing::toNTT(const PrimePoly &op, PrimePoly &res) const
{
	// cout << "[PrimeRing::toNTT(PrimePoly. PrimePoly)]" << endl;
	toNTT(op.data(), res.data());
}

void PrimeRing::toNTT(const LargePoly &op, PrimePoly &res) const
{
	// cout << "[PrimeRing::toNTT(LargePoly, PrimePoly)]" << endl;
	toNTT(op.data(), res.data());
}

void PrimeRing::toNTTInplace(PrimePoly &op) const
{
	// cout << "[PrimeRing::toNTTInPlace(PrimePoly)]" << endl;
	toNTTInplace(op.data());
}

void PrimeRing::toNTTLazyInplace(PrimePoly &op) const
{
	toNTTLazyInplace(op.data());
}

void PrimeRing::toNTT(const word64 *op, word64 *res) const
{
	std::copy(op, op + degree__, res);
	toNTTInplace(res);
}

void PrimeRing::toNTT(const BigInt *op, word64 *res) const
{
	// cout << "[PrimeRing::toNTT(BitInt*, word64*)]" << endl;
	std::copy(op, op + degree__, res);

	for(int i = 0; i < degree__; i++){
		if(res[i] > 2147483648){
			res[i] += prime__ - 4294967296;
		}
	}

	// for (int i = 0; i < degree__; i += 1) 
	// 	cout << "[PrimeRing::toNTT(BitInt*, word64*)] copy op[" << i << "] = " << op[i] << " to res[" << i << "] = " << res[i] << endl; 
	toNTTInplace(res);
}

void PrimeRing::toNTTInplace(word64 *op) const
{
	// cout << "[PrimeRing::toNTTInplac(word64*))]" << endl;
	////////////////////////////////////////
	word64 two_times_prime = prime__ * 2;
	size_t n = degree__;
	////////////////////////////////////////
	toNTTLazyInplace(op);
	////////////////////////////////////////
	for (; n--; op++)
	{
		if (*op >= two_times_prime)
		{
			*op -= two_times_prime;
		}
		while (*op >= prime__)
		{
			*op -= prime__;
		}
	}
}

void PrimeRing::toNTTLazyInplace(word64 *op) const
{
	////////////////////////////////////////
	word64 prime = prime__;
	word64 two_times_prime = prime__ * 2;
	size_t n = degree__;
	size_t t = n >> 1;
	////////////////////////////////////////
	// cout << "[PrimeRing::toNTTLazyInPlace] prime = " << prime << endl;

	for (size_t m = 1; m < n; m <<= 1)
	{
		if (t >= 4)
		{
			for (size_t i = 0; i < m; i++)
			{
				size_t j1 = 2 * i * t;
				size_t j2 = j1 + t;
				const word64 W = power_of_roots__[m + i];
				const word64 Wprime = scaled_power_of_roots__[m + i];
				word64 *X = op + j1;
				word64 *Y = X + t;
				word64 currX;
				unsigned long long Q;
				for (size_t j = j1; j < j2; j += 4)
				{
					currX = *X - (two_times_prime & static_cast<word64>(-static_cast<int64_t>(*X >= two_times_prime)));
					__heaan_MUL_64_64_high(Wprime, *Y, &Q);
					Q = *Y * W - Q * prime;
					*X++ = currX + Q;
					*Y++ = currX + (two_times_prime - Q);
					////////////////////////////////////////
					currX = *X - (two_times_prime & static_cast<word64>(-static_cast<int64_t>(*X >= two_times_prime)));
					__heaan_MUL_64_64_high(Wprime, *Y, &Q);
					Q = *Y * W - Q * prime;
					*X++ = currX + Q;
					*Y++ = currX + (two_times_prime - Q);
					////////////////////////////////////////
					currX = *X - (two_times_prime & static_cast<word64>(-static_cast<int64_t>(*X >= two_times_prime)));
					__heaan_MUL_64_64_high(Wprime, *Y, &Q);
					Q = *Y * W - Q * prime;
					*X++ = currX + Q;
					*Y++ = currX + (two_times_prime - Q);
					////////////////////////////////////////
					currX = *X - (two_times_prime & static_cast<word64>(-static_cast<int64_t>(*X >= two_times_prime)));
					__heaan_MUL_64_64_high(Wprime, *Y, &Q);
					Q = *Y * W - Q * prime;
					*X++ = currX + Q;
					*Y++ = currX + (two_times_prime - Q);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < m; i++)
			{
				size_t j1 = 2 * i * t;
				size_t j2 = j1 + t;
				const word64 W = power_of_roots__[m + i];
				const word64 Wprime = scaled_power_of_roots__[m + i];
				word64 *X = op + j1;
				word64 *Y = X + t;
				word64 currX;
				unsigned long long Q;
				for (size_t j = j1; j < j2; j++)
				{
					currX = *X - (two_times_prime & static_cast<word64>(-static_cast<int64_t>(*X >= two_times_prime)));
					__heaan_MUL_64_64_high(Wprime, *Y, &Q);
					Q = W * *Y - Q * prime;
					*X++ = currX + Q;
					*Y++ = currX + (two_times_prime - Q);
				}
			}
		}
		t >>= 1;
	}
}

void PrimeRing::fromNTT(const PrimePoly &op, PrimePoly &res) const
{
	fromNTT(op.data(), res.data());
}

void PrimeRing::fromNTT(const PrimePoly &op, LargePoly &res) const
{
	PrimePoly op_invNTT(degree__);
	fromNTT(op.data(), op_invNTT.data());
	for (int i = 0; i < degree__; i++) {
		// if ( i % 40 == 0 ) cout << "[PrimeRing::fromNTT] mod to p/2 i = " << i << endl;
		if (op_invNTT[i] > prime__ / 2) {
			// cout << "[PrimeRing::fromNTT] res[" << i << "] = " << op_invNTT[i] << " > " << prime__ / 2 << endl;
			op_invNTT[i] -= prime__;
			// cout << "[PrimeRing::fromNTT] moddown to res[" << i << "] = " << op_invNTT[i] << endl;
		} 
		// else {
		// 	// cout << "[PrimeRing::fromNTT] res[" << i << "] = " << op_invNTT[i] << " < " << prime__ / 2 << endl;
		// }
		// cout << "[PrimeRing::fromNTT] res[" << i << "] = " << res[i] << endl;
		// res[i] = op_invNTT[i] & 0xFFFFFFFF;
		res[i] = op_invNTT[i];
		res[i] <<= 32 - logq__;
		res[i] >>= 32 - logq__;
		// cout << "[PrimeRing::fromNTT] mod logq res[" << i << "] = " << res[i] << endl;
	}
	
	
	// for (int i = 0; i < degree__; i++) {

		
	// 	if (res[i] > prime__ / 2) {
	// 		res[i] -= prime__;
	// 	}
	// 	res[i] %= (1 << logq__);

	// 	// // if ( i % 40 == 0 ) cout << "[PrimeRing::fromNTT] mod to p/2 i = " << i << endl;
	// 	// if (res[i] > logq__ / 2) {
	// 	// 	cout << "res[i] = " << res[i] << ", logq = " << logq__ << endl;	 
	// 	// 	res[i] <<= 32 - logq__;
	// 	// 	res[i] >>= 32 - logq__;
	// 	// 	cout << "convert to " << res[i] << endl;
	// 	// }
	// }
}

void PrimeRing::fromNTTInplace(PrimePoly &op) const
{
	fromNTTInplace(op.data());
}

void PrimeRing::fromNTTLazyInplace(PrimePoly &op) const
{
	fromNTTLazyInplace(op.data());
}

void PrimeRing::fromNTT(const word64 *op, word64 *res) const
{
	std::copy(op, op + degree__, res);
	fromNTTInplace(res);
}

void PrimeRing::fromNTT(const word64 *op, BigInt *res) const
{
	// cout << "[PrimeRing::fromNTT(word64*, BigInt*)]" << endl;
	word64* tmp = new word64[degree__];
	std::copy(op, op + degree__, tmp);
	// for (int i = 0; i < degree__; i += 1) cout << "[PrimeRing::fromNTT(BitInt*, word64*)] copy op[" << i << "] = " << op[i] << " to tmp[" << i << "] = " << tmp[i] << endl; 
	fromNTTInplace(tmp);
	for (int i = 0; i < degree__; i++) {
		res[i] = tmp[i] & 0xFFFFFFFF;
		// if (i % 50 == 0) 
		// cout <<  "[PrimeRing::fromNTT(BitInt*, word64*)] res[" << i << "] = " << res[i] << endl;
	}
}

void PrimeRing::fromNTTInplace(word64 *op) const
{
	////////////////////////////////////////
	size_t n = degree__;
	////////////////////////////////////////
	fromNTTLazyInplace(op);
	////////////////////////////////////////
	for (; n--; op++)
	{
		if (*op >= prime__)
		{
			*op %= prime__;
		}
	}
}

void PrimeRing::fromNTTLazyInplace(word64 *op) const
{
	////////////////////////////////////////
	word64 prime = prime__;
	word64 two_times_prime = prime * 2;
	size_t n = degree__;
	size_t t = 1;
	////////////////////////////////////////
	for (size_t m = n; m > 1; m >>= 1)
	{
		size_t j1 = 0;
		size_t h = m >> 1;
		if (t >= 4)
		{
			for (size_t i = 0; i < h; i++)
			{
				size_t j2 = j1 + t;
				const word64 W = inverse_power_of_roots_div_two__[h + i];
				const word64 Wprime = inverse_scaled_power_of_roots_div_two__[h + i];
				word64 *U = op + j1;
				word64 *V = U + t;
				word64 currU;
				word64 T;
				unsigned long long H;
				for (size_t j = j1; j < j2; j += 4)
				{	
					T = two_times_prime - *V + *U;
					currU = *U + *V - (two_times_prime & static_cast<word64>(-static_cast<int64_t>((*U << 1) >= T)));
					*U++ = (currU + (prime & static_cast<word64>(-static_cast<int64_t>(T & 1)))) >> 1;
					__heaan_MUL_64_64_high(Wprime, T, &H);
					*V++ = T * W - H * prime;
					////////////////////////////////////////
					T = two_times_prime - *V + *U;
					currU = *U + *V - (two_times_prime & static_cast<word64>(-static_cast<int64_t>((*U << 1) >= T)));
					*U++ = (currU + (prime & static_cast<word64>(-static_cast<int64_t>(T & 1)))) >> 1;
					__heaan_MUL_64_64_high(Wprime, T, &H);
					*V++ = T * W - H * prime;
					////////////////////////////////////////
					T = two_times_prime - *V + *U;
					currU = *U + *V - (two_times_prime & static_cast<word64>(-static_cast<int64_t>((*U << 1) >= T)));
					*U++ = (currU + (prime & static_cast<word64>(-static_cast<int64_t>(T & 1)))) >> 1;
					__heaan_MUL_64_64_high(Wprime, T, &H);
					*V++ = T * W - H * prime;
					////////////////////////////////////////
					T = two_times_prime - *V + *U;
					currU = *U + *V - (two_times_prime & static_cast<word64>(-static_cast<int64_t>((*U << 1) >= T)));
					*U++ = (currU + (prime & static_cast<word64>(-static_cast<int64_t>(T & 1)))) >> 1;
					__heaan_MUL_64_64_high(Wprime, T, &H);
					*V++ = T * W - H * prime;
				}
				j1 += (t << 1);
			}
		}
		else
		{
			for (size_t i = 0; i < h; i++)
			{
				size_t j2 = j1 + t;
				const word64 W = inverse_power_of_roots_div_two__[h + i];
				const word64 Wprime = inverse_scaled_power_of_roots_div_two__[h + i];
				word64 *U = op + j1;
				word64 *V = U + t;
				word64 currU;
				word64 T;
				unsigned long long H;
				for (size_t j = j1; j < j2; j++)
				{
					T = two_times_prime - *V + *U;
					currU = *U + *V - (two_times_prime & static_cast<word64>(-static_cast<int64_t>((*U << 1) >= T)));
					*U++ = (currU + (prime & static_cast<word64>(-static_cast<int64_t>(T & 1)))) >> 1;
					__heaan_MUL_64_64_high(Wprime, T, &H);
					*V++ = W * T - H * prime;
				}
				j1 += (t << 1);
			}
		}
		t <<= 1;
	}	
}

void PrimeRing::hadamardMult(const PrimePoly &op1, const PrimePoly &op2, PrimePoly &res) const
{
	hadamardMult(op1.data(), op2.data(), res.data());
}

void PrimeRing::hadamardMultInplace(PrimePoly &op1, const PrimePoly &op2) const
{
	hadamardMultInplace(op1.data(), op2.data());
}

void PrimeRing::hadamardMult(const word64 *op1, const word64 *op2, word64 *res) const
{
	size_t n = degree__;
	long ratio = prime_modulus__->getBarrettRatio();
	long k = prime_modulus__->getBarrettK();
	word64 high, low;
	for (; n--; op1++, op2++, res++)
	{
		__heaan_MUL_64_64_128(*op1, *op2, &high, &low);
		// high = 0;
		// low = 0;
		__heaan_barrett_reduction(high, low, *res, prime__, ratio, k);
		// res = 0;
	}
}

void PrimeRing::hadamardMultInplace(word64 *op1, const word64 *op2) const
{
	size_t n = degree__;
	long ratio = prime_modulus__->getBarrettRatio();
	long k = prime_modulus__->getBarrettK();
	word64 high, low;
	for (; n--; op1++, op2++)
	{
		__heaan_MUL_64_64_128(*op1, *op2, &high, &low);
		__heaan_barrett_reduction(high, low, *op1, prime__, ratio, k);
	}
}

void PrimeRing::constMult(const PrimePoly &op1, const word64 op2, PrimePoly &res) const
{
	constMult(op1.data(), op2, res.data());
}

void PrimeRing::constMultInplace(PrimePoly &op1, const word64 op2) const
{
	constMultInplace(op1.data(), op2);
}

void PrimeRing::constMult(const word64 *op1, const word64 op2, word64 *res) const
{
	word64 op2_PsInv;
	prime_modulus__->pseudoInverse(op2, op2_PsInv);
	transform(op1, op1 + degree__, res, [&](auto coeff) {
		word64 _res;
		__heaan_mul_and_reduce_shoup(coeff, op2, op2_PsInv, prime__, _res);
		if (_res >= prime__)
			_res -= prime__;
		return _res;
	});
}

void PrimeRing::constMultInplace(word64 *op1, const word64 op2) const
{
	word64 op2_PsInv;
	prime_modulus__->pseudoInverse(op2, op2_PsInv);
	transform(op1, op1 + degree__, op1, [&](auto coeff) {
		word64 _res;
		__heaan_mul_and_reduce_shoup(coeff, op2, op2_PsInv, prime__, _res);
		if (_res >= prime__)
			_res -= prime__;
		return _res;
	});
}

void PrimeRing::frobeniusMap(const PrimePoly &op, long pow, PrimePoly &res) const
{
	frobeniusMap(op.data(), pow, res.data());
}

void PrimeRing::frobeniusMap(const word64 *op, long pow, word64 *res) const
{
	if (pow % 2 == 0)
	{
		throw invalid_argument("pow should be odd in FrobeniusMap");
	}
	else if (pow == -1)
	{
		*res = *op;
		size_t n = degree__ - 1;
		const word64 *ptr1 = op + degree__ - 1;
		word64 *ptr2 = res + 1;
		for (; n--; ptr1--, ptr2++)
		{
			if (*ptr1 != 0)
			{
				*res = prime__ - *ptr1;
			}
			else
			{
				*res = 0;
			}
		}
	}
	else
	{
		for (size_t i = 0; i < degree__; ++i, op++)
		{
			long ipow = i * pow;
			long shift = ipow % two_times_degree__;
			if (shift < degree__)
			{
				res[shift] = *op;
			}
			else
			{
				if (*op != 0)
				{
					res[shift - degree__] = prime__ - *op;
				}
				else
				{
					res[shift - degree__] = 0;
				}
			}
		}
	}
}

void PrimeRing::multMonomial(const PrimePoly &op, long pow, PrimePoly &res) const
{
	multMonomial(op.data(), pow, res.data());
}

void PrimeRing::multMonomial(const word64 *op, long pow, word64 *res) const
{
	long shift = pow % two_times_degree__;
	if (shift == 0)
	{
		std::copy(op, op + degree__, res);
		return;
	}
	else
	{
		word64 *tmp = new word64[degree__];
		if (shift < degree__)
		{
			std::copy(op, op + degree__, tmp);
		}
		else
		{
			negate(op, tmp);
		}
		shift %= degree__;
		for (size_t i = 0, j = degree__ - shift; i < shift; ++i, ++j)
		{
			if (tmp[j] != 0)
				res[i] = prime__ - tmp[j];
			else
				res[i] = 0;
		}
		for (size_t i = shift; i < degree__; ++i)
		{
			res[i] = tmp[i - shift];
		}
		return;
	}
}

void PrimeRing::computeNTTParameters__()
{
	// cout << "[PrimeRing::computeNTTParameters] degree = " << degree__ << ", prime = " << prime__ << endl; 
	////////////////////////////////////////
	two_times_degree__ = 2 * degree__;
	prime_modulus__->inverse(degree__, dimension_inverse__);
	////////////////////////////////////////
	power_of_roots__.resize(degree__);
	scaled_power_of_roots__.resize(degree__);
	inverse_power_of_roots__.resize(degree__);
	inverse_scaled_power_of_roots__.resize(degree__);
	inverse_power_of_roots_div_two__.resize(degree__);
	inverse_scaled_power_of_roots_div_two__.resize(degree__);
	////////////////////////////////////////
	word64 primitive_root = prime_modulus__->getPrimitiveRoot();
	// cout << "[PrimeRing::computeNTTParameters] primitive root = " << primitive_root << endl;
	word64 root, root_inverse, two_inv;
	prime_modulus__->power(primitive_root, (prime_modulus__->getPrime() - 1) / two_times_degree__, root);
	prime_modulus__->inverse(root, root_inverse);
	power_of_roots__[0] = 1;
	inverse_power_of_roots__[0] = 1;
	for (size_t i = 1; i < degree__; ++i)
	{
		prime_modulus__->mult(power_of_roots__[i - 1], root, power_of_roots__[i]);
		prime_modulus__->mult(inverse_power_of_roots__[i - 1], root_inverse, inverse_power_of_roots__[i]);
	}
	__heaan_bit_reverse_array(power_of_roots__);
	__heaan_bit_reverse_array(inverse_power_of_roots__);
	prime_modulus__->inverse(2, two_inv);
	for (size_t i = 0; i < degree__; ++i)
	{
		prime_modulus__->pseudoInverse(power_of_roots__[i], scaled_power_of_roots__[i]);
		prime_modulus__->pseudoInverse(inverse_power_of_roots__[i], inverse_scaled_power_of_roots__[i]);
		prime_modulus__->mult(two_inv, inverse_power_of_roots__[i], inverse_power_of_roots_div_two__[i]);
		prime_modulus__->pseudoInverse(inverse_power_of_roots_div_two__[i], inverse_scaled_power_of_roots_div_two__[i]);
		// cout << "[PrimeRing::computeNTTParameters] power of roots " << i << " : " << power_of_roots__[i] << ", " << scaled_power_of_roots__[i] << endl;
		// cout << "[PrimeRing::computeNTTParameters] power of roots " << i << " : " << inverse_power_of_roots__[i] << ", " << inverse_scaled_power_of_roots_div_two__[i] << endl;
	}
	prime_modulus__->pseudoInverse(dimension_inverse__, scaled_dimension_inverse__);
}

void PrimeRing::print(const PrimePoly &op, string str)
{
	cout << "modulus = " << prime__ << endl;
	if (degree__ > 100) {
        for (int i = 0; i < degree__; i += 20) {
        cout << str << "[" << i << "] = " << op[i] << endl;
        }
      } else {
        for (int i = 0; i < degree__; i += 1) {
        cout << str << "[" << i << "] = " << op[i] << endl;
        }
      }
}
} // namespace basic
} // namespace heaan
