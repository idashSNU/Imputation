/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_BASIC_PRIMERING_H_
#define HEAAN_BASIC_PRIMERING_H_

#include "PrimeModulus.h"
#include "PrimeRing.h"
#include "BigArith.h"

namespace heaan
{
namespace basic
{

class PrimeRing
{

	friend class Context;

public:
	PrimeRing() = default;

	PrimeRing(const word64 &prime, long degree);

	PrimeRing(const word64 &prime, long degree, long logq);

	PrimeRing(PrimeRing &other) = default;

	auto getPrime() const
	{
		return prime_modulus__->prime__;
	}

	auto getPrimeModulus() const
	{
		return prime_modulus__;
	}

	auto getBarrettRatio() const
	{
		return prime_modulus__->barrett_ratio__;
	}

	auto getBarrettK() const
	{
		return prime_modulus__->barrett_k__;
	}

	long getDegree()
	{
		return degree__;
	}

	void negate(const PrimePoly &op, PrimePoly &res) const;

	void negateInplace(PrimePoly &op) const;

	void negate(const word64 *op, word64 *res) const;

	void negateInplace(word64 *op) const;

	void add(const PrimePoly &op1, const PrimePoly &op2, PrimePoly &res) const;

	void addInplace(PrimePoly &op1, const PrimePoly &op2) const;

	void add(const word64 *op1, const word64 *op2, word64 *res) const;

	void addInplace(word64 *op1, const word64 *op2) const;

	void sub(const PrimePoly &op1, const PrimePoly &op2, PrimePoly &res) const;

	void subInplace(PrimePoly &op1, const PrimePoly &op2) const;

	void sub(const word64 *op1, const word64 *op2, word64 *res) const;

	void subInplace(word64 *op1, const word64 *op2) const;

	void mod(const PrimePoly &op, PrimePoly &res) const;

	void mod(const word64 *op, word64 *res) const;

	void toNTT(const PrimePoly &op, PrimePoly &res) const;
	void toNTT(const LargePoly &op, PrimePoly &res) const;

	void toNTTInplace(PrimePoly &op) const;

	void toNTTLazyInplace(PrimePoly &op) const;

	void toNTT(const word64 *op, word64 *res) const;
	void toNTT(const BigInt *op, word64 *res) const;

	void toNTTInplace(word64 *op) const;

	void toNTTLazyInplace(word64 *op) const;

	void fromNTT(const PrimePoly &op, PrimePoly &res) const;
	void fromNTT(const PrimePoly &op, LargePoly &res) const;

	void fromNTTInplace(PrimePoly &op) const;

	void fromNTTLazyInplace(PrimePoly &op) const;

	void fromNTT(const word64 *op, word64 *res) const;
	void fromNTT(const word64 *op, BigInt *res) const;

	void fromNTTInplace(word64 *op) const;

	void fromNTTLazyInplace(word64 *op) const;

	void hadamardMult(const PrimePoly &op1, const PrimePoly &op2, PrimePoly &res) const;

	void hadamardMultInplace(PrimePoly &op1, const PrimePoly &op2) const;

	void hadamardMult(const word64 *op1, const word64 *op2, word64 *res) const;

	void hadamardMultInplace(word64 *op1, const word64 *op2) const;

	void constMult(const PrimePoly &op1, const word64 op2, PrimePoly &res) const;

	void constMultInplace(PrimePoly &op1, const word64 op2) const;

	void constMult(const word64 *op1, const word64 op2, word64 *res) const;

	void constMultInplace(word64 *op1, const word64 op2) const;

	void frobeniusMap(const PrimePoly &op, long pow, PrimePoly &res) const;

	void frobeniusMap(const word64 *op, long pow, word64 *res) const;

	void multMonomial(const PrimePoly &op, long pow, PrimePoly &res) const;

	void multMonomial(const word64 *op, long pow, word64 *res) const;

	void print(const PrimePoly &op, string str);

private:
	long degree__ = 0;

	long two_times_degree__ = 0;

	long logq__ = 0;

	word64 prime__ = 0;

	shared_ptr<PrimeModulus> prime_modulus__;

	word64 dimension_inverse__ = 0;

	word64 scaled_dimension_inverse__ = 0;

	vector<word64> power_of_roots__; ///< Power of primitive Roots

	vector<word64> inverse_power_of_roots__; ///< Power of primitive inverse roots

	vector<word64> inverse_power_of_roots_div_two__; ///< Power of primitive inverse roots

	vector<word64> scaled_power_of_roots__; ///< Power of primitive Roots (scaled)

	vector<word64> inverse_scaled_power_of_roots__; ///< Power of primitive inverse roots (scaled)

	vector<word64> inverse_scaled_power_of_roots_div_two__; ///< Power of primitive inverse roots

	void computeNTTParameters__();

};

} // namespace basic
} // namespace heaan

#endif
