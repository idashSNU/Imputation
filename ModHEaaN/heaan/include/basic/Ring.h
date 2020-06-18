/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_BASIC_RING_H_
#define HEAAN_BASIC_RING_H_

#include "PrimeRing.h"
#include "BigArith.h"

using namespace std;

namespace heaan
{
namespace basic
{

class Ring
{

public:
	Ring() = default;
	Ring(long number_of_primes, long log_degree);

	void mult(const LargePoly &op1, const LargePoly &op2, const long logq, LargePoly &res) const;

	void mult(const LargePoly &op11, const LargePoly &op12, const LargePoly &op2, const long logq,
			  LargePoly &res1, LargePoly &res2) const;
	void normalize(const LargePoly &op, const long logq, vector<int32_t> &res) const;

	void add(const LargePoly &op1, const LargePoly &op2, const long logq, LargePoly &res) const;
	void addInplace(LargePoly &op1, const LargePoly &op2, const long logq) const;
	void addsmall(const LargePoly &op1, const LargePoly &op2, const long logq, LargePoly &res) const;
	void add(const PrimePolyVec &op1, const PrimePolyVec &op2, const long num_primes, PrimePolyVec &res) const;
	void multMonomial(const LargePoly &op, const long mon_deg, const long logq, LargePoly &res) const;

	void multConst(const LargePoly &op1, const BigInt &op2, const long logq, LargePoly &res) const;
	void multConstsmall(const LargePoly &op1, const BigInt &op2, const long logq, LargePoly &res) const;

private:

	long number_of_primes__ = 0; ///< the number of primes which are used to construct this ring

	long log_degree__ = 0; ///< log2(degree)

	long degree__ = 0; ///< degree (= power of two)

	vector<shared_ptr<PrimeRing>> prime_rings__; ///< shared pointer of PrimeRings

};

} // namespace basic
} // namespace heaan

#endif
