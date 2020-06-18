/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_BASIC_BIGARITH_H_
#define HEAAN_BASIC_BIGARITH_H_

#include "SmallArith.h"

typedef uint32_t BigInt;

typedef long double BigReal;

typedef vector<BigInt> LargePoly;

inline static void __heaan_random_bits(BigInt &res, int bits)
{
  res = rand();
  res <<= 32 - bits;
  res >>= 32 - bits;
}

inline static void __heaan_random_bits(word64 &res, int bits)
{
  res = rand();
  res <<= 64 - bits;
  res >>= 64 - bits;
}

inline static void __heaan_negate_mod_bigint(const BigInt &op, const BigInt &mod, BigInt &res)
{
  if(op == 0) {
    res = 0;
  } else {
    res = mod - op;
  }
}

inline static void __heaan_mod_bigint(const BigInt &op, const BigInt &mod, BigInt &res)
{
  res = op % mod;
}

inline static void __heaan_logmod_bigint(const BigInt &op, const long &logmod, BigInt &res)
{
  res = op << (32 - logmod);
  res >>= 32- logmod;
}

inline static void __heaan_add_mod_bigint(const BigInt &op1, const BigInt &op2, const BigInt &mod, BigInt &res)
{
  res = op1 + op2;
  res %= mod;
}
inline static void __heaan_add_logmod_bigint(const BigInt &op1, const BigInt &op2, const long &logmod, BigInt &res)
{
  res = op1 + op2;
  res <<= (32 - logmod);
  res >>= 32- logmod;
}

inline static void __heaan_add_logmod_bigint_inplace(BigInt &op1, const BigInt &op2, const long &logmod)
{
  op1 += op2;
  op1 <<= (32 - logmod);
  op1 >>= 32- logmod;
}

inline static void __heaan_add_mod_bigint_smallint(const BigInt &op1, const long &op2, const BigInt &mod, BigInt &res)
{
  res = op1 + op2;
  res %= mod;
}

inline static void __heaan_sub_mod_bigint(const BigInt &op1, const BigInt &op2, const BigInt &mod, BigInt &res)
{
  res = op1 - op2;
  res %= mod;
}

inline static void __heaan_sub_logmod_bigint(const BigInt &op1, const BigInt &op2, const long &logmod, BigInt &res)
{
  res = op1 - op2;
  res <<= (32 - logmod);
  res >>= 32- logmod;
}


inline static void __heaan_sub_mod_bigint_smallint(const BigInt &op1, const long &op2, const BigInt &mod, BigInt &res)
{
  res = op1 - op2;
  res %= mod;
}


inline static void __heaan_mult_mod_bigint(const BigInt &op1, const BigInt &op2, const BigInt &mod, BigInt &res)
{
  word64 tmp = (word64)op1 * op2;
  res = tmp % mod;
}

inline static void __heaan_mult_mod_bigint_smallint(const BigInt &op1, const BigInt &op2, const BigInt &mod, BigInt &res)
{
  word64 tmp = (word64)op1 * op2;
  res = tmp % mod;
}

inline static void __heaan_left_shift_mod_bigint(const BigInt &op, const long shift, const BigInt &mod, BigInt &res)
{
  word64 tmp = static_cast<word64>(op) << shift;
  tmp = tmp % static_cast<word64>(mod);
  res = (BigInt) tmp;
}

inline static void __heaan_right_shift_bigint(const BigInt &op, const long shift, BigInt &res)
{
  res = op >> shift;
}

inline static void __heaan_scale_up(const double &op, const int &shift, int32_t &res)
{
  res = (int32_t)(op * pow(2, shift));
}

inline static void __heaan_scale_down(const int32_t &op, const int &shift, double &res)
{
  res = (double) op / (double) pow(2,shift);
}

/// output maximun bit size of coefficients
inline static long __heaan_max_coeff_bit_size(const LargePoly &fx)
{
  long bit_size[fx.size()];
  long max_bit_size = 0;

    for (long i = 0; i < fx.size(); i++)
    {
      bit_size[i] = (int)log2(fx[i])+1;
    }

    for(long i = 0; i < fx.size(); i++)
      max_bit_size = max(max_bit_size, bit_size[i]);

  return max_bit_size;
}


#endif
