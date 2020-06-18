/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#include "basic/Ring.h"
#include "basic/BigArith.h"
#include <memory>

namespace heaan
{
  namespace basic
  {

    Ring::Ring(long number_of_primes, long log_degree) : number_of_primes__(number_of_primes)
    {
      log_degree__ = log_degree;
      degree__ = 1 << log_degree__;
    }

    void Ring::mult(const LargePoly &op1, const LargePoly &op2, const long logq, LargePoly &res) const
    {
      if (op1.size() != degree__ || op2.size() != degree__)
      {
        throw invalid_argument("[Ring::mult] input size is not degree!");
      }
      for(int i = 0; i < degree__; i++){
        res[i] = 0;
      }

      for(int i = 0; i < degree__; i++){
        for(int j = 0; j < degree__; j++){
          if(j < i+1)
            res[i] += op1[j] * op2[i-j]; 
          else
            res[i] -= op1[j] * op2[degree__+i-j];
        }
        res[i] <<= 32 - logq;
        res[i] >>= 32 - logq;
      }
    }

    void Ring::mult(const LargePoly &op11, const LargePoly &op12, const LargePoly &op2, const long logq,
        LargePoly &res1, LargePoly &res2) const
    { 
      if (op11.size() != degree__ || op12.size() != degree__ || op2.size() != degree__)
      {
        throw invalid_argument("[Ring::mult] input size is not degree!");
      }
      for(int i = 0; i < degree__; i++){
        res1[i] = 0;
        res2[i] = 0;
      }

      for(int i = 0; i < degree__; i++){
        for(int j = 0; j < i+1; j++){
          res1[i] += op11[j] * op2[i-j]; 
          res2[i] += op12[j] * op2[i-j]; 
          if(i==0){
          }
        }
        for(int j = i+1; j < degree__; j++){
          res1[i] -= op11[j] * op2[degree__+i-j];
          res2[i] -= op12[j] * op2[degree__+i-j];
          if(i==0){
          }
        }
        res1[i] <<= 32 - logq;
        res1[i] >>= 32 - logq;
        res2[i] <<= 32 - logq;
        res2[i] >>= 32 - logq;        
      }
    }

    void Ring::normalize(const LargePoly &op, const long logq, vector<int32_t> &res) const
    {
      BigInt q = BigInt(1) << logq;
      BigInt qh = BigInt(1) << (logq - 1);

      for (int i = 0; i < op.size(); ++i)
      {
        if (op[i] > qh)
          res[i] = op[i] - q;
        else
          res[i] = op[i];
      }
    }

    void Ring::add(const LargePoly &op1, const LargePoly &op2, const long logq, LargePoly &res) const
    {
      if (op1.size() != degree__ || op2.size() != degree__)
      {
        throw invalid_argument("[Ring::Add] input size is not degree!");
      }

      for (int i = 0; i < degree__; ++i)
      {
        __heaan_add_logmod_bigint(op1[i], op2[i], logq, res[i]);
      }
    }

    void Ring::addInplace(LargePoly &op1, const LargePoly &op2, const long logq) const
    {
      if (op1.size() != degree__ || op2.size() != degree__)
      {
        throw invalid_argument("[Ring::Add] input size is not degree!");
      }

      for (int i = 0; i < degree__; ++i)
      {
        __heaan_add_logmod_bigint_inplace(op1[i], op2[i], logq);
      }
    }

    void Ring::addsmall(const LargePoly &op1, const LargePoly &op2, const long logq, LargePoly &res) const
    {
      if (op1.size() != degree__ || op2.size() != degree__)
      {
        throw invalid_argument("[Ring::Add] input size is not degree!");
      }
      uint32_t tmp1, tmp2;
      for (int i = 0; i < op1.size(); ++i)
      {
        tmp1 = op1[i];
        tmp2 = op2[i];
        tmp1 += tmp2;
        tmp1 <<= (32 - logq);
        tmp1 >>= (32 - logq);
        res[i] = tmp1; 
      }
    }

    void Ring::add(const PrimePolyVec &op1, const PrimePolyVec &op2, const long num_primes, PrimePolyVec &res) const
    {
      for (size_t i = 0; i < num_primes; ++i)
      {
        prime_rings__[i]->add(op1[i], op2[i], res[i]);
      }
    }


    void Ring::multMonomial(const LargePoly &op1, const long mon_deg, const long logq, LargePoly &res) const
    {
      long shift = mon_deg % (2 * degree__);
      BigInt mod = BigInt(1) << logq;
      if(shift == 0)
      {
        for(long i = 0; i < degree__; ++i)
        {
          res[i] = op1[i];
        }
      }
      else
      {
        LargePoly tmp(degree__);
        if(shift < degree__)
        {
          for(long i = 0; i < degree__; ++i)
          {
            tmp[i] = op1[i];
          }
        }
        else
        {
          for(long i = 0; i < degree__; ++i)
          {
            tmp[i] = -op1[i];
          }
        }
        shift %= degree__;

        for(int i = 0; i < degree__; i++){
          if(i < shift){
            res[i] = -tmp[degree__ - shift + i];
          }
          else{
            res[i] = tmp[i - shift];
          }
        }
      }
    }

    void Ring::multConst(const LargePoly &op1, const BigInt &op2, const long logq, LargePoly &res) const
    {
      BigInt q = BigInt(1) << logq;

      for (int i = 0; i < op1.size(); ++i)
      {
        __heaan_mult_mod_bigint(op1[i], op2, q, res[i]);
      }
    }

    void Ring::multConstsmall(const LargePoly &op1, const BigInt &op2, const long logq, LargePoly &res) const
    {
      for (int i = 0; i < op1.size(); ++i)
      {
        res[i] = op1[i];
        res[i] *= op2;
        res[i] <<= (32 - logq);
        res[i] >>= (32 - logq);
      }
    }

  } // namespace basic
} // namespace heaan
