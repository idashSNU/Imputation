/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Context.h"

namespace heaan
{

  Context::Context(Parameters params)
  {
    params__ = make_shared<Parameters>(params);

    genRing__();
    genPrimeRing__();

    // Compute rotation group
    rot_group__.resize(params__->getDegree() / 2);
    long fivePows = 1;
    long M = params__->getDegree() * 2;
    for (long i = 0; i < params__->getDegree() / 2; ++i)
    {
      rot_group__[i] = fivePows;
      fivePows *= 5;
      fivePows %= M;
    }

    // Compute powers of ksi which is M-th root of unity in complex field
    ksi_pows__.resize(M + 1);
    for (long j = 0; j < M; ++j)
    {
      double angle = 2.0 * M_PI * j / M;
      ksi_pows__[j].real(cos(angle));
      ksi_pows__[j].imag(sin(angle));
    }
    ksi_pows__[M] = ksi_pows__[0];

  }

  void Context::encode(const Message &msg, const long quantize_bits, Plaintext &ptxt) const
  {

    // Parameters
    long degree = params__->getDegree();
    long number_of_slots = msg.size();

    // Embedding
    Message tmp(number_of_slots);
    fromEMB__(msg, tmp);

 
    // Make a polynomial with integer coefficient
    long i, jdx, idx;
    long gap = degree / 2 / number_of_slots;
    ptxt.mx__.resize(degree);
    for (i = 0, jdx = degree / 2, idx = 0; i < number_of_slots; ++i, jdx += gap, idx += gap)
    {
      __heaan_scale_up(tmp[i].real(), quantize_bits, ptxt.mx__[idx]);
      __heaan_scale_up(tmp[i].imag(), quantize_bits, ptxt.mx__[jdx]);
    }

    // Set plaintext parameters
    ptxt.number_of_slots__ = number_of_slots;
    ptxt.quantize_bits__ = quantize_bits;

  }

  void Context::decode(const Plaintext &ptxt, Message &msg) const
  {

    // Parameters
    long degree = params__->getDegree();
    long number_of_slots = ptxt.number_of_slots__;
    long quantize_bits = ptxt.quantize_bits__;

    // Make a vector of complex from a polynomial with integer coeffcients
    long i, jdx, idx;
    long gap = degree / 2 / number_of_slots;
    double tmp;
    
    for (i = 0, jdx = degree / 2, idx = 0; i < number_of_slots; ++i, jdx += gap, idx += gap)
    {
      __heaan_scale_down(ptxt.mx__[idx], quantize_bits, tmp); msg[i].real(tmp);
      __heaan_scale_down(ptxt.mx__[jdx], quantize_bits, tmp); msg[i].imag(tmp);
    }

    toEMBInplace__(msg);
  }

  void Context::genRing__()
  {
    // The bit-size of primes in CRT modulus
    // double bit_size_of_prime = log2(HEAAN_SP_BOUND) - 1.;

    // Compute the maximum bit-size that we need
    long max_modulus_bit_size = max(2 * params__->getSpModulusBits() + params__->getModulusBits(), params__->getModulusBits() + params__->getModulusBits());
    // long number_of_primes = ceil((max_modulus_bit_size + params__->getLogDegree() + 4) / bit_size_of_prime) + 1; //assign one more  prime so that it is sufficient for rotate and sum 
    long number_of_primes = 1;

    // Set a shared pointer of ring
    ring__ = make_shared<Ring>(number_of_primes, params__->getLogDegree());
    
  }

  void Context::genPrimeRing__()
  {
    word64 prime = 1099511592961; // 2^40 - 17 * 2^11 + 1
    primeRing__ = make_shared<PrimeRing>(prime, params__->getDegree(), params__->getModulusBits());
  }
  

  void Context::toEMB__(const Message &msg, Message &res) const
  {
    res = msg;
    __heaan_bit_reverse_array(res);
    int size = res.size();
    int M = params__->getDegree() << 1;
    for (long len = 2; len <= size; len <<= 1) // fft-style optimize
    {
      for (long i = 0; i < size; i += len)
      {
        long lenh = len >> 1;
        long lenq = len << 2;
        long gap = M / lenq;
        for (long j = 0; j < lenh; ++j)
        {
          long idx = ((rot_group__[j] % lenq)) * gap;
          auto u = res[i + j];
          auto v = res[i + j + lenh];
          v *= ksi_pows__[idx];
          res[i + j] = u + v;
          res[i + j + lenh] = u - v;
        }
      }
    }
  }

void Context::toEMBInplace__(Message &msg) const
  {
    __heaan_bit_reverse_array(msg);
    int size = msg.size();
    int M = params__->getDegree() << 1;
    for (long len = 2; len <= size; len <<= 1) // fft-style optimize
    {
      for (long i = 0; i < size; i += len)
      {
        long lenh = len >> 1;
        long lenq = len << 2;
        long gap = M / lenq;
        for (long j = 0; j < lenh; ++j)
        {
          long idx = ((rot_group__[j] % lenq)) * gap;
          auto u = msg[i + j];
          auto v = msg[i + j + lenh];
          v *= ksi_pows__[idx];
          msg[i + j] = u + v;
          msg[i + j + lenh] = u - v;
        }
      }
    }
  }

  void Context::fromEMB__(const Message &msg, Message &res) const
  {
    int size = msg.size();
    if (msg != res)
    {
      res.resize(size);
      res = msg;
    }
    int M = params__->getDegree() << 1;
    for (long len = size; len >= 1; len >>= 1) // fft-style optimize
    {
      for (long i = 0; i < size; i += len)
      {
        long lenh = len >> 1;
        long lenq = len << 2;
        long gap = M / lenq;
        for (long j = 0; j < lenh; ++j)
        {
          long idx = (lenq - (rot_group__[j] % lenq)) * gap;
          auto u = res[i + j] + res[i + j + lenh];
          auto v = res[i + j] - res[i + j + lenh];
          v *= ksi_pows__[idx];
          res[i + j] = u;
          res[i + j + lenh] = v;
        }
      }
    }
    __heaan_bit_reverse_array(res);
    for (long i = 0; i < size; ++i)
    {
      res[i] /= size;
    }
  }


  void Context::subFromGauss__(const LargePoly &op, const long logq, LargePoly &res) const
  {
    BigInt q = BigInt(1) << logq;
    double sigma = params__->getSigma();
    for (int i = 0; i < op.size(); i += 2)
    {
      double r1 = (double) rand() / RAND_MAX;
      double r2 = (double) rand() / RAND_MAX;
      
      double theta = 2 * M_PI * r1;
      double rr = sqrt(-2.0 * log(r2)) * sigma;
      BigInt tmp1 = -op[i];
      BigInt tmp2 = -op[i + 1];
      long e1 = floor(rr * cos(theta) + 0.5);
      long e2 = floor(rr * sin(theta) + 0.5);
      BigInt c1 = e1;
      BigInt c2 = e2;
      __heaan_add_logmod_bigint(tmp1, c1, logq, res[i]);
      __heaan_add_logmod_bigint(tmp2, c2, logq, res[i + 1]);
    }
  }

  

  void Context::addGauss__(const LargePoly &op, const long logq, LargePoly &res) const
  {
    BigInt q = BigInt(1) << logq;
    double sigma = params__->getSigma();

    unsigned int seed = (omp_get_thread_num() + 1) ^ (unsigned int)time(NULL);
    
    for (int i = 0; i < op.size(); i += 2)
    {
      double r1, r2;
      r1 = (double)rand_r(&seed) / RAND_MAX;
      r2 = (double)rand_r(&seed) / RAND_MAX;

      double theta = 2 * M_PI * r1;
      double rr = sqrt(-2.0 * log(r2)) * sigma;
      long e1 = floor(rr * cos(theta) + 0.5);
      long e2 = floor(rr * sin(theta) + 0.5);
      BigInt c1 = e1;
      BigInt c2 = e2;

      __heaan_add_logmod_bigint(op[i], c1, logq, res[i]);
      __heaan_add_logmod_bigint(op[i + 1], c2, logq, res[i + 1]);
    }
  }

  void Context::sampleHWT__(LargePoly &res) const
  {
    long degree = params__->getDegree();
    res.resize(degree);
    int idx = 0;
    int hwt = params__->getHemmingWeight();
    while (idx < hwt)
    {
      int i = rand() % degree;
      if (res[i] == 0)
      {
        res[i] = (rand() % 2) ? 1 : -1;
        idx++;
      }
    }
  }


  void Context::sampleZO__(LargePoly &res) const
  {
    unsigned int seed;

    long degree = params__->getDegree();
    res.resize(degree);
    seed = (omp_get_thread_num() + 1) ^ (unsigned int)time(NULL);
    for (int i = 0; i < degree; ++i)
    {
      res[i] = (rand_r(&seed) % 2) ? 0 : ((rand_r(&seed)>>1) % 2) ? 1 : -1;
    }
  }


  void Context::sampleUniform__(const long modulus_bits, LargePoly &res) const
  {
    long degree = params__->getDegree();
    res.resize(degree);
    for (size_t i = 0; i < degree; ++i)
    {
      __heaan_random_bits(res[i], modulus_bits);
    }
  }

  void Context::sampleUniform__(const long modulus_bits, PrimePoly &res) const
  {
    long degree = params__->getDegree();
    res.resize(degree);
    for (size_t i = 0; i < degree; ++i)
    {
      __heaan_random_bits(res[i], primeRing__->getPrime());
    }
  }

} // namespace heaan
