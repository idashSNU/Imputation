/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_CONTEXT_H_
#define HEAAN_CONTEXT_H_

#include "Parameters.h"
#include "Plaintext.h"
#include "Ciphertext.h"
#include <omp.h>

using namespace std;
using namespace heaan::basic;

namespace heaan
{
  class Context
  {


    friend class SecretKey;

    friend class PublicKeyPack;

    friend class Encryptor;

    friend class Decryptor;

    friend class HomEvaluator;

    public:

    Context() = default;

    Context(Parameters params);

    void encode(const Message &msg, const long quantize_bits, Plaintext &ptxt) const;

    void decode(const Plaintext &ptxt, Message &msg) const;

    inline long getDegree() const
    {
      return params__->getDegree();
    }

    // private:
    shared_ptr<Ring> ring__; ///< pointer to the underlying ring
    shared_ptr<PrimeRing> primeRing__; ///< pointer to the underlying ring

    shared_ptr<Parameters> params__; ///< pointer to the underlying parameters

    vector<long> rot_group__; ///< powers of 5(generator in Z_M^*) roots

    vector<complex<double> > ksi_pows__; ///< powers of M-th roots of unity

    void genRing__();
    void genPrimeRing__();
    void toEMB__(const Message &msg, Message &res) const;
    void toEMBInplace__(Message &msg) const;

    void fromEMB__(const Message &msg, Message &res) const;

    void subFromGauss__(const LargePoly &op, const long modulus_bits, LargePoly &res) const;
    void addGauss__(const LargePoly &op, const long modulus_bits, LargePoly &res) const;
    void sampleHWT__(LargePoly &res) const;
    void sampleZO__(LargePoly &res) const;
    void sampleUniform__(const long modulus_bits, LargePoly &res) const;
    void sampleUniform__(const long modulus_bits, PrimePoly &res) const;

  };

} // namespace heaan

 #endif
