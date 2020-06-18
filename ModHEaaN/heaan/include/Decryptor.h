/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_DECRYPTOR_H_
#define HEAAN_DECRYPTOR_H_

#include "PublicKeyPack.h"

namespace heaan
{
  class Decryptor
  {

    public:
      Decryptor() = default;

      Decryptor(Context context)
      {
        context__ = make_shared<Context>(context);  ///< shared pointer of context
      }

      void decrypt(const Ciphertext &ctxt, const SecretKey &secret_key, Plaintext &ptxt) const;

      void decrypt_alloc(const Ciphertext &ctxt, const SecretKey &secret_key, PrimePoly &pp1, PrimePoly &pp2, LargePoly &lp, Plaintext &ptxt) const;

      void decrypt(const Ciphertext &ctxt, const SecretKey &secret_key, Message &msg) const
      {
        Plaintext ptxt;
        decrypt(ctxt, secret_key, ptxt);
        context__->decode(ptxt, msg);
      }

      void decrypt_alloc(const Ciphertext &ctxt, const SecretKey &secret_key, PrimePoly &pp1, PrimePoly &pp2, LargePoly &lp, Plaintext &ptxt, Message &msg) const
      {
        decrypt_alloc(ctxt, secret_key, pp1, pp2, lp, ptxt);
        context__->decode(ptxt, msg);
      }

      shared_ptr<Context> context__;
  };

} // namespace heaan

#endif /* DECRYPTOR_H_ */
