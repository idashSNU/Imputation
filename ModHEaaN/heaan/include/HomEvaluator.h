/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_HOMEVALUATOR_H_
#define HEAAN_HOMEVALUATOR_H_

#include "PublicKeyPack.h"
#include "Encryptor.h"

using namespace std;

namespace heaan
{
  class HomEvaluator
  {

    public:
      HomEvaluator() = default;
      HomEvaluator(Context context);

      void add(const Ciphertext &ctxt1, const Ciphertext &ctxt2, Ciphertext &ctxt_out) const;
      
      void constmultWithoutRescale(const Ciphertext &ctxt1, const double &const2, const long quantize_bits, Ciphertext &ctxt_out) const;

      void monomialmultWithoutRescale(const Ciphertext &ctxt1, const long deg2, const double &const2, const long quantize_bits, Ciphertext &ctxt_out) const;
      
    private:
      shared_ptr<Context> context__; ///< shared pointer of context
  };

} // namespace heaan

#endif
