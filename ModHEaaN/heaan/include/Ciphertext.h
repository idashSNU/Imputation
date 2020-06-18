/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_CIPHERTEXT_H_
#define HEAAN_CIPHERTEXT_H_

#include "basic/Ring.h"

namespace heaan
{
  class Ciphertext
  {

    friend class HomEvaluator;

    friend class Encryptor;

    public:
    Ciphertext() = default;

    Ciphertext(const Ciphertext &copy) = default;

    Ciphertext(long degree, long number_of_slots, long modulus_bits, long quantize_bits);

    void copyParams(const Ciphertext &o);

    void allocate(const long degree) {
      if (ax__.size() != degree) ax__.resize(degree);
      if (bx__.size() != degree) bx__.resize(degree);
    }

    const long getNumberOfSlots() const
    {
      return number_of_slots__;
    }

    const long getModulusBits() const
    {
      return modulus_bits__;
    }

    const long getQuantizeBits() const
    {
      return quantize_bits__;
    }

    const LargePoly &getAxData() const
    {
      return ax__;
    }

    const LargePoly &getBxData() const
    {
      return bx__;
    }

    void save(ostream &stream) const;

    void load(istream &stream);


    void save(const string& path)const
    {
      ofstream fout;
      fout.open(path);
      save(fout);
      fout.close();
    }

    void load(const string& path)
    {
      ifstream fin;
      fin.open(path);
      load(fin);
      fin.close();
    }

    private:

    long number_of_slots__ = 0; ///< number of slots ciphertext encrypts

    long modulus_bits__ = 0; ///< number of ciphertext modulus bits

    long quantize_bits__ = 0; ///< number of ciphertext quantization bits

    LargePoly ax__; ///< ring element, part of ciphertext

    LargePoly bx__; ///< ring element, part of ciphertext
  };

} // namespace heaan

#endif
