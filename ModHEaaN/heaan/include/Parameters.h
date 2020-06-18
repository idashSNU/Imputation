/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_PARAMETERS_H_
#define HEAAN_PARAMETERS_H_

#include "basic/Define.h"

using namespace std;

namespace heaan
{

  class Parameters
  {

    public:

      Parameters() = default;

      Parameters(long degree, long modulus_bits, long quantize_bits, long sp_num)
      {
        modulus_bits__ = modulus_bits;
        quantize_bits__ = quantize_bits;
        sp_num__ = sp_num;
        sp_modulus_bits__ = ceil(modulus_bits__ / (double)sp_num__);
        degree__ = degree;
        h__ = (long) ((2.0 * (double) degree__) / 3.0);
        log_degree__ = log2(degree__);
      }

      long getQuantizeBits() const
      {
        return quantize_bits__;
      }

      long getDegree() const
      {
        return degree__;
      }

      long getLogDegree() const
      {
        return log_degree__;
      }

      long getModulusBits() const
      {
        return modulus_bits__;
      }

      long getSpModulusBits() const
      {
        return sp_modulus_bits__;
      }

      double getSigma() const
      {
        return sigma__;
      }

      long getHemmingWeight() const
      {
        return h__;
      }

      void save(ostream& stream) const;

      void save(const string& path)const
      {
        ofstream fout;
        fout.open(path);
        save(fout);
        fout.close();
      }


      /**
        Load a parameter set
        */
      void load(istream& stream);

      void load(const string& path)
      {
        ifstream fin;
        fin.open(path);
        load(fin);
        fin.close();
      }

    private:

      long quantize_bits__ = 0; ///< standart quantization bits used in ciphertexts

      long log_degree__ = 0; ///< log2 of ring dimension

      long degree__ = 0; ///< ring dimension

      long modulus_bits__ = 0; ///< number of maximum modulus bits used by ciphertexts

      long sp_modulus_bits__ = 0; ///< number of special modulus bits which is used in key-swithcing technique

      long sp_num__ = 0; ///< number of public keys for key-swithcing

      double sigma__ = 3.2; ///< Standard Deviation of Discrete Gaussian Distribution

      long h__ = 0; ///< hamming weight of secret key
  };

} // namespace heaan

#endif
