/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Ciphertext.h"

namespace heaan
{

  Ciphertext::Ciphertext(long degree, long number_of_slots, long modulus_bits, long quantize_bits)
    : number_of_slots__(number_of_slots), modulus_bits__(modulus_bits), quantize_bits__(quantize_bits)
  {
    ax__.resize(degree);
    bx__.resize(degree);
  }

  void Ciphertext::copyParams(const Ciphertext &o)
  {
    number_of_slots__ = o.number_of_slots__;
    quantize_bits__ = o.quantize_bits__;
    modulus_bits__ = o.modulus_bits__;
  }

  void Ciphertext::save(ostream &stream) const
  {
    long degree = ax__.size();
    stream.write(reinterpret_cast<const char *>(&degree), sizeof(long));
    stream.write(reinterpret_cast<const char *>(&number_of_slots__), sizeof(long));
    stream.write(reinterpret_cast<const char *>(&modulus_bits__), sizeof(long));
    stream.write(reinterpret_cast<const char *>(&quantize_bits__), sizeof(long));

    
    stream.write(reinterpret_cast<const char *>(ax__.data()), sizeof(BigInt) * degree);
    stream.write(reinterpret_cast<const char *>(bx__.data()), sizeof(BigInt) * degree);
  }

  void Ciphertext::load(istream &stream)
  {
    long degree;
    stream.read(reinterpret_cast<char *>(&degree), sizeof(long));
    stream.read(reinterpret_cast<char *>(&number_of_slots__), sizeof(long));
    stream.read(reinterpret_cast<char *>(&modulus_bits__), sizeof(long));
    stream.read(reinterpret_cast<char *>(&quantize_bits__), sizeof(long));
    stream.read(reinterpret_cast<char *>(ax__.data()), sizeof(BigInt) * degree);
    stream.read(reinterpret_cast<char *>(bx__.data()), sizeof(BigInt) * degree);
  }

} // namespace heaan
