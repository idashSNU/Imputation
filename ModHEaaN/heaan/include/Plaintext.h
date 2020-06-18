/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_PLAINTEXT_H_
#define HEAAN_PLAINTEXT_H_

#include "basic/Ring.h"

namespace heaan
{

class Plaintext
{

	friend class Context;

	friend class HomEvaluator;

	friend class Encryptor;

	friend class Decryptor;

public:
	Plaintext() = default;

	long getNumberOfSlots() const
	{
		return number_of_slots__;
	}

	long getQuantizeBits() const
	{
		return quantize_bits__;
	}

	const vector<int32_t> &getMxData() const
	{
		return mx__;
	}

    void save(ostream& stream) const
    {
        long degree = mx__.size();
        stream.write(reinterpret_cast<const char*>(&degree), sizeof(long));
        stream.write(reinterpret_cast<const char*>(&number_of_slots__), sizeof(long));
        stream.write(reinterpret_cast<const char*>(&quantize_bits__), sizeof(long));
        stream.write(reinterpret_cast<const char *>(mx__.data()), sizeof(int32_t) * degree);
    }
 
    void load(istream& stream)
    {
        long degree;
        stream.read(reinterpret_cast<char*>(&degree), sizeof(long));
        stream.read(reinterpret_cast<char*>(&number_of_slots__), sizeof(long));
        stream.read(reinterpret_cast<char*>(&quantize_bits__), sizeof(long));
        mx__.resize(degree);
        stream.read(reinterpret_cast<char *>(mx__.data()), sizeof(int32_t) * degree);
    }
	
    long number_of_slots__ = 0; ///< number of slots plaintext encodes

	long quantize_bits__ = 0; ///< number of plaintext quantization bits

	vector<int32_t> mx__; ///< ring element

	PrimePolyVec mx_rns__; ///< ring element in RNS form

	bool in_rns__ = false; ///< is ring element in RNS form exist
};

} // namespace heaan

#endif
