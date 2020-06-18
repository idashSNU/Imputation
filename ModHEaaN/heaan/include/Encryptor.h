/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAAN_ENCRYPTOR_H_
#define HEAAN_ENCRYPTOR_H_

#include "PublicKeyPack.h"

namespace heaan
{

class Encryptor
{

public:
	Encryptor() = default;
	Encryptor(Context context)
	{
		context__ = make_shared<Context>(context);
	}

	void encrypt(const Message &msg, const PublicKey &enc_key, Ciphertext &ctxt) const
	{
		long quantize_bits = context__->params__->getQuantizeBits();
		Plaintext ptxt;
		context__->encode(msg, quantize_bits, ptxt);
		encrypt(ptxt, context__->params__->getModulusBits(), enc_key, ctxt);
	}

	void encrypt(const Plaintext &ptxt, const long modulus_bits, const PublicKey &enc_key, Ciphertext &ctxt) const;

private:
	shared_ptr<Context> context__{nullptr}; ///< shared pointer of context
};

} // namespace heaan

#endif
