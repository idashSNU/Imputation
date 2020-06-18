/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Encryptor.h"

namespace heaan
{

void Encryptor::encrypt(const Plaintext &ptxt, const long modulus_bits, const PublicKey &enc_key, Ciphertext &ctxt) const
{
	auto params = context__->params__;
	auto ring = context__->ring__;
	auto primeRing = context__->primeRing__;
	long logq = modulus_bits;

	ctxt.quantize_bits__ = ptxt.quantize_bits__;
	ctxt.modulus_bits__ = modulus_bits;
	ctxt.number_of_slots__ = ptxt.number_of_slots__;

	long degree = params->getDegree();
	ctxt.allocate(degree);

	LargePoly vx(degree);

	context__->sampleZO__(vx);

	PrimePoly vx_ntt(degree);
	primeRing->toNTT(vx, vx_ntt);

	PrimePoly ctxt_ax_ntt(degree);
	primeRing->hadamardMult(enc_key.getAxData(), vx_ntt, ctxt_ax_ntt);
	
	PrimePoly ctxt_bx_ntt(degree);
	primeRing->hadamardMult(enc_key.getBxData(), vx_ntt, ctxt_bx_ntt);
	
	primeRing->fromNTT(ctxt_ax_ntt, ctxt.ax__);
	primeRing->fromNTT(ctxt_bx_ntt, ctxt.bx__);

	context__->addGauss__(ctxt.ax__, logq, ctxt.ax__);
	context__->addGauss__(ctxt.bx__, logq, ctxt.bx__);
	LargePoly unsigned_mx;
	unsigned_mx.resize(params->getDegree());

	for(int i = 0; i < degree; i++){
		unsigned_mx[i] = ptxt.mx__[i];
	}
	context__->ring__->add(ctxt.bx__, unsigned_mx, logq, ctxt.bx__);
}

} // namespace heaan
