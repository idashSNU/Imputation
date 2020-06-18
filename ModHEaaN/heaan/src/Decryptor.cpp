/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Decryptor.h"

namespace heaan
{

void Decryptor::decrypt(const Ciphertext &ctxt, const SecretKey &secret_key, Plaintext &ptxt) const
{
	// get params and ring
	auto params = context__->params__;
	auto ring = context__->ring__;
	auto primeRing = context__->primeRing__;
	long degree = context__->getDegree();

	// get ctxt parameters
	long quantize_bits = ctxt.getQuantizeBits();
	long modulus_bits = ctxt.getModulusBits();
	long number_of_slots = ctxt.getNumberOfSlots();

	PrimePoly pp1(degree);
	PrimePoly pp2(degree);
	LargePoly lp(degree);

	// set parameters
	ptxt.number_of_slots__ = number_of_slots;
	ptxt.quantize_bits__ = quantize_bits;
	ptxt.mx__.resize(degree);

	primeRing->toNTT(ctxt.getAxData(), pp1); // pp1 = ax_ntt
	primeRing->toNTT(secret_key.getSxData(), pp2); // pp2 = sx_ntt

	// primeRing->hadamardMult(pp1, pp2, pp1); // pp1 <- as_ntt
	primeRing->hadamardMultInplace(pp1, pp2);

	primeRing->fromNTT(pp1, lp); // lp = as

	// ring->add(ctxt.getBxData(), lp, modulus_bits, lp); // lp = mx
	ring->addInplace(lp, ctxt.getBxData(), modulus_bits);

	for (int i = 0; i < degree; i++) {
		if (lp[i] > (1 << 31)) {
			lp[i] -= (1 << 31);
		}
	}

	ring->normalize(lp, modulus_bits, ptxt.mx__);
}

void Decryptor::decrypt_alloc(const Ciphertext &ctxt, const SecretKey &secret_key, PrimePoly &pp1, PrimePoly &pp2, LargePoly &lp, Plaintext &ptxt) const
{
	// get params and ring
	auto params = context__->params__;
	auto ring = context__->ring__;
	auto primeRing = context__->primeRing__;
	long degree = context__->getDegree();

	// get ctxt parameters
	long quantize_bits = ctxt.getQuantizeBits();
	long modulus_bits = ctxt.getModulusBits();
	long number_of_slots = ctxt.getNumberOfSlots();

	// set parameters
	ptxt.number_of_slots__ = number_of_slots;
	ptxt.quantize_bits__ = quantize_bits;

	primeRing->toNTT(ctxt.getAxData(), pp1); // pp1 = ax_ntt
	primeRing->toNTT(secret_key.getSxData(), pp2); // pp2 = sx_ntt

	primeRing->hadamardMultInplace(pp1, pp2);
	
	primeRing->fromNTT(pp1, lp); // lp = as

	ring->addInplace(lp, ctxt.getBxData(), modulus_bits);

	for (int i = 0; i < degree; i++) {
		if (lp[i] > (1 << 31)) {
			lp[i] -= (1 << 31);
		}
	}

	ring->normalize(lp, modulus_bits, ptxt.mx__);
}

} // namespace heaan
