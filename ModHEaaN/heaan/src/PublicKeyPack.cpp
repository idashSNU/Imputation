/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "PublicKeyPack.h"

namespace heaan
{
	void PublicKeyPack::genEncryptionKey__(SecretKey &secret_key)
	{
		auto params = context__->params__;
		auto ring = context__->ring__;
		auto primeRing = context__->primeRing__;
		long logq = params->getModulusBits();
		long degree = params->getDegree();

		// No NTT version
		LargePoly ax(degree);
		LargePoly bx(degree);

		context__->sampleUniform__(logq, ax);
		ring->mult(ax, secret_key.getSxData(), logq, bx);

		context__->subFromGauss__(bx, logq, bx);

		// save as NTT
		PublicKey encryption_key;
		encryption_key.ax_ntt__.resize(params->getDegree());
		encryption_key.bx_ntt__.resize(params->getDegree());
		primeRing->toNTT(ax, encryption_key.ax_ntt__);
		primeRing->toNTT(bx, encryption_key.bx_ntt__);

		encryption_key.public_key_id__ = 0;
	  
	  encryption_key.save(getKeyDirPath() + "/PK/EncKey.bin");
	}
} // namespace heaan
