/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#include "PublicKey.h"

namespace heaan
{

void PublicKey::save(ostream &stream) const
{
	long degree = ax_ntt__.size();
	stream.write(reinterpret_cast<const char *>(&public_key_id__), sizeof(long));
	stream.write(reinterpret_cast<const char *>(&degree), sizeof(long));
	// save ax part
	stream.write(reinterpret_cast<const char *>(ax_ntt__.data()), sizeof(word64) * degree);

	// save bx part
	stream.write(reinterpret_cast<const char *>(bx_ntt__.data()), sizeof(word64) * degree);

}

void PublicKey::load(istream &stream)
{
	long degree;
	stream.read(reinterpret_cast<char *>(&public_key_id__), sizeof(long));
	stream.read(reinterpret_cast<char *>(&degree), sizeof(long));
	
	// read ax part
	ax_ntt__.resize(degree);
	stream.read(reinterpret_cast<char *>(ax_ntt__.data()), sizeof(word64) * degree);

	// read bx part
	bx_ntt__.resize(degree);
	stream.read(reinterpret_cast<char *>(bx_ntt__.data()), sizeof(word64) * degree);

}

} // namespace heaan