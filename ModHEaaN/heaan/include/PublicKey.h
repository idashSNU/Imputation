/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_PUBLICKEY_H_
#define HEAAN_PUBLICKEY_H_

#include "Context.h"

namespace heaan
{
  class PublicKey
  {

    friend class PublicKeyPack;

    public:

    PublicKey() = default;

    PublicKey(const PublicKey &copy) = default;

    void save(ostream &stream) const;
    void save(const string& path)const
    {
      ofstream fout;
      fout.open(path);
      save(fout);
      fout.close();
    }
    /**
      Load a public key
      */
    void load(istream &stream);
    void load(const string& path)
    {
      ifstream fin;
      fin.open(path);
      load(fin);
      fin.close();
    }

    const long getPublicKeyID() const
    {
      return public_key_id__;
    }

    const PrimePoly &getAxData() const
    {
      return ax_ntt__;
    }

    const PrimePoly &getBxData() const
    {
      return bx_ntt__;
    }

    private:


    long public_key_id__ = 0; ///< id of public key

    PrimePoly ax_ntt__;
    PrimePoly bx_ntt__;

  };

} // namespace heaan

#endif
