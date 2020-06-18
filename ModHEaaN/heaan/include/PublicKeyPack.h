/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_KEYGENERATOR_H_
#define HEAAN_KEYGENERATOR_H_

#include "SecretKey.h"
#include "PublicKey.h"
#include <sys/stat.h>
namespace heaan
{
  class PublicKeyPack
  {

    public:
      PublicKeyPack() = default;

      PublicKeyPack(Context context)
      {
        context__ = make_shared<Context>(context);

      }

      PublicKeyPack(Context context, SecretKey &secret_key, const string& key_dir_path)
      {
        setKeyDirPath(key_dir_path);
        context__ = make_shared<Context>(context);
        genEncryptionKey__(secret_key);
      }

      void load(const string& key_dir_path)
      {
        setKeyDirPath(key_dir_path);
      }

      const PublicKey getEncKey() const
      {
        PublicKey enc_key;
        enc_key.load(getKeyDirPath() + "/PK/EncKey.bin");
        return enc_key;
      }


      void setKeyDirPath(const string& key_dir_path)
      {
        mkdir((key_dir_path + "/PK/").data(), 0775);
        this->key_dir_path = key_dir_path;
      }

      string getKeyDirPath() const
      {
        return this->key_dir_path;
      }



    private:
      void genEncryptionKey__(SecretKey &secret_key);


      shared_ptr<Context> context__ = nullptr; ///< shared pointer of context
  
      string key_dir_path;

  };

} // namespace heaan

#endif
