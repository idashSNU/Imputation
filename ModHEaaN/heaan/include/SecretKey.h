/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#ifndef HEAAN_SECRETKEY_H_
#define HEAAN_SECRETKEY_H_

#include "Context.h"

namespace heaan
{

  class SecretKey
  {

    public:
      SecretKey() = default;

      SecretKey(Context &context)
      {
        context.sampleHWT__(sx__);
      }
	
      const LargePoly &getSxData() const
      {
        return sx__;
      }

      void save(ostream &stream) const
      {
        long degree = sx__.size();
        stream.write(reinterpret_cast<const char *>(&degree), sizeof(long));
        for (long i = 0; i < degree; ++i)
        {
          long tmp = sx__[i];
          stream.write(reinterpret_cast<const char *>(&tmp), sizeof(long));
        }
      }
      void save(const string& path)const
      {
        ofstream fout;
        fout.open(path);
        save(fout);
        fout.close();
      }

      void load(istream &stream)
      {
        long degree;
        stream.read(reinterpret_cast<char *>(&degree), sizeof(long));
        sx__.resize(degree);
        for (long i = 0; i < degree; ++i)
        {
          long tmp;
          stream.read(reinterpret_cast<char *>(&tmp), sizeof(long));
          sx__[i] = tmp;
        }
      }

      void load(const string& path)
      {
        ifstream fin;
        fin.open(path);
        load(fin);
        fin.close();
      }

    private:
      LargePoly sx__; ///< ring element
  };

}; // namespace heaan

#endif
