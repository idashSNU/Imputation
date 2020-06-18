/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/


#include "HomEvaluator.h"

namespace heaan
{
  HomEvaluator::HomEvaluator(Context context)
  {
    context__ = make_shared<Context>(context);
  }
  
  void HomEvaluator::add(const Ciphertext &ctxt1, const Ciphertext &ctxt2, Ciphertext &ctxt_out) const
  {
    if (ctxt1.quantize_bits__ != ctxt2.quantize_bits__)
    {
      throw invalid_argument("ctxt1 and ctxt2 should have same quantize_bits__ at " + string(__FILE__) + ":" + to_string(__LINE__));
    }
    if (ctxt1.number_of_slots__ != ctxt2.number_of_slots__)
    {
      cerr << "!!!Warning!!!: ctxt1 and ctxt2 have different number_of_slots__" << endl;
    }

    long degree = context__->getDegree();
    ctxt_out.allocate(degree);

    ctxt_out.copyParams(ctxt1);

    context__->ring__->add(ctxt1.ax__, ctxt2.ax__, ctxt_out.modulus_bits__, ctxt_out.ax__);
    context__->ring__->add(ctxt1.bx__, ctxt2.bx__, ctxt_out.modulus_bits__, ctxt_out.bx__);
  }


  void HomEvaluator::constmultWithoutRescale(const Ciphertext &ctxt1, const double &const2, const long quantize_bits, Ciphertext &ctxt_out) const{
    long degree = context__->getDegree();
    long logq = ctxt1.modulus_bits__;
    ctxt_out.allocate(degree);

    int32_t const2_scaleup;
    __heaan_scale_up(const2, quantize_bits, const2_scaleup);

    BigInt const2_scaleup_unsigned = const2_scaleup;;

    ctxt_out.copyParams(ctxt1);
    context__->ring__->multConstsmall(ctxt1.ax__, const2_scaleup, logq, ctxt_out.ax__);
    context__->ring__->multConstsmall(ctxt1.bx__, const2_scaleup, logq, ctxt_out.bx__);
    ctxt_out.quantize_bits__ += quantize_bits;
  }
  
  void HomEvaluator::monomialmultWithoutRescale(const Ciphertext &ctxt1, const long deg2, const double &const2, const long quantize_bits, Ciphertext &ctxt_out) const{
    long degree = context__->getDegree();
    long logq = ctxt1.modulus_bits__;
    ctxt_out.allocate(degree);

    int32_t const2_scaleup;
    __heaan_scale_up(const2, quantize_bits, const2_scaleup);

    BigInt const2_scaleup_unsigned = const2_scaleup;

    ctxt_out.copyParams(ctxt1);
    context__->ring__->multConstsmall(ctxt1.ax__, const2_scaleup, logq, ctxt_out.ax__);
    context__->ring__->multConstsmall(ctxt1.bx__, const2_scaleup, logq, ctxt_out.bx__);
    context__->ring__->multMonomial(ctxt_out.ax__, deg2, logq, ctxt_out.ax__);
    context__->ring__->multMonomial(ctxt_out.bx__, deg2, logq, ctxt_out.bx__);
    ctxt_out.quantize_bits__ += quantize_bits;   

  }

} // namespace heaan
