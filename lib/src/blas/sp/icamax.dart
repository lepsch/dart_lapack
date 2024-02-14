      import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

int icamax(final int N,final Array<Complex> CX,_final int INCX,) {
  final CX, = CX,_.dim();

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double SMAX;
      int     I,IX;

      var result = 0;
      if (N < 1 || INCX <= 0) return result;
      result = 1;
      if (N == 1) return result;
      if (INCX == 1) {

         // code for increment equal to 1

         SMAX = SCABS1(CX(1));
         for (I = 2; I <= N; I++) {
            if (SCABS1(CX(I)) > SMAX) {
               result = I;
               SMAX = SCABS1(CX(I));
            }
         }
      } else {

         // code for increment not equal to 1

         IX = 1;
         SMAX = SCABS1(CX(1));
         IX = IX + INCX;
         for (I = 2; I <= N; I++) {
            if (SCABS1(CX(IX)) > SMAX) {
               result = I;
               SMAX = SCABS1(CX(IX));
            }
            IX = IX + INCX;
         }
      }
      return result;
      }
