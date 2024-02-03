      SUBROUTINE ZLA_WWADDW( N, X, Y, W );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * ), W( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      COMPLEX*16         S;
      int                I;
      // ..
      // .. Executable Statements ..
      for (I = 1; I <= N; I++) { // 10
        S = X(I) + W(I);
        S = (S + S) - S;
        Y(I) = ((X(I) - S) + W(I)) + Y(I);
        X(I) = S;
      } // 10
      RETURN;

      // End of ZLA_WWADDW

      }
