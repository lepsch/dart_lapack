      SUBROUTINE ZLAIPD( N, A, INDA, VINDA );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INDA, N, VINDA;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, IA, IXA;
      double             BIGNUM;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      BIGNUM = DLAMCH( 'Epsilon' ) / DLAMCH( 'Safe minimum' );
      IA = 1;
      IXA = INDA;
      for (I = 1; I <= N; I++) { // 10
         A( IA ) = DCMPLX( DBLE( A( IA ) ), BIGNUM );
         IA = IA + IXA;
         IXA = IXA + VINDA;
      } // 10
      RETURN;
      }
