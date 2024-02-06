      void cqrt13(SCALE, M, N, A, LDA, NORMA, ISEED ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, SCALE;
      double               NORMA;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      Complex            A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      double               BIGNUM, SMLNUM;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SCASUM, SLAMCH;
      // EXTERNAL CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARNV, CLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL, SIGN
      // ..
      // .. Local Arrays ..
      double               DUMMY( 1 );
      // ..
      // .. Executable Statements ..

      if (M <= 0 || N <= 0) return;

      // benign matrix

      for (J = 1; J <= N; J++) { // 10
         clarnv(2, ISEED, M, A( 1, J ) );
         if ( J <= M ) {
            A[J][J] = A( J, J ) + CMPLX( sign( SCASUM( M, A( 1, J ), 1 ), double( A( J, J ) ) ) );
         }
      } // 10

      // scaled versions

      if ( SCALE != 1 ) {
         NORMA = CLANGE( 'Max', M, N, A, LDA, DUMMY );
         SMLNUM = SLAMCH( 'Safe minimum' );
         BIGNUM = ONE / SMLNUM;
         SMLNUM = SMLNUM / SLAMCH( 'Epsilon' );
         BIGNUM = ONE / SMLNUM;

         if ( SCALE == 2 ) {

            // matrix scaled up

            clascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
         } else if ( SCALE == 3 ) {

            // matrix scaled down

            clascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
         }
      }

      NORMA = CLANGE( 'One-norm', M, N, A, LDA, DUMMY );
      return;
      }
