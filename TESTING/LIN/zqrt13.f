      void zqrt13(SCALE, M, N, A, LDA, NORMA, ISEED ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, SCALE;
      double             NORMA;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX*16         A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      double             BIGNUM, SMLNUM;
      // ..
      // .. External Functions ..
      double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL DLAMCH, DZASUM, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARNV, ZLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, SIGN
      // ..
      // .. Local Arrays ..
      double             DUMMY( 1 );
      // ..
      // .. Executable Statements ..

      if (M <= 0 || N <= 0) return;

      // benign matrix

      for (J = 1; J <= N; J++) { // 10
         zlarnv(2, ISEED, M, A( 1, J ) );
         if ( J <= M ) {
            A( J, J ) = A( J, J ) + DCMPLX( SIGN( DZASUM( M, A( 1, J ), 1 ), DBLE( A( J, J ) ) ) );
         }
      } // 10

      // scaled versions

      if ( SCALE != 1 ) {
         NORMA = ZLANGE( 'Max', M, N, A, LDA, DUMMY );
         SMLNUM = DLAMCH( 'Safe minimum' );
         BIGNUM = ONE / SMLNUM;
         SMLNUM = SMLNUM / DLAMCH( 'Epsilon' );
         BIGNUM = ONE / SMLNUM;

         if ( SCALE == 2 ) {

            // matrix scaled up

            zlascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
         } else if ( SCALE == 3 ) {

            // matrix scaled down

            zlascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
         }
      }

      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, DUMMY );
      return;

      // End of ZQRT13

      }
