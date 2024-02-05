      void dqrt13(SCALE, M, N, A, LDA, NORMA, ISEED ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, SCALE;
      double             NORMA;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * );
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
      //- double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARNV, DLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SIGN
      // ..
      // .. Local Arrays ..
      double             DUMMY( 1 );
      // ..
      // .. Executable Statements ..

      if (M <= 0 || N <= 0) return;

      // benign matrix

      for (J = 1; J <= N; J++) { // 10
         dlarnv(2, ISEED, M, A( 1, J ) );
         if ( J <= M ) {
            A[J, J] = A( J, J ) + sign( dasum( M, A( 1, J ), 1 ), A( J, J ) );
         }
      } // 10

      // scaled versions

      if ( SCALE != 1 ) {
         NORMA = dlange( 'Max', M, N, A, LDA, DUMMY );
         SMLNUM = dlamch( 'Safe minimum' );
         BIGNUM = ONE / SMLNUM;
         SMLNUM = SMLNUM / dlamch( 'Epsilon' );
         BIGNUM = ONE / SMLNUM;

         if ( SCALE == 2 ) {

            // matrix scaled up

            dlascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
         } else if ( SCALE == 3 ) {

            // matrix scaled down

            dlascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
         }
      }

      NORMA = dlange( 'One-norm', M, N, A, LDA, DUMMY );
      return;
      }
