      void dqrt13(final int SCALE, final int M, final int N, final Matrix<double> A_, final int LDA, final int NORMA, final int ISEED,) {
  final A = A_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, M, N, SCALE;
      double             NORMA;
      int                ISEED( 4 );
      double             A( LDA, * );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
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
      double             DUMMY( 1 );

      if (M <= 0 || N <= 0) return;

      // benign matrix

      for (J = 1; J <= N; J++) { // 10
         dlarnv(2, ISEED, M, A( 1, J ) );
         if ( J <= M ) {
            A[J][J] = A( J, J ) + sign( dasum( M, A( 1, J ), 1 ), A( J, J ) );
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
      }
