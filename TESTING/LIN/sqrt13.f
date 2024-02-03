      SUBROUTINE SQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, SCALE;
      REAL               NORMA
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      REAL               A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      REAL               BIGNUM, SMLNUM
      // ..
      // .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE
      // EXTERNAL SASUM, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARNV, SLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SIGN
      // ..
      // .. Local Arrays ..
      REAL               DUMMY( 1 )
      // ..
      // .. Executable Statements ..

      if (M <= 0 || N <= 0) RETURN;

      // benign matrix

      for (J = 1; J <= N; J++) { // 10
         slarnv(2, ISEED, M, A( 1, J ) );
         if ( J <= M ) {
            A( J, J ) = A( J, J ) + SIGN( SASUM( M, A( 1, J ), 1 ), A( J, J ) )
         }
      } // 10

      // scaled versions

      if ( SCALE != 1 ) {
         NORMA = SLANGE( 'Max', M, N, A, LDA, DUMMY )
         SMLNUM = SLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
         SMLNUM = SMLNUM / SLAMCH( 'Epsilon' )
         BIGNUM = ONE / SMLNUM

         if ( SCALE == 2 ) {

            // matrix scaled up

            slascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
         } else if ( SCALE == 3 ) {

            // matrix scaled down

            slascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
         }
      }

      NORMA = SLANGE( 'One-norm', M, N, A, LDA, DUMMY )
      RETURN

      // End of SQRT13

      }
