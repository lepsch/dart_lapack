      SUBROUTINE DQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, SCALE;
      double             NORMA;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      double             BIGNUM, SMLNUM;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE;
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

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      // benign matrix

      for (J = 1; J <= N; J++) { // 10
         dlarnv(2, ISEED, M, A( 1, J ) );
         if ( J.LE.M ) {
            A( J, J ) = A( J, J ) + SIGN( DASUM( M, A( 1, J ), 1 ), A( J, J ) )
         }
      } // 10

      // scaled versions

      if ( SCALE.NE.1 ) {
         NORMA = DLANGE( 'Max', M, N, A, LDA, DUMMY )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
         SMLNUM = SMLNUM / DLAMCH( 'Epsilon' )
         BIGNUM = ONE / SMLNUM

         if ( SCALE.EQ.2 ) {

            // matrix scaled up

            dlascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
         } else if ( SCALE.EQ.3 ) {

            // matrix scaled down

            dlascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
         }
      }

      NORMA = DLANGE( 'One-norm', M, N, A, LDA, DUMMY )
      RETURN

      // End of DQRT13

      }
