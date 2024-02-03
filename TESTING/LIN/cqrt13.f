      SUBROUTINE CQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, SCALE;
      REAL               NORMA
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX            A( LDA, * )
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
      REAL               CLANGE, SCASUM, SLAMCH
      // EXTERNAL CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARNV, CLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL, SIGN
      // ..
      // .. Local Arrays ..
      REAL               DUMMY( 1 )
      // ..
      // .. Executable Statements ..

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      // benign matrix

      DO 10 J = 1, N
         clarnv(2, ISEED, M, A( 1, J ) );
         if ( J.LE.M ) {
            A( J, J ) = A( J, J ) + CMPLX( SIGN( SCASUM( M, A( 1, J ), 1 ), REAL( A( J, J ) ) ) )
         }
   10 CONTINUE

      // scaled versions

      if ( SCALE.NE.1 ) {
         NORMA = CLANGE( 'Max', M, N, A, LDA, DUMMY )
         SMLNUM = SLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
         SMLNUM = SMLNUM / SLAMCH( 'Epsilon' )
         BIGNUM = ONE / SMLNUM

         if ( SCALE.EQ.2 ) {

            // matrix scaled up

            clascl('General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO );
         } else if ( SCALE.EQ.3 ) {

            // matrix scaled down

            clascl('General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO );
         }
      }

      NORMA = CLANGE( 'One-norm', M, N, A, LDA, DUMMY )
      RETURN

      // End of CQRT13

      }
