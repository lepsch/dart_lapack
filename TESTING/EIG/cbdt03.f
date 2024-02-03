      SUBROUTINE CBDT03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDU, LDVT, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), S( * )
      COMPLEX            U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

* ======================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               BNORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SCASUM, SLAMCH
      // EXTERNAL LSAME, ISAMAX, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      RESID = ZERO
      IF( N.LE.0 ) RETURN

      // Compute B - U * S * V' one column at a time.

      BNORM = ZERO
      if ( KD.GE.1 ) {

         // B is bidiagonal.

         if ( LSAME( UPLO, 'U' ) ) {

            // B is upper bidiagonal.

            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= N; I++) { // 10
                  WORK( N+I ) = S( I )*VT( I, J )
               } // 10
               cgemv('No transpose', N, N, -CMPLX( ONE ), U, LDU, WORK( N+1 ), 1, CMPLX( ZERO ), WORK, 1 );
               WORK( J ) = WORK( J ) + D( J )
               if ( J.GT.1 ) {
                  WORK( J-1 ) = WORK( J-1 ) + E( J-1 )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J-1 ) ) )
               } else {
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               }
               RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
            } // 20
         } else {

            // B is lower bidiagonal.

            for (J = 1; J <= N; J++) { // 40
               for (I = 1; I <= N; I++) { // 30
                  WORK( N+I ) = S( I )*VT( I, J )
               } // 30
               cgemv('No transpose', N, N, -CMPLX( ONE ), U, LDU, WORK( N+1 ), 1, CMPLX( ZERO ), WORK, 1 );
               WORK( J ) = WORK( J ) + D( J )
               if ( J.LT.N ) {
                  WORK( J+1 ) = WORK( J+1 ) + E( J )
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J ) ) )
               } else {
                  BNORM = MAX( BNORM, ABS( D( J ) ) )
               }
               RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
            } // 40
         }
      } else {

         // B is diagonal.

         for (J = 1; J <= N; J++) { // 60
            for (I = 1; I <= N; I++) { // 50
               WORK( N+I ) = S( I )*VT( I, J )
            } // 50
            cgemv('No transpose', N, N, -CMPLX( ONE ), U, LDU, WORK( N+1 ), 1, CMPLX( ZERO ), WORK, 1 );
            WORK( J ) = WORK( J ) + D( J )
            RESID = MAX( RESID, SCASUM( N, WORK, 1 ) )
         } // 60
         J = ISAMAX( N, D, 1 )
         BNORM = ABS( D( J ) )
      }

      // Compute norm(B - U * S * V') / ( n * norm(B) * EPS )

      EPS = SLAMCH( 'Precision' )

      if ( BNORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         if ( BNORM.GE.RESID ) {
            RESID = ( RESID / BNORM ) / ( REAL( N )*EPS )
         } else {
            if ( BNORM.LT.ONE ) {
               RESID = ( MIN( RESID, REAL( N )*BNORM ) / BNORM ) / ( REAL( N )*EPS )
            } else {
               RESID = MIN( RESID / BNORM, REAL( N ) ) / ( REAL( N )*EPS )
            }
         }
      }

      RETURN

      // End of CBDT03

      }
