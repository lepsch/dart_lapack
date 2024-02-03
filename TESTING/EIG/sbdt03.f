      SUBROUTINE SBDT03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK, RESID );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDU, LDVT, N;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               BNORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SASUM, SLAMCH;
      // EXTERNAL LSAME, ISAMAX, SASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      RESID = ZERO;
      if (N <= 0) return;

      // Compute B - U * S * V' one column at a time.

      BNORM = ZERO;
      if ( KD >= 1 ) {

         // B is bidiagonal.

         if ( LSAME( UPLO, 'U' ) ) {

            // B is upper bidiagonal.

            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= N; I++) { // 10
                  WORK( N+I ) = S( I )*VT( I, J );
               } // 10
               sgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
               WORK( J ) = WORK( J ) + D( J );
               if ( J > 1 ) {
                  WORK( J-1 ) = WORK( J-1 ) + E( J-1 );
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J-1 ) ) );
               } else {
                  BNORM = MAX( BNORM, ABS( D( J ) ) );
               }
               RESID = MAX( RESID, SASUM( N, WORK, 1 ) );
            } // 20
         } else {

            // B is lower bidiagonal.

            for (J = 1; J <= N; J++) { // 40
               for (I = 1; I <= N; I++) { // 30
                  WORK( N+I ) = S( I )*VT( I, J );
               } // 30
               sgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
               WORK( J ) = WORK( J ) + D( J );
               if ( J < N ) {
                  WORK( J+1 ) = WORK( J+1 ) + E( J );
                  BNORM = MAX( BNORM, ABS( D( J ) )+ABS( E( J ) ) );
               } else {
                  BNORM = MAX( BNORM, ABS( D( J ) ) );
               }
               RESID = MAX( RESID, SASUM( N, WORK, 1 ) );
            } // 40
         }
      } else {

         // B is diagonal.

         for (J = 1; J <= N; J++) { // 60
            for (I = 1; I <= N; I++) { // 50
               WORK( N+I ) = S( I )*VT( I, J );
            } // 50
            sgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
            WORK( J ) = WORK( J ) + D( J );
            RESID = MAX( RESID, SASUM( N, WORK, 1 ) );
         } // 60
         J = ISAMAX( N, D, 1 );
         BNORM = ABS( D( J ) );
      }

      // Compute norm(B - U * S * V') / ( n * norm(B) * EPS )

      EPS = SLAMCH( 'Precision' );

      if ( BNORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( BNORM >= RESID ) {
            RESID = ( RESID / BNORM ) / ( REAL( N )*EPS );
         } else {
            if ( BNORM < ONE ) {
               RESID = ( MIN( RESID, REAL( N )*BNORM ) / BNORM ) / ( REAL( N )*EPS );
            } else {
               RESID = MIN( RESID / BNORM, REAL( N ) ) / ( REAL( N )*EPS );
            }
         }
      }

      return;

      // End of SBDT03

      }
