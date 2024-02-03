      void sbdt04(UPLO, N, D, E, S, NS, U, LDU, VT, LDVT, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDU, LDVT, N, NS;
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
      int                I, J, K;
      REAL               BNORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SASUM, SLAMCH;
      // EXTERNAL LSAME, ISAMAX, SASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      RESID = ZERO;
      if (N <= 0 || NS <= 0) return;

      EPS = SLAMCH( 'Precision' );

      // Compute S - U' * B * V.

      BNORM = ZERO;

      if ( LSAME( UPLO, 'U' ) ) {

         // B is upper bidiagonal.

         K = 0;
         for (I = 1; I <= NS; I++) { // 20
            for (J = 1; J <= N-1; J++) { // 10
               K = K + 1;
               WORK( K ) = D( J )*VT( I, J ) + E( J )*VT( I, J+1 );
            } // 10
            K = K + 1;
            WORK( K ) = D( N )*VT( I, N );
         } // 20
         BNORM = ABS( D( 1 ) );
         for (I = 2; I <= N; I++) { // 30
            BNORM = max( BNORM, ABS( D( I ) )+ABS( E( I-1 ) ) );
         } // 30
      } else {

         // B is lower bidiagonal.

         K = 0;
         for (I = 1; I <= NS; I++) { // 50
            K = K + 1;
            WORK( K ) = D( 1 )*VT( I, 1 );
            for (J = 1; J <= N-1; J++) { // 40
               K = K + 1;
               WORK( K ) = E( J )*VT( I, J ) + D( J+1 )*VT( I, J+1 );
            } // 40
         } // 50
         BNORM = ABS( D( N ) );
         for (I = 1; I <= N-1; I++) { // 60
            BNORM = max( BNORM, ABS( D( I ) )+ABS( E( I ) ) );
         } // 60
      }

      sgemm('T', 'N', NS, NS, N, -ONE, U, LDU, WORK( 1 ), N, ZERO, WORK( 1+N*NS ), NS );

      // norm(S - U' * B * V)

      K = N*NS;
      for (I = 1; I <= NS; I++) { // 70
         WORK( K+I ) =  WORK( K+I ) + S( I );
         RESID = max( RESID, SASUM( NS, WORK( K+1 ), 1 ) );
         K = K + NS;
      } // 70

      if ( BNORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( BNORM >= RESID ) {
            RESID = ( RESID / BNORM ) / ( REAL( N )*EPS );
         } else {
            if ( BNORM < ONE ) {
               RESID = ( min( RESID, REAL( N )*BNORM ) / BNORM ) / ( REAL( N )*EPS );
            } else {
               RESID = min( RESID / BNORM, REAL( N ) ) / ( REAL( N )*EPS );
            }
         }
      }

      return;
      }
