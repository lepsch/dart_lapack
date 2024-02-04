      void dbdt03(UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                KD, LDU, LDVT, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             BNORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                idamax;
      //- double             DASUM, DLAMCH;
      // EXTERNAL lsame, idamax, DASUM, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      RESID = ZERO;
      if (N <= 0) return;

      // Compute B - U * S * V' one column at a time.

      BNORM = ZERO;
      if ( KD >= 1 ) {

         // B is bidiagonal.

         if ( lsame( UPLO, 'U' ) ) {

            // B is upper bidiagonal.

            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= N; I++) { // 10
                  WORK[N+I] = S( I )*VT( I, J );
               } // 10
               dgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
               WORK[J] = WORK( J ) + D( J );
               if ( J > 1 ) {
                  WORK[J-1] = WORK( J-1 ) + E( J-1 );
                  BNORM = max( BNORM, ( D( J ) ).abs()+( E( J-1 ) ) ).abs();
               } else {
                  BNORM = max( BNORM, ( D( J ) ) ).abs();
               }
               RESID = max( RESID, DASUM( N, WORK, 1 ) );
            } // 20
         } else {

            // B is lower bidiagonal.

            for (J = 1; J <= N; J++) { // 40
               for (I = 1; I <= N; I++) { // 30
                  WORK[N+I] = S( I )*VT( I, J );
               } // 30
               dgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
               WORK[J] = WORK( J ) + D( J );
               if ( J < N ) {
                  WORK[J+1] = WORK( J+1 ) + E( J );
                  BNORM = max( BNORM, ( D( J ) ).abs()+( E( J ) ) ).abs();
               } else {
                  BNORM = max( BNORM, ( D( J ) ) ).abs();
               }
               RESID = max( RESID, DASUM( N, WORK, 1 ) );
            } // 40
         }
      } else {

         // B is diagonal.

         for (J = 1; J <= N; J++) { // 60
            for (I = 1; I <= N; I++) { // 50
               WORK[N+I] = S( I )*VT( I, J );
            } // 50
            dgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
            WORK[J] = WORK( J ) + D( J );
            RESID = max( RESID, DASUM( N, WORK, 1 ) );
         } // 60
         J = idamax( N, D, 1 );
         BNORM = ( D( J ) ).abs();
      }

      // Compute norm(B - U * S * V') / ( n * norm(B) * EPS )

      EPS = DLAMCH( 'Precision' );

      if ( BNORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( BNORM >= RESID ) {
            RESID = ( RESID / BNORM ) / ( N.toDouble()*EPS );
         } else {
            if ( BNORM < ONE ) {
               RESID = ( min( RESID, (N).toDouble()*BNORM ) / BNORM ) / ( N.toDouble()*EPS );
            } else {
               RESID = min( RESID / BNORM, (N).toDouble() ) / ( N.toDouble()*EPS );
            }
         }
      }

      return;
      }