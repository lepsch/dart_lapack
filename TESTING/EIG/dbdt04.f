      void dbdt04(UPLO, N, D, E, S, NS, U, LDU, VT, LDVT, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDU, LDVT, N, NS;
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
      int                I, J, K;
      double             BNORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                IDAMAX;
      //- double             DASUM, DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DASUM, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      RESID = ZERO;
      if (N <= 0 || NS <= 0) return;

      EPS = DLAMCH( 'Precision' );

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
         BNORM = ( D( 1 ) ).abs();
         for (I = 2; I <= N; I++) { // 30
            BNORM = max( BNORM, ( D( I ) ).abs()+( E( I-1 ) ) ).abs();
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
         BNORM = ( D( N ) ).abs();
         for (I = 1; I <= N-1; I++) { // 60
            BNORM = max( BNORM, ( D( I ) ).abs()+( E( I ) ) ).abs();
         } // 60
      }

      dgemm('T', 'N', NS, NS, N, -ONE, U, LDU, WORK( 1 ), N, ZERO, WORK( 1+N*NS ), NS );

      // norm(S - U' * B * V)

      K = N*NS;
      for (I = 1; I <= NS; I++) { // 70
         WORK( K+I ) =  WORK( K+I ) + S( I );
         RESID = max( RESID, DASUM( NS, WORK( K+1 ), 1 ) );
         K = K + NS;
      } // 70

      if ( BNORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( BNORM >= RESID ) {
            RESID = ( RESID / BNORM ) / ( DBLE( N )*EPS );
         } else {
            if ( BNORM < ONE ) {
               RESID = ( min( RESID, DBLE( N )*BNORM ) / BNORM ) / ( DBLE( N )*EPS );
            } else {
               RESID = min( RESID / BNORM, DBLE( N ) ) / ( DBLE( N )*EPS );
            }
         }
      }

      return;
      }
