      SUBROUTINE ZGTT05( TRANS, N, NRHS, DL, D, DU, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), FERR( * ), RESLTS( * );
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * ), XACT( LDXACT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IMAX, J, K, NZ;
      double             AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM;
      COMPLEX*16         ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IZAMAX, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N.LE.0 || NRHS.LE.0 ) {
         RESLTS( 1 ) = ZERO
         RESLTS( 2 ) = ZERO
         RETURN
      }

      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      NOTRAN = LSAME( TRANS, 'N' )
      NZ = 4

      // Test 1:  Compute the maximum of
         // norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO
      for (J = 1; J <= NRHS; J++) { // 30
         IMAX = IZAMAX( N, X( 1, J ), 1 )
         XNORM = MAX( CABS1( X( IMAX, J ) ), UNFL )
         DIFF = ZERO
         for (I = 1; I <= N; I++) { // 10
            DIFF = MAX( DIFF, CABS1( X( I, J )-XACT( I, J ) ) )
         } // 10

         if ( XNORM > ONE ) {
            GO TO 20
         } else if ( DIFF.LE.OVFL*XNORM ) {
            GO TO 20
         } else {
            ERRBND = ONE / EPS
            GO TO 30
         }

         } // 20
         if ( DIFF / XNORM.LE.FERR( J ) ) {
            ERRBND = MAX( ERRBND, ( DIFF / XNORM ) / FERR( J ) )
         } else {
            ERRBND = ONE / EPS
         }
      } // 30
      RESLTS( 1 ) = ERRBND

      // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
      // (*) = NZ*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )

      for (K = 1; K <= NRHS; K++) { // 60
         if ( NOTRAN ) {
            if ( N == 1 ) {
               AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 ) )*CABS1( X( 1, K ) )
            } else {
               AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 ) )*CABS1( X( 1, K ) ) + CABS1( DU( 1 ) )*CABS1( X( 2, K ) )
               for (I = 2; I <= N - 1; I++) { // 40
                  TMP = CABS1( B( I, K ) ) + CABS1( DL( I-1 ) )*CABS1( X( I-1, K ) ) + CABS1( D( I ) )*CABS1( X( I, K ) ) + CABS1( DU( I ) )*CABS1( X( I+1, K ) )
                  AXBI = MIN( AXBI, TMP )
               } // 40
               TMP = CABS1( B( N, K ) ) + CABS1( DL( N-1 ) )* CABS1( X( N-1, K ) ) + CABS1( D( N ) )* CABS1( X( N, K ) )
               AXBI = MIN( AXBI, TMP )
            }
         } else {
            if ( N == 1 ) {
               AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 ) )*CABS1( X( 1, K ) )
            } else {
               AXBI = CABS1( B( 1, K ) ) + CABS1( D( 1 ) )*CABS1( X( 1, K ) ) + CABS1( DL( 1 ) )*CABS1( X( 2, K ) )
               for (I = 2; I <= N - 1; I++) { // 50
                  TMP = CABS1( B( I, K ) ) + CABS1( DU( I-1 ) )*CABS1( X( I-1, K ) ) + CABS1( D( I ) )*CABS1( X( I, K ) ) + CABS1( DL( I ) )*CABS1( X( I+1, K ) )
                  AXBI = MIN( AXBI, TMP )
               } // 50
               TMP = CABS1( B( N, K ) ) + CABS1( DU( N-1 ) )* CABS1( X( N-1, K ) ) + CABS1( D( N ) )* CABS1( X( N, K ) )
               AXBI = MIN( AXBI, TMP )
            }
         }
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / MAX( AXBI, NZ*UNFL ) )
         if ( K == 1 ) {
            RESLTS( 2 ) = TMP
         } else {
            RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
         }
      } // 60

      RETURN

      // End of ZGTT05

      }
