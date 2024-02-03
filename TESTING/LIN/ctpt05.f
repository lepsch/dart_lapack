      SUBROUTINE CTPT05( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               BERR( * ), FERR( * ), RESLTS( * )
      COMPLEX            AP( * ), B( LDB, * ), X( LDX, * ), XACT( LDXACT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, UNIT, UPPER;
      int                I, IFU, IMAX, J, JC, K;
      REAL               AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM
      COMPLEX            ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ICAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESLTS( 1 ) = ZERO
         RESLTS( 2 ) = ZERO
         RETURN
      }

      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      UNIT = LSAME( DIAG, 'U' )

      // Test 1:  Compute the maximum of
         // norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO
      for (J = 1; J <= NRHS; J++) { // 30
         IMAX = ICAMAX( N, X( 1, J ), 1 )
         XNORM = MAX( CABS1( X( IMAX, J ) ), UNFL )
         DIFF = ZERO
         for (I = 1; I <= N; I++) { // 10
            DIFF = MAX( DIFF, CABS1( X( I, J )-XACT( I, J ) ) )
         } // 10

         if ( XNORM > ONE ) {
            GO TO 20
         } else if ( DIFF <= OVFL*XNORM ) {
            GO TO 20
         } else {
            ERRBND = ONE / EPS
            GO TO 30
         }

         } // 20
         if ( DIFF / XNORM <= FERR( J ) ) {
            ERRBND = MAX( ERRBND, ( DIFF / XNORM ) / FERR( J ) )
         } else {
            ERRBND = ONE / EPS
         }
      } // 30
      RESLTS( 1 ) = ERRBND

      // Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
      // (*) = (n+1)*UNFL / (min_i (abs(A)*abs(X) +abs(b))_i )

      IFU = 0
      if (UNIT) IFU = 1;
      for (K = 1; K <= NRHS; K++) { // 90
         for (I = 1; I <= N; I++) { // 80
            TMP = CABS1( B( I, K ) )
            if ( UPPER ) {
               JC = ( ( I-1 )*I ) / 2
               if ( !NOTRAN ) {
                  for (J = 1; J <= I - IFU; J++) { // 40
                     TMP = TMP + CABS1( AP( JC+J ) )*CABS1( X( J, K ) )
                  } // 40
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
               } else {
                  JC = JC + I
                  if ( UNIT ) {
                     TMP = TMP + CABS1( X( I, K ) )
                     JC = JC + I
                  }
                  for (J = I + IFU; J <= N; J++) { // 50
                     TMP = TMP + CABS1( AP( JC ) )*CABS1( X( J, K ) )
                     JC = JC + J
                  } // 50
               }
            } else {
               if ( NOTRAN ) {
                  JC = I
                  for (J = 1; J <= I - IFU; J++) { // 60
                     TMP = TMP + CABS1( AP( JC ) )*CABS1( X( J, K ) )
                     JC = JC + N - J
                  } // 60
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
               } else {
                  JC = ( I-1 )*( N-I ) + ( I*( I+1 ) ) / 2
                  if (UNIT) TMP = TMP + CABS1( X( I, K ) );
                  for (J = I + IFU; J <= N; J++) { // 70
                     TMP = TMP + CABS1( AP( JC+J-I ) )* CABS1( X( J, K ) )
                  } // 70
               }
            }
            if ( I == 1 ) {
               AXBI = TMP
            } else {
               AXBI = MIN( AXBI, TMP )
            }
         } // 80
         TMP = BERR( K ) / ( ( N+1 )*EPS+( N+1 )*UNFL / MAX( AXBI, ( N+1 )*UNFL ) )
         if ( K == 1 ) {
            RESLTS( 2 ) = TMP
         } else {
            RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
         }
      } // 90

      RETURN

      // End of CTPT05

      }
