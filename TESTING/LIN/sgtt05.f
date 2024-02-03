      SUBROUTINE SGTT05( TRANS, N, NRHS, DL, D, DU, B, LDB, X, LDX, XACT, LDXACT, FERR, BERR, RESLTS )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, LDXACT, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), BERR( * ), D( * ), DL( * ), DU( * ), FERR( * ), RESLTS( * ), X( LDX, * ), XACT( LDXACT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IMAX, J, K, NZ;
      REAL               AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ISAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0 or NRHS = 0.

      if ( N.LE.0 .OR. NRHS.LE.0 ) {
         RESLTS( 1 ) = ZERO
         RESLTS( 2 ) = ZERO
         RETURN
      }

      EPS = SLAMCH( 'Epsilon' )
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      NOTRAN = LSAME( TRANS, 'N' )
      NZ = 4

      // Test 1:  Compute the maximum of
         // norm(X - XACT) / ( norm(X) * FERR )
      // over all the vectors X and XACT using the infinity-norm.

      ERRBND = ZERO
      for (J = 1; J <= NRHS; J++) { // 30
         IMAX = ISAMAX( N, X( 1, J ), 1 )
         XNORM = MAX( ABS( X( IMAX, J ) ), UNFL )
         DIFF = ZERO
         for (I = 1; I <= N; I++) { // 10
            DIFF = MAX( DIFF, ABS( X( I, J )-XACT( I, J ) ) )
   10    CONTINUE

         if ( XNORM.GT.ONE ) {
            GO TO 20
         } else if ( DIFF.LE.OVFL*XNORM ) {
            GO TO 20
         } else {
            ERRBND = ONE / EPS
            GO TO 30
         }

   20    CONTINUE
         if ( DIFF / XNORM.LE.FERR( J ) ) {
            ERRBND = MAX( ERRBND, ( DIFF / XNORM ) / FERR( J ) )
         } else {
            ERRBND = ONE / EPS
         }
   30 CONTINUE
      RESLTS( 1 ) = ERRBND

      // Test 2:  Compute the maximum of BERR / ( NZ*EPS + (*) ), where
      // (*) = NZ*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )

      for (K = 1; K <= NRHS; K++) { // 60
         if ( NOTRAN ) {
            if ( N.EQ.1 ) {
               AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) )
            } else {
               AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) ) + ABS( DU( 1 )*X( 2, K ) )
               DO 40 I = 2, N - 1
                  TMP = ABS( B( I, K ) ) + ABS( DL( I-1 )*X( I-1, K ) ) + ABS( D( I )*X( I, K ) ) + ABS( DU( I )*X( I+1, K ) )
                  AXBI = MIN( AXBI, TMP )
   40          CONTINUE
               TMP = ABS( B( N, K ) ) + ABS( DL( N-1 )*X( N-1, K ) ) + ABS( D( N )*X( N, K ) )
               AXBI = MIN( AXBI, TMP )
            }
         } else {
            if ( N.EQ.1 ) {
               AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) )
            } else {
               AXBI = ABS( B( 1, K ) ) + ABS( D( 1 )*X( 1, K ) ) + ABS( DL( 1 )*X( 2, K ) )
               DO 50 I = 2, N - 1
                  TMP = ABS( B( I, K ) ) + ABS( DU( I-1 )*X( I-1, K ) ) + ABS( D( I )*X( I, K ) ) + ABS( DL( I )*X( I+1, K ) )
                  AXBI = MIN( AXBI, TMP )
   50          CONTINUE
               TMP = ABS( B( N, K ) ) + ABS( DU( N-1 )*X( N-1, K ) ) + ABS( D( N )*X( N, K ) )
               AXBI = MIN( AXBI, TMP )
            }
         }
         TMP = BERR( K ) / ( NZ*EPS+NZ*UNFL / MAX( AXBI, NZ*UNFL ) )
         if ( K.EQ.1 ) {
            RESLTS( 2 ) = TMP
         } else {
            RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
         }
   60 CONTINUE

      RETURN

      // End of SGTT05

      }
