      SUBROUTINE SGET34( RMAX, LMAX, NINFO, KNT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX;
      REAL               RMAX
      // ..
      // .. Array Arguments ..
      int                NINFO( 2 );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0 ;
      REAL               TWO, THREE
      const              TWO = 2.0E0, THREE = 3.0E0 ;
      int                LWORK;
      const              LWORK = 32 ;
      // ..
      // .. Local Scalars ..
      int                I, IA, IA11, IA12, IA21, IA22, IAM, IB, IC, IC11, IC12, IC21, IC22, ICM, INFO, J;
      REAL               BIGNUM, EPS, RES, SMLNUM, TNRM
      // ..
      // .. Local Arrays ..
      REAL               Q( 4, 4 ), RESULT( 2 ), T( 4, 4 ), T1( 4, 4 ), VAL( 9 ), VM( 2 ), WORK( LWORK )
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLAEXC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Get machine parameters

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // Set up test case parameters

      VAL( 1 ) = ZERO
      VAL( 2 ) = SQRT( SMLNUM )
      VAL( 3 ) = ONE
      VAL( 4 ) = TWO
      VAL( 5 ) = SQRT( BIGNUM )
      VAL( 6 ) = -SQRT( SMLNUM )
      VAL( 7 ) = -ONE
      VAL( 8 ) = -TWO
      VAL( 9 ) = -SQRT( BIGNUM )
      VM( 1 ) = ONE
      VM( 2 ) = ONE + TWO*EPS
      scopy(16, VAL( 4 ), 0, T( 1, 1 ), 1 );

      NINFO( 1 ) = 0
      NINFO( 2 ) = 0
      KNT = 0
      LMAX = 0
      RMAX = ZERO

      // Begin test loop

      for (IA = 1; IA <= 9; IA++) { // 40
         for (IAM = 1; IAM <= 2; IAM++) { // 30
            for (IB = 1; IB <= 9; IB++) { // 20
               for (IC = 1; IC <= 9; IC++) { // 10
                  T( 1, 1 ) = VAL( IA )*VM( IAM )
                  T( 2, 2 ) = VAL( IC )
                  T( 1, 2 ) = VAL( IB )
                  T( 2, 1 ) = ZERO
                  TNRM = MAX( ABS( T( 1, 1 ) ), ABS( T( 2, 2 ) ), ABS( T( 1, 2 ) ) )
                  scopy(16, T, 1, T1, 1 );
                  scopy(16, VAL( 1 ), 0, Q, 1 );
                  scopy(4, VAL( 3 ), 0, Q, 5 );
                  slaexc(.TRUE., 2, T, 4, Q, 4, 1, 1, 1, WORK, INFO )                   IF( INFO.NE.0 ) NINFO( INFO ) = NINFO( INFO ) + 1                   CALL SHST01( 2, 1, 2, T1, 4, T, 4, Q, 4, WORK, LWORK, RESULT );
                  RES = RESULT( 1 ) + RESULT( 2 )
                  IF( INFO.NE.0 ) RES = RES + ONE / EPS                   IF( T( 1, 1 ).NE.T1( 2, 2 ) ) RES = RES + ONE / EPS                   IF( T( 2, 2 ).NE.T1( 1, 1 ) ) RES = RES + ONE / EPS                   IF( T( 2, 1 ).NE.ZERO ) RES = RES + ONE / EPS
                  KNT = KNT + 1
                  if ( RES.GT.RMAX ) {
                     LMAX = KNT
                     RMAX = RES
                  }
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE

      for (IA = 1; IA <= 5; IA++) { // 110
         for (IAM = 1; IAM <= 2; IAM++) { // 100
            for (IB = 1; IB <= 5; IB++) { // 90
               for (IC11 = 1; IC11 <= 5; IC11++) { // 80
                  for (IC12 = 2; IC12 <= 5; IC12++) { // 70
                     for (IC21 = 2; IC21 <= 4; IC21++) { // 60
                        DO 50 IC22 = -1, 1, 2
                           T( 1, 1 ) = VAL( IA )*VM( IAM )
                           T( 1, 2 ) = VAL( IB )
                           T( 1, 3 ) = -TWO*VAL( IB )
                           T( 2, 1 ) = ZERO
                           T( 2, 2 ) = VAL( IC11 )
                           T( 2, 3 ) = VAL( IC12 )
                           T( 3, 1 ) = ZERO
                           T( 3, 2 ) = -VAL( IC21 )
                           T( 3, 3 ) = VAL( IC11 )*REAL( IC22 )
                           TNRM = MAX( ABS( T( 1, 1 ) ), ABS( T( 1, 2 ) ), ABS( T( 1, 3 ) ), ABS( T( 2, 2 ) ), ABS( T( 2, 3 ) ), ABS( T( 3, 2 ) ), ABS( T( 3, 3 ) ) )
                           scopy(16, T, 1, T1, 1 );
                           scopy(16, VAL( 1 ), 0, Q, 1 );
                           scopy(4, VAL( 3 ), 0, Q, 5 );
                           slaexc(.TRUE., 3, T, 4, Q, 4, 1, 1, 2, WORK, INFO )                            IF( INFO.NE.0 ) NINFO( INFO ) = NINFO( INFO ) + 1                            CALL SHST01( 3, 1, 3, T1, 4, T, 4, Q, 4, WORK, LWORK, RESULT );
                           RES = RESULT( 1 ) + RESULT( 2 )
                           if ( INFO.EQ.0 ) {
                              IF( T1( 1, 1 ).NE.T( 3, 3 ) ) RES = RES + ONE / EPS                               IF( T( 3, 1 ).NE.ZERO ) RES = RES + ONE / EPS                               IF( T( 3, 2 ).NE.ZERO ) RES = RES + ONE / EPS                               IF( T( 2, 1 ).NE.0 .AND. ( T( 1, 1 ).NE.T( 2, 2 ) .OR. SIGN( ONE, T( 1, 2 ) ).EQ.SIGN( ONE, T( 2, 1 ) ) ) ) RES = RES + ONE / EPS
                           }
                           KNT = KNT + 1
                           if ( RES.GT.RMAX ) {
                              LMAX = KNT
                              RMAX = RES
                           }
   50                   CONTINUE
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE

      for (IA11 = 1; IA11 <= 5; IA11++) { // 180
         for (IA12 = 2; IA12 <= 5; IA12++) { // 170
            for (IA21 = 2; IA21 <= 4; IA21++) { // 160
               DO 150 IA22 = -1, 1, 2
                  for (ICM = 1; ICM <= 2; ICM++) { // 140
                     for (IB = 1; IB <= 5; IB++) { // 130
                        for (IC = 1; IC <= 5; IC++) { // 120
                           T( 1, 1 ) = VAL( IA11 )
                           T( 1, 2 ) = VAL( IA12 )
                           T( 1, 3 ) = -TWO*VAL( IB )
                           T( 2, 1 ) = -VAL( IA21 )
                           T( 2, 2 ) = VAL( IA11 )*REAL( IA22 )
                           T( 2, 3 ) = VAL( IB )
                           T( 3, 1 ) = ZERO
                           T( 3, 2 ) = ZERO
                           T( 3, 3 ) = VAL( IC )*VM( ICM )
                           TNRM = MAX( ABS( T( 1, 1 ) ), ABS( T( 1, 2 ) ), ABS( T( 1, 3 ) ), ABS( T( 2, 2 ) ), ABS( T( 2, 3 ) ), ABS( T( 3, 2 ) ), ABS( T( 3, 3 ) ) )
                           scopy(16, T, 1, T1, 1 );
                           scopy(16, VAL( 1 ), 0, Q, 1 );
                           scopy(4, VAL( 3 ), 0, Q, 5 );
                           slaexc(.TRUE., 3, T, 4, Q, 4, 1, 2, 1, WORK, INFO )                            IF( INFO.NE.0 ) NINFO( INFO ) = NINFO( INFO ) + 1                            CALL SHST01( 3, 1, 3, T1, 4, T, 4, Q, 4, WORK, LWORK, RESULT );
                           RES = RESULT( 1 ) + RESULT( 2 )
                           if ( INFO.EQ.0 ) {
                              IF( T1( 3, 3 ).NE.T( 1, 1 ) ) RES = RES + ONE / EPS                               IF( T( 2, 1 ).NE.ZERO ) RES = RES + ONE / EPS                               IF( T( 3, 1 ).NE.ZERO ) RES = RES + ONE / EPS                               IF( T( 3, 2 ).NE.0 .AND. ( T( 2, 2 ).NE.T( 3, 3 ) .OR. SIGN( ONE, T( 2, 3 ) ).EQ.SIGN( ONE, T( 3, 2 ) ) ) ) RES = RES + ONE / EPS
                           }
                           KNT = KNT + 1
                           if ( RES.GT.RMAX ) {
                              LMAX = KNT
                              RMAX = RES
                           }
  120                   CONTINUE
  130                CONTINUE
  140             CONTINUE
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE

      for (IA11 = 1; IA11 <= 5; IA11++) { // 300
         for (IA12 = 2; IA12 <= 5; IA12++) { // 290
            for (IA21 = 2; IA21 <= 4; IA21++) { // 280
               DO 270 IA22 = -1, 1, 2
                  for (IB = 1; IB <= 5; IB++) { // 260
                     for (IC11 = 3; IC11 <= 4; IC11++) { // 250
                        for (IC12 = 3; IC12 <= 4; IC12++) { // 240
                           for (IC21 = 3; IC21 <= 4; IC21++) { // 230
                              DO 220 IC22 = -1, 1, 2
                                 for (ICM = 5; ICM <= 7; ICM++) { // 210
                                    IAM = 1
                                    T( 1, 1 ) = VAL( IA11 )*VM( IAM )
                                    T( 1, 2 ) = VAL( IA12 )*VM( IAM )
                                    T( 1, 3 ) = -TWO*VAL( IB )
                                    T( 1, 4 ) = HALF*VAL( IB )
                                    T( 2, 1 ) = -T( 1, 2 )*VAL( IA21 )
                                    T( 2, 2 ) = VAL( IA11 )* REAL( IA22 )*VM( IAM )
                                    T( 2, 3 ) = VAL( IB )
                                    T( 2, 4 ) = THREE*VAL( IB )
                                    T( 3, 1 ) = ZERO
                                    T( 3, 2 ) = ZERO
                                    T( 3, 3 ) = VAL( IC11 )* ABS( VAL( ICM ) )                                     T( 3, 4 ) = VAL( IC12 )* ABS( VAL( ICM ) )
                                    T( 4, 1 ) = ZERO
                                    T( 4, 2 ) = ZERO
                                    T( 4, 3 ) = -T( 3, 4 )*VAL( IC21 )* ABS( VAL( ICM ) )                                     T( 4, 4 ) = VAL( IC11 )* REAL( IC22 )* ABS( VAL( ICM ) )
                                    TNRM = ZERO
                                    for (I = 1; I <= 4; I++) { // 200
                                       for (J = 1; J <= 4; J++) { // 190
                                          TNRM = MAX( TNRM, ABS( T( I, J ) ) )
  190                                  CONTINUE
  200                               CONTINUE
                                    scopy(16, T, 1, T1, 1 );
                                    scopy(16, VAL( 1 ), 0, Q, 1 );
                                    scopy(4, VAL( 3 ), 0, Q, 5 );
                                    slaexc(.TRUE., 4, T, 4, Q, 4, 1, 2, 2, WORK, INFO )                                     IF( INFO.NE.0 ) NINFO( INFO ) = NINFO( INFO ) + 1                                     CALL SHST01( 4, 1, 4, T1, 4, T, 4, Q, 4, WORK, LWORK, RESULT );
                                    RES = RESULT( 1 ) + RESULT( 2 )
                                    if ( INFO.EQ.0 ) {
                                       IF( T( 3, 1 ).NE.ZERO ) RES = RES + ONE / EPS                                        IF( T( 4, 1 ).NE.ZERO ) RES = RES + ONE / EPS                                        IF( T( 3, 2 ).NE.ZERO ) RES = RES + ONE / EPS                                        IF( T( 4, 2 ).NE.ZERO ) RES = RES + ONE / EPS                                        IF( T( 2, 1 ).NE.0 .AND. ( T( 1, 1 ).NE.T( 2, 2 ) .OR. SIGN( ONE, T( 1, 2 ) ).EQ.SIGN( ONE, T( 2, 1 ) ) ) )RES = RES + ONE / EPS                                        IF( T( 4, 3 ).NE.0 .AND. ( T( 3, 3 ).NE.T( 4, 4 ) .OR. SIGN( ONE, T( 3, 4 ) ).EQ.SIGN( ONE, T( 4, 3 ) ) ) )RES = RES + ONE / EPS
                                    }
                                    KNT = KNT + 1
                                    if ( RES.GT.RMAX ) {
                                       LMAX = KNT
                                       RMAX = RES
                                    }
  210                            CONTINUE
  220                         CONTINUE
  230                      CONTINUE
  240                   CONTINUE
  250                CONTINUE
  260             CONTINUE
  270          CONTINUE
  280       CONTINUE
  290    CONTINUE
  300 CONTINUE

      RETURN

      // End of SGET34

      }
