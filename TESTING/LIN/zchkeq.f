      SUBROUTINE ZCHKEQ( THRESH, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NOUT;
      double             THRESH;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0D0, ONE = 1.0D+0, TEN = 1.0D1 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D0, 0.0D0 ) ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D0, 0.0D0 ) ;
      int                NSZ, NSZB;
      const              NSZ = 5, NSZB = 3*NSZ-2 ;
      int                NSZP, NPOW;
      const              NSZP = ( NSZ*( NSZ+1 ) ) / 2, NPOW = 2*NSZ+1 ;
      // ..
      // .. Local Scalars ..
      bool               OK;
      String             PATH;
      int                I, INFO, J, KL, KU, M, N;
      double             CCOND, EPS, NORM, RATIO, RCMAX, RCMIN, RCOND;
      // ..
      // .. Local Arrays ..
      double             C( NSZ ), POW( NPOW ), R( NSZ ), RESLTS( 5 ), RPOW( NPOW );
      COMPLEX*16         A( NSZ, NSZ ), AB( NSZB, NSZ ), AP( NSZP )
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGBEQU, ZGEEQU, ZPBEQU, ZPOEQU, ZPPEQU
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'EQ'

      EPS = DLAMCH( 'P' )
      for (I = 1; I <= 5; I++) { // 10
         RESLTS( I ) = ZERO
   10 CONTINUE
      for (I = 1; I <= NPOW; I++) { // 20
         POW( I ) = TEN**( I-1 )
         RPOW( I ) = ONE / POW( I )
   20 CONTINUE

      // Test ZGEEQU

      for (N = 0; N <= NSZ; N++) { // 80
         for (M = 0; M <= NSZ; M++) { // 70

            for (J = 1; J <= NSZ; J++) { // 40
               for (I = 1; I <= NSZ; I++) { // 30
                  if ( I.LE.M .AND. J.LE.N ) {
                     A( I, J ) = POW( I+J+1 )*( -1 )**( I+J )
                  } else {
                     A( I, J ) = CZERO
                  }
   30          CONTINUE
   40       CONTINUE

            zgeequ(M, N, A, NSZ, R, C, RCOND, CCOND, NORM, INFO );

            if ( INFO.NE.0 ) {
               RESLTS( 1 ) = ONE
            } else {
               if ( N.NE.0 .AND. M.NE.0 ) {
                  RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( RCOND-RPOW( M ) ) / RPOW( M ) ) )                   RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( CCOND-RPOW( N ) ) / RPOW( N ) ) )                   RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( NORM-POW( N+M+1 ) ) / POW( N+M+ 1 ) ) )
                  for (I = 1; I <= M; I++) { // 50
                     RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( R( I )-RPOW( I+N+1 ) ) / RPOW( I+N+1 ) ) )
   50             CONTINUE
                  for (J = 1; J <= N; J++) { // 60
                     RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( C( J )-POW( N-J+1 ) ) / POW( N-J+1 ) ) )
   60             CONTINUE
               }
            }

   70    CONTINUE
   80 CONTINUE

      // Test with zero rows and columns

      for (J = 1; J <= NSZ; J++) { // 90
         A( MAX( NSZ-1, 1 ), J ) = CZERO
   90 CONTINUE
      zgeequ(NSZ, NSZ, A, NSZ, R, C, RCOND, CCOND, NORM, INFO );
      IF( INFO.NE.MAX( NSZ-1, 1 ) ) RESLTS( 1 ) = ONE

      for (J = 1; J <= NSZ; J++) { // 100
         A( MAX( NSZ-1, 1 ), J ) = CONE
  100 CONTINUE
      for (I = 1; I <= NSZ; I++) { // 110
         A( I, MAX( NSZ-1, 1 ) ) = CZERO
  110 CONTINUE
      zgeequ(NSZ, NSZ, A, NSZ, R, C, RCOND, CCOND, NORM, INFO );
      IF( INFO.NE.NSZ+MAX( NSZ-1, 1 ) ) RESLTS( 1 ) = ONE
      RESLTS( 1 ) = RESLTS( 1 ) / EPS

      // Test ZGBEQU

      for (N = 0; N <= NSZ; N++) { // 250
         for (M = 0; M <= NSZ; M++) { // 240
            DO 230 KL = 0, MAX( M-1, 0 )
               DO 220 KU = 0, MAX( N-1, 0 )

                  for (J = 1; J <= NSZ; J++) { // 130
                     for (I = 1; I <= NSZB; I++) { // 120
                        AB( I, J ) = CZERO
  120                CONTINUE
  130             CONTINUE
                  for (J = 1; J <= N; J++) { // 150
                     for (I = 1; I <= M; I++) { // 140
                        IF( I.LE.MIN( M, J+KL ) .AND. I.GE. MAX( 1, J-KU ) .AND. J.LE.N ) THEN                            AB( KU+1+I-J, J ) = POW( I+J+1 )* ( -1 )**( I+J )
                        }
  140                CONTINUE
  150             CONTINUE

                  zgbequ(M, N, KL, KU, AB, NSZB, R, C, RCOND, CCOND, NORM, INFO );

                  if ( INFO.NE.0 ) {
                     if ( .NOT.( ( N+KL.LT.M .AND. INFO.EQ.N+KL+1 ) .OR. ( M+KU.LT.N .AND. INFO.EQ.2*M+KU+1 ) ) ) {
                        RESLTS( 2 ) = ONE
                     }
                  } else {
                     if ( N.NE.0 .AND. M.NE.0 ) {

                        RCMIN = R( 1 )
                        RCMAX = R( 1 )
                        for (I = 1; I <= M; I++) { // 160
                           RCMIN = MIN( RCMIN, R( I ) )
                           RCMAX = MAX( RCMAX, R( I ) )
  160                   CONTINUE
                        RATIO = RCMIN / RCMAX
                        RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ( RCOND-RATIO ) / RATIO ) )

                        RCMIN = C( 1 )
                        RCMAX = C( 1 )
                        for (J = 1; J <= N; J++) { // 170
                           RCMIN = MIN( RCMIN, C( J ) )
                           RCMAX = MAX( RCMAX, C( J ) )
  170                   CONTINUE
                        RATIO = RCMIN / RCMAX
                        RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ( CCOND-RATIO ) / RATIO ) )

                        RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ( NORM-POW( N+M+1 ) ) / POW( N+M+1 ) ) )
                        for (I = 1; I <= M; I++) { // 190
                           RCMAX = ZERO
                           for (J = 1; J <= N; J++) { // 180
                              if ( I.LE.J+KL .AND. I.GE.J-KU ) {
                                 RATIO = ABS( R( I )*POW( I+J+1 )* C( J ) )
                                 RCMAX = MAX( RCMAX, RATIO )
                              }
  180                      CONTINUE
                           RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ONE-RCMAX ) )
  190                   CONTINUE

                        for (J = 1; J <= N; J++) { // 210
                           RCMAX = ZERO
                           for (I = 1; I <= M; I++) { // 200
                              if ( I.LE.J+KL .AND. I.GE.J-KU ) {
                                 RATIO = ABS( R( I )*POW( I+J+1 )* C( J ) )
                                 RCMAX = MAX( RCMAX, RATIO )
                              }
  200                      CONTINUE
                           RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ONE-RCMAX ) )
  210                   CONTINUE
                     }
                  }

  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RESLTS( 2 ) = RESLTS( 2 ) / EPS

      // Test ZPOEQU

      for (N = 0; N <= NSZ; N++) { // 290

         for (I = 1; I <= NSZ; I++) { // 270
            for (J = 1; J <= NSZ; J++) { // 260
               if ( I.LE.N .AND. J.EQ.I ) {
                  A( I, J ) = POW( I+J+1 )*( -1 )**( I+J )
               } else {
                  A( I, J ) = CZERO
               }
  260       CONTINUE
  270    CONTINUE

         zpoequ(N, A, NSZ, R, RCOND, NORM, INFO );

         if ( INFO.NE.0 ) {
            RESLTS( 3 ) = ONE
         } else {
            if ( N.NE.0 ) {
               RESLTS( 3 ) = MAX( RESLTS( 3 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                RESLTS( 3 ) = MAX( RESLTS( 3 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) )
               for (I = 1; I <= N; I++) { // 280
                  RESLTS( 3 ) = MAX( RESLTS( 3 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ 1 ) ) )
  280          CONTINUE
            }
         }
  290 CONTINUE
      A( MAX( NSZ-1, 1 ), MAX( NSZ-1, 1 ) ) = -CONE
      zpoequ(NSZ, A, NSZ, R, RCOND, NORM, INFO );
      IF( INFO.NE.MAX( NSZ-1, 1 ) ) RESLTS( 3 ) = ONE
      RESLTS( 3 ) = RESLTS( 3 ) / EPS

      // Test ZPPEQU

      for (N = 0; N <= NSZ; N++) { // 360

         // Upper triangular packed storage

         DO 300 I = 1, ( N*( N+1 ) ) / 2
            AP( I ) = CZERO
  300    CONTINUE
         for (I = 1; I <= N; I++) { // 310
            AP( ( I*( I+1 ) ) / 2 ) = POW( 2*I+1 )
  310    CONTINUE

         zppequ('U', N, AP, R, RCOND, NORM, INFO );

         if ( INFO.NE.0 ) {
            RESLTS( 4 ) = ONE
         } else {
            if ( N.NE.0 ) {
               RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) )
               for (I = 1; I <= N; I++) { // 320
                  RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ 1 ) ) )
  320          CONTINUE
            }
         }

         // Lower triangular packed storage

         DO 330 I = 1, ( N*( N+1 ) ) / 2
            AP( I ) = CZERO
  330    CONTINUE
         J = 1
         for (I = 1; I <= N; I++) { // 340
            AP( J ) = POW( 2*I+1 )
            J = J + ( N-I+1 )
  340    CONTINUE

         zppequ('L', N, AP, R, RCOND, NORM, INFO );

         if ( INFO.NE.0 ) {
            RESLTS( 4 ) = ONE
         } else {
            if ( N.NE.0 ) {
               RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) )
               for (I = 1; I <= N; I++) { // 350
                  RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ 1 ) ) )
  350          CONTINUE
            }
         }

  360 CONTINUE
      I = ( NSZ*( NSZ+1 ) ) / 2 - 2
      AP( I ) = -CONE
      zppequ('L', NSZ, AP, R, RCOND, NORM, INFO );
      IF( INFO.NE.MAX( NSZ-1, 1 ) ) RESLTS( 4 ) = ONE
      RESLTS( 4 ) = RESLTS( 4 ) / EPS

      // Test ZPBEQU

      for (N = 0; N <= NSZ; N++) { // 460
         DO 450 KL = 0, MAX( N-1, 0 )

            // Test upper triangular storage

            for (J = 1; J <= NSZ; J++) { // 380
               for (I = 1; I <= NSZB; I++) { // 370
                  AB( I, J ) = CZERO
  370          CONTINUE
  380       CONTINUE
            for (J = 1; J <= N; J++) { // 390
               AB( KL+1, J ) = POW( 2*J+1 )
  390       CONTINUE

            zpbequ('U', N, KL, AB, NSZB, R, RCOND, NORM, INFO );

            if ( INFO.NE.0 ) {
               RESLTS( 5 ) = ONE
            } else {
               if ( N.NE.0 ) {
                  RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                   RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) )
                  for (I = 1; I <= N; I++) { // 400
                     RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+1 ) ) )
  400             CONTINUE
               }
            }
            if ( N.NE.0 ) {
               AB( KL+1, MAX( N-1, 1 ) ) = -CONE
               zpbequ('U', N, KL, AB, NSZB, R, RCOND, NORM, INFO );
               IF( INFO.NE.MAX( N-1, 1 ) ) RESLTS( 5 ) = ONE
            }

            // Test lower triangular storage

            for (J = 1; J <= NSZ; J++) { // 420
               for (I = 1; I <= NSZB; I++) { // 410
                  AB( I, J ) = CZERO
  410          CONTINUE
  420       CONTINUE
            for (J = 1; J <= N; J++) { // 430
               AB( 1, J ) = POW( 2*J+1 )
  430       CONTINUE

            zpbequ('L', N, KL, AB, NSZB, R, RCOND, NORM, INFO );

            if ( INFO.NE.0 ) {
               RESLTS( 5 ) = ONE
            } else {
               if ( N.NE.0 ) {
                  RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                   RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) )
                  for (I = 1; I <= N; I++) { // 440
                     RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+1 ) ) )
  440             CONTINUE
               }
            }
            if ( N.NE.0 ) {
               AB( 1, MAX( N-1, 1 ) ) = -CONE
               zpbequ('L', N, KL, AB, NSZB, R, RCOND, NORM, INFO );
               IF( INFO.NE.MAX( N-1, 1 ) ) RESLTS( 5 ) = ONE
            }
  450    CONTINUE
  460 CONTINUE
      RESLTS( 5 ) = RESLTS( 5 ) / EPS
      OK = ( RESLTS( 1 ).LE.THRESH ) .AND. ( RESLTS( 2 ).LE.THRESH ) .AND. ( RESLTS( 3 ).LE.THRESH ) .AND. ( RESLTS( 4 ).LE.THRESH ) .AND. ( RESLTS( 5 ).LE.THRESH )
      WRITE( NOUT, FMT = * )
      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH
      } else {
         IF( RESLTS( 1 ).GT.THRESH ) WRITE( NOUT, FMT = 9998 )RESLTS( 1 ), THRESH          IF( RESLTS( 2 ).GT.THRESH ) WRITE( NOUT, FMT = 9997 )RESLTS( 2 ), THRESH          IF( RESLTS( 3 ).GT.THRESH ) WRITE( NOUT, FMT = 9996 )RESLTS( 3 ), THRESH          IF( RESLTS( 4 ).GT.THRESH ) WRITE( NOUT, FMT = 9995 )RESLTS( 4 ), THRESH          IF( RESLTS( 5 ).GT.THRESH ) WRITE( NOUT, FMT = 9994 )RESLTS( 5 ), THRESH
      }
 9999 FORMAT( 1X, 'All tests for ', A3, ' routines passed the threshold' )
 9998 FORMAT( ' ZGEEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 )
 9997 FORMAT( ' ZGBEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 )
 9996 FORMAT( ' ZPOEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 )
 9995 FORMAT( ' ZPPEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 )
 9994 FORMAT( ' ZPBEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 )
      RETURN

      // End of ZCHKEQ

      }
