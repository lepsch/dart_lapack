      SUBROUTINE DCHKEQ( THRESH, NOUT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NOUT;
      double             THRESH;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 1.0e1 ;
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
      double             A( NSZ, NSZ ), AB( NSZB, NSZ ), AP( NSZP ), C( NSZ ), POW( NPOW ), R( NSZ ), RESLTS( 5 ), RPOW( NPOW );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGBEQU, DGEEQU, DPBEQU, DPOEQU, DPPEQU
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'EQ';

      EPS = DLAMCH( 'P' );
      for (I = 1; I <= 5; I++) { // 10
         RESLTS( I ) = ZERO;
      } // 10
      for (I = 1; I <= NPOW; I++) { // 20
         POW( I ) = TEN**( I-1 );
         RPOW( I ) = ONE / POW( I );
      } // 20

      // Test DGEEQU

      for (N = 0; N <= NSZ; N++) { // 80
         for (M = 0; M <= NSZ; M++) { // 70

            for (J = 1; J <= NSZ; J++) { // 40
               for (I = 1; I <= NSZ; I++) { // 30
                  if ( I <= M && J <= N ) {
                     A( I, J ) = POW( I+J+1 )*( -1 )**( I+J );
                  } else {
                     A( I, J ) = ZERO;
                  }
               } // 30
            } // 40

            dgeequ(M, N, A, NSZ, R, C, RCOND, CCOND, NORM, INFO );

            if ( INFO != 0 ) {
               RESLTS( 1 ) = ONE;
            } else {
               if ( N != 0 && M != 0 ) {
                  RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( RCOND-RPOW( M ) ) / RPOW( M ) ) )                   RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( CCOND-RPOW( N ) ) / RPOW( N ) ) )                   RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( NORM-POW( N+M+1 ) ) / POW( N+M+ 1 ) ) );
                  for (I = 1; I <= M; I++) { // 50
                     RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( R( I )-RPOW( I+N+1 ) ) / RPOW( I+N+1 ) ) );
                  } // 50
                  for (J = 1; J <= N; J++) { // 60
                     RESLTS( 1 ) = MAX( RESLTS( 1 ), ABS( ( C( J )-POW( N-J+1 ) ) / POW( N-J+1 ) ) );
                  } // 60
               }
            }

         } // 70
      } // 80

      // Test with zero rows and columns

      for (J = 1; J <= NSZ; J++) { // 90
         A( MAX( NSZ-1, 1 ), J ) = ZERO;
      } // 90
      dgeequ(NSZ, NSZ, A, NSZ, R, C, RCOND, CCOND, NORM, INFO );
      if( INFO != MAX( NSZ-1, 1 ) ) RESLTS( 1 ) = ONE;

      for (J = 1; J <= NSZ; J++) { // 100
         A( MAX( NSZ-1, 1 ), J ) = ONE;
      } // 100
      for (I = 1; I <= NSZ; I++) { // 110
         A( I, MAX( NSZ-1, 1 ) ) = ZERO;
      } // 110
      dgeequ(NSZ, NSZ, A, NSZ, R, C, RCOND, CCOND, NORM, INFO );
      if( INFO != NSZ+MAX( NSZ-1, 1 ) ) RESLTS( 1 ) = ONE;
      RESLTS( 1 ) = RESLTS( 1 ) / EPS;

      // Test DGBEQU

      for (N = 0; N <= NSZ; N++) { // 250
         for (M = 0; M <= NSZ; M++) { // 240
            DO 230 KL = 0, MAX( M-1, 0 );
               DO 220 KU = 0, MAX( N-1, 0 );

                  for (J = 1; J <= NSZ; J++) { // 130
                     for (I = 1; I <= NSZB; I++) { // 120
                        AB( I, J ) = ZERO;
                     } // 120
                  } // 130
                  for (J = 1; J <= N; J++) { // 150
                     for (I = 1; I <= M; I++) { // 140
                        if( I <= MIN( M, J+KL ) && I >= MAX( 1, J-KU ) && J <= N ) {
                           AB( KU+1+I-J, J ) = POW( I+J+1 )* ( -1 )**( I+J );
                        }
                     } // 140
                  } // 150

                  dgbequ(M, N, KL, KU, AB, NSZB, R, C, RCOND, CCOND, NORM, INFO );

                  if ( INFO != 0 ) {
                     if ( !( ( N+KL < M && INFO == N+KL+1 ) || ( M+KU < N && INFO == 2*M+KU+1 ) ) ) {
                        RESLTS( 2 ) = ONE;
                     }
                  } else {
                     if ( N != 0 && M != 0 ) {

                        RCMIN = R( 1 );
                        RCMAX = R( 1 );
                        for (I = 1; I <= M; I++) { // 160
                           RCMIN = MIN( RCMIN, R( I ) );
                           RCMAX = MAX( RCMAX, R( I ) );
                        } // 160
                        RATIO = RCMIN / RCMAX;
                        RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ( RCOND-RATIO ) / RATIO ) );

                        RCMIN = C( 1 );
                        RCMAX = C( 1 );
                        for (J = 1; J <= N; J++) { // 170
                           RCMIN = MIN( RCMIN, C( J ) );
                           RCMAX = MAX( RCMAX, C( J ) );
                        } // 170
                        RATIO = RCMIN / RCMAX;
                        RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ( CCOND-RATIO ) / RATIO ) );

                        RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ( NORM-POW( N+M+1 ) ) / POW( N+M+1 ) ) );
                        for (I = 1; I <= M; I++) { // 190
                           RCMAX = ZERO;
                           for (J = 1; J <= N; J++) { // 180
                              if ( I <= J+KL && I >= J-KU ) {
                                 RATIO = ABS( R( I )*POW( I+J+1 )* C( J ) );
                                 RCMAX = MAX( RCMAX, RATIO );
                              }
                           } // 180
                           RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ONE-RCMAX ) );
                        } // 190

                        for (J = 1; J <= N; J++) { // 210
                           RCMAX = ZERO;
                           for (I = 1; I <= M; I++) { // 200
                              if ( I <= J+KL && I >= J-KU ) {
                                 RATIO = ABS( R( I )*POW( I+J+1 )* C( J ) );
                                 RCMAX = MAX( RCMAX, RATIO );
                              }
                           } // 200
                           RESLTS( 2 ) = MAX( RESLTS( 2 ), ABS( ONE-RCMAX ) );
                        } // 210
                     }
                  }

               } // 220
            } // 230
         } // 240
      } // 250
      RESLTS( 2 ) = RESLTS( 2 ) / EPS;

      // Test DPOEQU

      for (N = 0; N <= NSZ; N++) { // 290

         for (I = 1; I <= NSZ; I++) { // 270
            for (J = 1; J <= NSZ; J++) { // 260
               if ( I <= N && J == I ) {
                  A( I, J ) = POW( I+J+1 )*( -1 )**( I+J );
               } else {
                  A( I, J ) = ZERO;
               }
            } // 260
         } // 270

         dpoequ(N, A, NSZ, R, RCOND, NORM, INFO );

         if ( INFO != 0 ) {
            RESLTS( 3 ) = ONE;
         } else {
            if ( N != 0 ) {
               RESLTS( 3 ) = MAX( RESLTS( 3 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                RESLTS( 3 ) = MAX( RESLTS( 3 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) );
               for (I = 1; I <= N; I++) { // 280
                  RESLTS( 3 ) = MAX( RESLTS( 3 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ 1 ) ) );
               } // 280
            }
         }
      } // 290
      A( MAX( NSZ-1, 1 ), MAX( NSZ-1, 1 ) ) = -ONE;
      dpoequ(NSZ, A, NSZ, R, RCOND, NORM, INFO );
      if( INFO != MAX( NSZ-1, 1 ) ) RESLTS( 3 ) = ONE;
      RESLTS( 3 ) = RESLTS( 3 ) / EPS;

      // Test DPPEQU

      for (N = 0; N <= NSZ; N++) { // 360

         // Upper triangular packed storage

         for (I = 1; I <= ( N*( N+1 ) ) / 2; I++) { // 300
            AP( I ) = ZERO;
         } // 300
         for (I = 1; I <= N; I++) { // 310
            AP( ( I*( I+1 ) ) / 2 ) = POW( 2*I+1 );
         } // 310

         dppequ('U', N, AP, R, RCOND, NORM, INFO );

         if ( INFO != 0 ) {
            RESLTS( 4 ) = ONE;
         } else {
            if ( N != 0 ) {
               RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) );
               for (I = 1; I <= N; I++) { // 320
                  RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ 1 ) ) );
               } // 320
            }
         }

         // Lower triangular packed storage

         for (I = 1; I <= ( N*( N+1 ) ) / 2; I++) { // 330
            AP( I ) = ZERO;
         } // 330
         J = 1;
         for (I = 1; I <= N; I++) { // 340
            AP( J ) = POW( 2*I+1 );
            J = J + ( N-I+1 );
         } // 340

         dppequ('L', N, AP, R, RCOND, NORM, INFO );

         if ( INFO != 0 ) {
            RESLTS( 4 ) = ONE;
         } else {
            if ( N != 0 ) {
               RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) );
               for (I = 1; I <= N; I++) { // 350
                  RESLTS( 4 ) = MAX( RESLTS( 4 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+ 1 ) ) );
               } // 350
            }
         }

      } // 360
      I = ( NSZ*( NSZ+1 ) ) / 2 - 2;
      AP( I ) = -ONE;
      dppequ('L', NSZ, AP, R, RCOND, NORM, INFO );
      if( INFO != MAX( NSZ-1, 1 ) ) RESLTS( 4 ) = ONE;
      RESLTS( 4 ) = RESLTS( 4 ) / EPS;

      // Test DPBEQU

      for (N = 0; N <= NSZ; N++) { // 460
         DO 450 KL = 0, MAX( N-1, 0 );

            // Test upper triangular storage

            for (J = 1; J <= NSZ; J++) { // 380
               for (I = 1; I <= NSZB; I++) { // 370
                  AB( I, J ) = ZERO;
               } // 370
            } // 380
            for (J = 1; J <= N; J++) { // 390
               AB( KL+1, J ) = POW( 2*J+1 );
            } // 390

            dpbequ('U', N, KL, AB, NSZB, R, RCOND, NORM, INFO );

            if ( INFO != 0 ) {
               RESLTS( 5 ) = ONE;
            } else {
               if ( N != 0 ) {
                  RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                   RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) );
                  for (I = 1; I <= N; I++) { // 400
                     RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+1 ) ) );
                  } // 400
               }
            }
            if ( N != 0 ) {
               AB( KL+1, MAX( N-1, 1 ) ) = -ONE;
               dpbequ('U', N, KL, AB, NSZB, R, RCOND, NORM, INFO );
               if( INFO != MAX( N-1, 1 ) ) RESLTS( 5 ) = ONE;
            }

            // Test lower triangular storage

            for (J = 1; J <= NSZ; J++) { // 420
               for (I = 1; I <= NSZB; I++) { // 410
                  AB( I, J ) = ZERO;
               } // 410
            } // 420
            for (J = 1; J <= N; J++) { // 430
               AB( 1, J ) = POW( 2*J+1 );
            } // 430

            dpbequ('L', N, KL, AB, NSZB, R, RCOND, NORM, INFO );

            if ( INFO != 0 ) {
               RESLTS( 5 ) = ONE;
            } else {
               if ( N != 0 ) {
                  RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( RCOND-RPOW( N ) ) / RPOW( N ) ) )                   RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( NORM-POW( 2*N+1 ) ) / POW( 2*N+ 1 ) ) );
                  for (I = 1; I <= N; I++) { // 440
                     RESLTS( 5 ) = MAX( RESLTS( 5 ), ABS( ( R( I )-RPOW( I+1 ) ) / RPOW( I+1 ) ) );
                  } // 440
               }
            }
            if ( N != 0 ) {
               AB( 1, MAX( N-1, 1 ) ) = -ONE;
               dpbequ('L', N, KL, AB, NSZB, R, RCOND, NORM, INFO );
               if( INFO != MAX( N-1, 1 ) ) RESLTS( 5 ) = ONE;
            }
         } // 450
      } // 460
      RESLTS( 5 ) = RESLTS( 5 ) / EPS;
      OK = ( RESLTS( 1 ) <= THRESH ) && ( RESLTS( 2 ) <= THRESH ) && ( RESLTS( 3 ) <= THRESH ) && ( RESLTS( 4 ) <= THRESH ) && ( RESLTS( 5 ) <= THRESH );
      WRITE( NOUT, FMT = * );
      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH;
      } else {
         if( RESLTS( 1 ) > THRESH ) WRITE( NOUT, FMT = 9998 )RESLTS( 1 ), THRESH;
         if( RESLTS( 2 ) > THRESH ) WRITE( NOUT, FMT = 9997 )RESLTS( 2 ), THRESH;
         if( RESLTS( 3 ) > THRESH ) WRITE( NOUT, FMT = 9996 )RESLTS( 3 ), THRESH;
         if( RESLTS( 4 ) > THRESH ) WRITE( NOUT, FMT = 9995 )RESLTS( 4 ), THRESH;
         IF( RESLTS( 5 ) > THRESH ) WRITE( NOUT, FMT = 9994 )RESLTS( 5 ), THRESH;
      }
 9999 FORMAT( 1X, 'All tests for ', A3, ' routines passed the threshold' );
 9998 FORMAT( ' DGEEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 );
 9997 FORMAT( ' DGBEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 );
 9996 FORMAT( ' DPOEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 );
 9995 FORMAT( ' DPPEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 );
 9994 FORMAT( ' DPBEQU failed test with value ', D10.3, ' exceeding', ' threshold ', D10.3 );
      return;

      // End of DCHKEQ

      }
