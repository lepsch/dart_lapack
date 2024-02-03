      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIRECT, PIVOT, SIDE;
      int                LDA, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), C( * ), S( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double             CTEMP, STEMP, TEMP;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      if ( .NOT.( LSAME( SIDE, 'L' ) || LSAME( SIDE, 'R' ) ) ) {
         INFO = 1
      } else if ( .NOT.( LSAME( PIVOT, 'V' ) || LSAME( PIVOT, 'T' ) || LSAME( PIVOT, 'B' ) ) ) {
         INFO = 2
      } else if ( .NOT.( LSAME( DIRECT, 'F' ) || LSAME( DIRECT, 'B' ) ) ) {
         INFO = 3
      } else if ( M < 0 ) {
         INFO = 4
      } else if ( N < 0 ) {
         INFO = 5
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = 9
      }
      if ( INFO != 0 ) {
         xerbla('DLASR ', INFO );
         RETURN
      }

      // Quick return if possible

      IF( ( M == 0 ) || ( N == 0 ) ) RETURN
      if ( LSAME( SIDE, 'L' ) ) {

         // Form  P * A

         if ( LSAME( PIVOT, 'V' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               for (J = 1; J <= M - 1; J++) { // 20
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= N; I++) { // 10
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                     } // 10
                  }
               } // 20
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= N; I++) { // 30
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                     } // 30
                  }
               } // 40
            }
         } else if ( LSAME( PIVOT, 'T' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               for (J = 2; J <= M; J++) { // 60
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= N; I++) { // 50
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                     } // 50
                  }
               } // 60
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= N; I++) { // 70
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                     } // 70
                  }
               } // 80
            }
         } else if ( LSAME( PIVOT, 'B' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               for (J = 1; J <= M - 1; J++) { // 100
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= N; I++) { // 90
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                     } // 90
                  }
               } // 100
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= N; I++) { // 110
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                     } // 110
                  }
               } // 120
            }
         }
      } else if ( LSAME( SIDE, 'R' ) ) {

         // Form A * P**T

         if ( LSAME( PIVOT, 'V' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               for (J = 1; J <= N - 1; J++) { // 140
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= M; I++) { // 130
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                     } // 130
                  }
               } // 140
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= M; I++) { // 150
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                     } // 150
                  }
               } // 160
            }
         } else if ( LSAME( PIVOT, 'T' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               for (J = 2; J <= N; J++) { // 180
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= M; I++) { // 170
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                     } // 170
                  }
               } // 180
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= M; I++) { // 190
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                     } // 190
                  }
               } // 200
            }
         } else if ( LSAME( PIVOT, 'B' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               for (J = 1; J <= N - 1; J++) { // 220
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= M; I++) { // 210
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                     } // 210
                  }
               } // 220
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP != ONE ) || ( STEMP != ZERO ) ) {
                     for (I = 1; I <= M; I++) { // 230
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                     } // 230
                  }
               } // 240
            }
         }
      }

      RETURN

      // End of DLASR

      }
