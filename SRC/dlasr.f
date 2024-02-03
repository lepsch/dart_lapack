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
      if ( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) {
         INFO = 1
      } else if ( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, 'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) {
         INFO = 2
      } else if ( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) {
         INFO = 3
      } else if ( M.LT.0 ) {
         INFO = 4
      } else if ( N.LT.0 ) {
         INFO = 5
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = 9
      }
      if ( INFO.NE.0 ) {
         xerbla('DLASR ', INFO );
         RETURN
      }

      // Quick return if possible

      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN
      if ( LSAME( SIDE, 'L' ) ) {

         // Form  P * A

         if ( LSAME( PIVOT, 'V' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  }
   20          CONTINUE
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  }
   40          CONTINUE
            }
         } else if ( LSAME( PIVOT, 'T' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  }
   60          CONTINUE
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  }
   80          CONTINUE
            }
         } else if ( LSAME( PIVOT, 'B' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  }
  100          CONTINUE
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  }
  120          CONTINUE
            }
         }
      } else if ( LSAME( SIDE, 'R' ) ) {

         // Form A * P**T

         if ( LSAME( PIVOT, 'V' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  }
  140          CONTINUE
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  }
  160          CONTINUE
            }
         } else if ( LSAME( PIVOT, 'T' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  }
  180          CONTINUE
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  }
  200          CONTINUE
            }
         } else if ( LSAME( PIVOT, 'B' ) ) {
            if ( LSAME( DIRECT, 'F' ) ) {
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  }
  220          CONTINUE
            } else if ( LSAME( DIRECT, 'B' ) ) {
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  if ( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) {
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  }
  240          CONTINUE
            }
         }
      }

      RETURN

      // End of DLASR

      }
