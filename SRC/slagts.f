      SUBROUTINE SLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, JOB, N;
      REAL               TOL
      // ..
      // .. Array Arguments ..
      int                IN( * );
      REAL               A( * ), B( * ), C( * ), D( * ), Y( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                K;
      REAL               ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0
      if ( ( ABS( JOB ).GT.2 ) || ( JOB == 0 ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      }
      if ( INFO != 0 ) {
         xerbla('SLAGTS', -INFO );
         RETURN
      }

      if (N == 0) RETURN;

      EPS = SLAMCH( 'Epsilon' )
      SFMIN = SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SFMIN

      if ( JOB < 0 ) {
         if ( TOL.LE.ZERO ) {
            TOL = ABS( A( 1 ) )
            if (N.GT.1) TOL = MAX( TOL, ABS( A( 2 ) ), ABS( B( 1 ) ) );
            for (K = 3; K <= N; K++) { // 10
               TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ), ABS( D( K-2 ) ) )
            } // 10
            TOL = TOL*EPS
            if (TOL == ZERO) TOL = EPS;
         }
      }

      if ( ABS( JOB ) == 1 ) {
         for (K = 2; K <= N; K++) { // 20
            if ( IN( K-1 ) == 0 ) {
               Y( K ) = Y( K ) - C( K-1 )*Y( K-1 )
            } else {
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            }
         } // 20
         if ( JOB == 1 ) {
            DO 30 K = N, 1, -1
               if ( K.LE.N-2 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               } else if ( K == N-1 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               } else {
                  TEMP = Y( K )
               }
               AK = A( K )
               ABSAK = ABS( AK )
               if ( ABSAK < ONE ) {
                  if ( ABSAK < SFMIN ) {
                     if ( ABSAK == ZERO || ABS( TEMP )*SFMIN.GT.ABSAK ) {
                        INFO = K
                        RETURN
                     } else {
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     }
                  } else if ( ABS( TEMP ).GT.ABSAK*BIGNUM ) {
                     INFO = K
                     RETURN
                  }
               }
               Y( K ) = TEMP / AK
            } // 30
         } else {
            DO 50 K = N, 1, -1
               if ( K.LE.N-2 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               } else if ( K == N-1 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               } else {
                  TEMP = Y( K )
               }
               AK = A( K )
               PERT = SIGN( TOL, AK )
               } // 40
               ABSAK = ABS( AK )
               if ( ABSAK < ONE ) {
                  if ( ABSAK < SFMIN ) {
                     if ( ABSAK == ZERO || ABS( TEMP )*SFMIN.GT.ABSAK ) {
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 40
                     } else {
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     }
                  } else if ( ABS( TEMP ).GT.ABSAK*BIGNUM ) {
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 40
                  }
               }
               Y( K ) = TEMP / AK
            } // 50
         }
      } else {

         // Come to here if  JOB = 2 or -2

         if ( JOB == 2 ) {
            for (K = 1; K <= N; K++) { // 60
               if ( K.GE.3 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               } else if ( K == 2 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               } else {
                  TEMP = Y( K )
               }
               AK = A( K )
               ABSAK = ABS( AK )
               if ( ABSAK < ONE ) {
                  if ( ABSAK < SFMIN ) {
                     if ( ABSAK == ZERO || ABS( TEMP )*SFMIN.GT.ABSAK ) {
                        INFO = K
                        RETURN
                     } else {
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     }
                  } else if ( ABS( TEMP ).GT.ABSAK*BIGNUM ) {
                     INFO = K
                     RETURN
                  }
               }
               Y( K ) = TEMP / AK
            } // 60
         } else {
            for (K = 1; K <= N; K++) { // 80
               if ( K.GE.3 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               } else if ( K == 2 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               } else {
                  TEMP = Y( K )
               }
               AK = A( K )
               PERT = SIGN( TOL, AK )
               } // 70
               ABSAK = ABS( AK )
               if ( ABSAK < ONE ) {
                  if ( ABSAK < SFMIN ) {
                     if ( ABSAK == ZERO || ABS( TEMP )*SFMIN.GT.ABSAK ) {
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 70
                     } else {
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     }
                  } else if ( ABS( TEMP ).GT.ABSAK*BIGNUM ) {
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 70
                  }
               }
               Y( K ) = TEMP / AK
            } // 80
         }

         DO 90 K = N, 2, -1
            if ( IN( K-1 ) == 0 ) {
               Y( K-1 ) = Y( K-1 ) - C( K-1 )*Y( K )
            } else {
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            }
         } // 90
      }

      // End of SLAGTS

      }
