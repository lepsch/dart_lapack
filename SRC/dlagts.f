      SUBROUTINE DLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, JOB, N;
      double             TOL;
      // ..
      // .. Array Arguments ..
      int                IN( * );
      double             A( * ), B( * ), C( * ), D( * ), Y( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                K;
      double             ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SIGN
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0
      if ( ( ABS( JOB ).GT.2 ) .OR. ( JOB.EQ.0 ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         xerbla('DLAGTS', -INFO );
         RETURN
      }

      IF( N.EQ.0 ) RETURN

      EPS = DLAMCH( 'Epsilon' )
      SFMIN = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SFMIN

      if ( JOB.LT.0 ) {
         if ( TOL.LE.ZERO ) {
            TOL = ABS( A( 1 ) )
            IF( N.GT.1 ) TOL = MAX( TOL, ABS( A( 2 ) ), ABS( B( 1 ) ) )
            DO 10 K = 3, N
               TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ), ABS( D( K-2 ) ) )
   10       CONTINUE
            TOL = TOL*EPS
            IF( TOL.EQ.ZERO ) TOL = EPS
         }
      }

      if ( ABS( JOB ).EQ.1 ) {
         DO 20 K = 2, N
            if ( IN( K-1 ).EQ.0 ) {
               Y( K ) = Y( K ) - C( K-1 )*Y( K-1 )
            } else {
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            }
   20    CONTINUE
         if ( JOB.EQ.1 ) {
            DO 30 K = N, 1, -1
               if ( K.LE.N-2 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               } else if ( K.EQ.N-1 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               } else {
                  TEMP = Y( K )
               }
               AK = A( K )
               ABSAK = ABS( AK )
               if ( ABSAK.LT.ONE ) {
                  if ( ABSAK.LT.SFMIN ) {
                     if ( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) {
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
   30       CONTINUE
         } else {
            DO 50 K = N, 1, -1
               if ( K.LE.N-2 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               } else if ( K.EQ.N-1 ) {
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               } else {
                  TEMP = Y( K )
               }
               AK = A( K )
               PERT = SIGN( TOL, AK )
   40          CONTINUE
               ABSAK = ABS( AK )
               if ( ABSAK.LT.ONE ) {
                  if ( ABSAK.LT.SFMIN ) {
                     if ( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) {
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
   50       CONTINUE
         }
      } else {

         // Come to here if  JOB = 2 or -2

         if ( JOB.EQ.2 ) {
            DO 60 K = 1, N
               if ( K.GE.3 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               } else if ( K.EQ.2 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               } else {
                  TEMP = Y( K )
               }
               AK = A( K )
               ABSAK = ABS( AK )
               if ( ABSAK.LT.ONE ) {
                  if ( ABSAK.LT.SFMIN ) {
                     if ( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) {
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
   60       CONTINUE
         } else {
            DO 80 K = 1, N
               if ( K.GE.3 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               } else if ( K.EQ.2 ) {
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               } else {
                  TEMP = Y( K )
               }
               AK = A( K )
               PERT = SIGN( TOL, AK )
   70          CONTINUE
               ABSAK = ABS( AK )
               if ( ABSAK.LT.ONE ) {
                  if ( ABSAK.LT.SFMIN ) {
                     if ( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK ) {
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
   80       CONTINUE
         }

         DO 90 K = N, 2, -1
            if ( IN( K-1 ).EQ.0 ) {
               Y( K-1 ) = Y( K-1 ) - C( K-1 )*Y( K )
            } else {
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            }
   90    CONTINUE
      }

      // End of DLAGTS

      }
