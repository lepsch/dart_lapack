      SUBROUTINE DTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDT, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      double             T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLV, BOTHV, LEFTV, OVER, PAIR, RIGHTV, SOMEV;
      int                I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, N2;
      double             BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL, VCRIT, VMAX, WI, WR, XNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DDOT, DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DDOT, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMV, DLALN2, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Local Arrays ..
      double             X( 2, 2 );
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV

      ALLV = LSAME( HOWMNY, 'A' )
      OVER = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )

      INFO = 0
      if ( .NOT.RIGHTV .AND. .NOT.LEFTV ) {
         INFO = -1
      } else if ( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDT.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) {
         INFO = -8
      } else if ( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) {
         INFO = -10
      } else {

         // Set M to the number of columns required to store the selected
         // eigenvectors, standardize the array SELECT if necessary, and
         // test MM.

         if ( SOMEV ) {
            M = 0
            PAIR = .FALSE.
            for (J = 1; J <= N; J++) { // 10
               if ( PAIR ) {
                  PAIR = .FALSE.
                  SELECT( J ) = .FALSE.
               } else {
                  if ( J.LT.N ) {
                     if ( T( J+1, J ).EQ.ZERO ) {
                        IF( SELECT( J ) ) M = M + 1
                     } else {
                        PAIR = .TRUE.
                        if ( SELECT( J ) .OR. SELECT( J+1 ) ) {
                           SELECT( J ) = .TRUE.
                           M = M + 2
                        }
                     }
                  } else {
                     IF( SELECT( N ) ) M = M + 1
                  }
               }
   10       CONTINUE
         } else {
            M = N
         }

         if ( MM.LT.M ) {
            INFO = -11
         }
      }
      if ( INFO.NE.0 ) {
         xerbla('DTREVC', -INFO );
         RETURN
      }

      // Quick return if possible.

      IF( N.EQ.0 ) RETURN

      // Set the constants to control overflow.

      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
      BIGNUM = ( ONE-ULP ) / SMLNUM

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      WORK( 1 ) = ZERO
      for (J = 2; J <= N; J++) { // 30
         WORK( J ) = ZERO
         DO 20 I = 1, J - 1
            WORK( J ) = WORK( J ) + ABS( T( I, J ) )
   20    CONTINUE
   30 CONTINUE

      // Index IP is used to specify the real or complex eigenvalue:
        // IP = 0, real eigenvalue,
             // 1, first of conjugate complex pair: (wr,wi)
            // -1, second of conjugate complex pair: (wr,wi)

      N2 = 2*N

      if ( RIGHTV ) {

         // Compute right eigenvectors.

         IP = 0
         IS = M
         DO 140 KI = N, 1, -1

            IF( IP.EQ.1 ) GO TO 130             IF( KI.EQ.1 ) GO TO 40             IF( T( KI, KI-1 ).EQ.ZERO ) GO TO 40
            IP = -1

   40       CONTINUE
            if ( SOMEV ) {
               if ( IP.EQ.0 ) {
                  IF( .NOT.SELECT( KI ) ) GO TO 130
               } else {
                  IF( .NOT.SELECT( KI-1 ) ) GO TO 130
               }
            }

            // Compute the KI-th eigenvalue (WR,WI).

            WR = T( KI, KI )
            WI = ZERO
            IF( IP.NE.0 ) WI = SQRT( ABS( T( KI, KI-1 ) ) )* SQRT( ABS( T( KI-1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )

            if ( IP.EQ.0 ) {

               // Real right eigenvector

               WORK( KI+N ) = ONE

               // Form right-hand side

               DO 50 K = 1, KI - 1
                  WORK( K+N ) = -T( K, KI )
   50          CONTINUE

               // Solve the upper quasi-triangular system:
                  // (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.

               JNXT = KI - 1
               DO 60 J = KI - 1, 1, -1
                  IF( J.GT.JNXT ) GO TO 60
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  if ( J.GT.1 ) {
                     if ( T( J, J-1 ).NE.ZERO ) {
                        J1 = J - 1
                        JNXT = J - 2
                     }
                  }

                  if ( J1.EQ.J2 ) {

                     // 1-by-1 diagonal block

                     dlaln2(.FALSE., 1, 1, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+N ), N, WR, ZERO, X, 2, SCALE, XNORM, IERR );

                     // Scale X(1,1) to avoid overflow when updating
                     // the right-hand side.

                     if ( XNORM.GT.ONE ) {
                        if ( WORK( J ).GT.BIGNUM / XNORM ) {
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        }
                     }

                     // Scale if necessary

                     IF( SCALE.NE.ONE ) CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )

                     // Update right-hand side

                     daxpy(J-1, -X( 1, 1 ), T( 1, J ), 1, WORK( 1+N ), 1 );

                  } else {

                     // 2-by-2 diagonal block

                     dlaln2(.FALSE., 2, 1, SMIN, ONE, T( J-1, J-1 ), LDT, ONE, ONE, WORK( J-1+N ), N, WR, ZERO, X, 2, SCALE, XNORM, IERR );

                     // Scale X(1,1) and X(2,1) to avoid overflow when
                     // updating the right-hand side.

                     if ( XNORM.GT.ONE ) {
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        if ( BETA.GT.BIGNUM / XNORM ) {
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 2, 1 ) = X( 2, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        }
                     }

                     // Scale if necessary

                     IF( SCALE.NE.ONE ) CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )

                     // Update right-hand side

                     daxpy(J-2, -X( 1, 1 ), T( 1, J-1 ), 1, WORK( 1+N ), 1 )                      CALL DAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, WORK( 1+N ), 1 );
                  }
   60          CONTINUE

               // Copy the vector x or Q*x to VR and normalize.

               if ( .NOT.OVER ) {
                  dcopy(KI, WORK( 1+N ), 1, VR( 1, IS ), 1 );

                  II = IDAMAX( KI, VR( 1, IS ), 1 )
                  REMAX = ONE / ABS( VR( II, IS ) )
                  dscal(KI, REMAX, VR( 1, IS ), 1 );

                  DO 70 K = KI + 1, N
                     VR( K, IS ) = ZERO
   70             CONTINUE
               } else {
                  IF( KI.GT.1 ) CALL DGEMV( 'N', N, KI-1, ONE, VR, LDVR, WORK( 1+N ), 1, WORK( KI+N ), VR( 1, KI ), 1 )

                  II = IDAMAX( N, VR( 1, KI ), 1 )
                  REMAX = ONE / ABS( VR( II, KI ) )
                  dscal(N, REMAX, VR( 1, KI ), 1 );
               }

            } else {

               // Complex right eigenvector.

               // Initial solve
                 // [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
                 // [ (T(KI,KI-1)   T(KI,KI)   )               ]

               if ( ABS( T( KI-1, KI ) ).GE.ABS( T( KI, KI-1 ) ) ) {
                  WORK( KI-1+N ) = ONE
                  WORK( KI+N2 ) = WI / T( KI-1, KI )
               } else {
                  WORK( KI-1+N ) = -WI / T( KI, KI-1 )
                  WORK( KI+N2 ) = ONE
               }
               WORK( KI+N ) = ZERO
               WORK( KI-1+N2 ) = ZERO

               // Form right-hand side

               DO 80 K = 1, KI - 2
                  WORK( K+N ) = -WORK( KI-1+N )*T( K, KI-1 )
                  WORK( K+N2 ) = -WORK( KI+N2 )*T( K, KI )
   80          CONTINUE

               // Solve upper quasi-triangular system:
               // (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)

               JNXT = KI - 2
               DO 90 J = KI - 2, 1, -1
                  IF( J.GT.JNXT ) GO TO 90
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  if ( J.GT.1 ) {
                     if ( T( J, J-1 ).NE.ZERO ) {
                        J1 = J - 1
                        JNXT = J - 2
                     }
                  }

                  if ( J1.EQ.J2 ) {

                     // 1-by-1 diagonal block

                     dlaln2(.FALSE., 1, 2, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+N ), N, WR, WI, X, 2, SCALE, XNORM, IERR );

                     // Scale X(1,1) and X(1,2) to avoid overflow when
                     // updating the right-hand side.

                     if ( XNORM.GT.ONE ) {
                        if ( WORK( J ).GT.BIGNUM / XNORM ) {
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 1, 2 ) = X( 1, 2 ) / XNORM
                           SCALE = SCALE / XNORM
                        }
                     }

                     // Scale if necessary

                     if ( SCALE.NE.ONE ) {
                        dscal(KI, SCALE, WORK( 1+N ), 1 );
                        dscal(KI, SCALE, WORK( 1+N2 ), 1 );
                     }
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )

                     // Update the right-hand side

                     daxpy(J-1, -X( 1, 1 ), T( 1, J ), 1, WORK( 1+N ), 1 )                      CALL DAXPY( J-1, -X( 1, 2 ), T( 1, J ), 1, WORK( 1+N2 ), 1 );

                  } else {

                     // 2-by-2 diagonal block

                     dlaln2(.FALSE., 2, 2, SMIN, ONE, T( J-1, J-1 ), LDT, ONE, ONE, WORK( J-1+N ), N, WR, WI, X, 2, SCALE, XNORM, IERR );

                     // Scale X to avoid overflow when updating
                     // the right-hand side.

                     if ( XNORM.GT.ONE ) {
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        if ( BETA.GT.BIGNUM / XNORM ) {
                           REC = ONE / XNORM
                           X( 1, 1 ) = X( 1, 1 )*REC
                           X( 1, 2 ) = X( 1, 2 )*REC
                           X( 2, 1 ) = X( 2, 1 )*REC
                           X( 2, 2 ) = X( 2, 2 )*REC
                           SCALE = SCALE*REC
                        }
                     }

                     // Scale if necessary

                     if ( SCALE.NE.ONE ) {
                        dscal(KI, SCALE, WORK( 1+N ), 1 );
                        dscal(KI, SCALE, WORK( 1+N2 ), 1 );
                     }
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
                     WORK( J-1+N2 ) = X( 1, 2 )
                     WORK( J+N2 ) = X( 2, 2 )

                     // Update the right-hand side

                     daxpy(J-2, -X( 1, 1 ), T( 1, J-1 ), 1, WORK( 1+N ), 1 )                      CALL DAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, WORK( 1+N ), 1 )                      CALL DAXPY( J-2, -X( 1, 2 ), T( 1, J-1 ), 1, WORK( 1+N2 ), 1 )                      CALL DAXPY( J-2, -X( 2, 2 ), T( 1, J ), 1, WORK( 1+N2 ), 1 );
                  }
   90          CONTINUE

               // Copy the vector x or Q*x to VR and normalize.

               if ( .NOT.OVER ) {
                  dcopy(KI, WORK( 1+N ), 1, VR( 1, IS-1 ), 1 );
                  dcopy(KI, WORK( 1+N2 ), 1, VR( 1, IS ), 1 );

                  EMAX = ZERO
                  for (K = 1; K <= KI; K++) { // 100
                     EMAX = MAX( EMAX, ABS( VR( K, IS-1 ) )+ ABS( VR( K, IS ) ) )
  100             CONTINUE

                  REMAX = ONE / EMAX
                  dscal(KI, REMAX, VR( 1, IS-1 ), 1 );
                  dscal(KI, REMAX, VR( 1, IS ), 1 );

                  DO 110 K = KI + 1, N
                     VR( K, IS-1 ) = ZERO
                     VR( K, IS ) = ZERO
  110             CONTINUE

               } else {

                  if ( KI.GT.2 ) {
                     dgemv('N', N, KI-2, ONE, VR, LDVR, WORK( 1+N ), 1, WORK( KI-1+N ), VR( 1, KI-1 ), 1 )                      CALL DGEMV( 'N', N, KI-2, ONE, VR, LDVR, WORK( 1+N2 ), 1, WORK( KI+N2 ), VR( 1, KI ), 1 );
                  } else {
                     dscal(N, WORK( KI-1+N ), VR( 1, KI-1 ), 1 );
                     dscal(N, WORK( KI+N2 ), VR( 1, KI ), 1 );
                  }

                  EMAX = ZERO
                  for (K = 1; K <= N; K++) { // 120
                     EMAX = MAX( EMAX, ABS( VR( K, KI-1 ) )+ ABS( VR( K, KI ) ) )
  120             CONTINUE
                  REMAX = ONE / EMAX
                  dscal(N, REMAX, VR( 1, KI-1 ), 1 );
                  dscal(N, REMAX, VR( 1, KI ), 1 );
               }
            }

            IS = IS - 1
            IF( IP.NE.0 ) IS = IS - 1
  130       CONTINUE
            IF( IP.EQ.1 ) IP = 0             IF( IP.EQ.-1 ) IP = 1
  140    CONTINUE
      }

      if ( LEFTV ) {

         // Compute left eigenvectors.

         IP = 0
         IS = 1
         for (KI = 1; KI <= N; KI++) { // 260

            IF( IP.EQ.-1 ) GO TO 250             IF( KI.EQ.N ) GO TO 150             IF( T( KI+1, KI ).EQ.ZERO ) GO TO 150
            IP = 1

  150       CONTINUE
            if ( SOMEV ) {
               IF( .NOT.SELECT( KI ) ) GO TO 250
            }

            // Compute the KI-th eigenvalue (WR,WI).

            WR = T( KI, KI )
            WI = ZERO
            IF( IP.NE.0 ) WI = SQRT( ABS( T( KI, KI+1 ) ) )* SQRT( ABS( T( KI+1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )

            if ( IP.EQ.0 ) {

               // Real left eigenvector.

               WORK( KI+N ) = ONE

               // Form right-hand side

               DO 160 K = KI + 1, N
                  WORK( K+N ) = -T( KI, K )
  160          CONTINUE

               // Solve the quasi-triangular system:
                  // (T(KI+1:N,KI+1:N) - WR)**T*X = SCALE*WORK

               VMAX = ONE
               VCRIT = BIGNUM

               JNXT = KI + 1
               DO 170 J = KI + 1, N
                  IF( J.LT.JNXT ) GO TO 170
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  if ( J.LT.N ) {
                     if ( T( J+1, J ).NE.ZERO ) {
                        J2 = J + 1
                        JNXT = J + 2
                     }
                  }

                  if ( J1.EQ.J2 ) {

                     // 1-by-1 diagonal block

                     // Scale if necessary to avoid overflow when forming
                     // the right-hand side.

                     if ( WORK( J ).GT.VCRIT ) {
                        REC = ONE / VMAX
                        dscal(N-KI+1, REC, WORK( KI+N ), 1 );
                        VMAX = ONE
                        VCRIT = BIGNUM
                     }

                     WORK( J+N ) = WORK( J+N ) - DDOT( J-KI-1, T( KI+1, J ), 1, WORK( KI+1+N ), 1 )

                     // Solve (T(J,J)-WR)**T*X = WORK

                     dlaln2(.FALSE., 1, 1, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+N ), N, WR, ZERO, X, 2, SCALE, XNORM, IERR );

                     // Scale if necessary

                     IF( SCALE.NE.ONE ) CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     VMAX = MAX( ABS( WORK( J+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX

                  } else {

                     // 2-by-2 diagonal block

                     // Scale if necessary to avoid overflow when forming
                     // the right-hand side.

                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     if ( BETA.GT.VCRIT ) {
                        REC = ONE / VMAX
                        dscal(N-KI+1, REC, WORK( KI+N ), 1 );
                        VMAX = ONE
                        VCRIT = BIGNUM
                     }

                     WORK( J+N ) = WORK( J+N ) - DDOT( J-KI-1, T( KI+1, J ), 1, WORK( KI+1+N ), 1 )

                     WORK( J+1+N ) = WORK( J+1+N ) - DDOT( J-KI-1, T( KI+1, J+1 ), 1, WORK( KI+1+N ), 1 )

                     // Solve
                       // [T(J,J)-WR   T(J,J+1)     ]**T * X = SCALE*( WORK1 )
                       // [T(J+1,J)    T(J+1,J+1)-WR]                ( WORK2 )

                     dlaln2(.TRUE., 2, 1, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+N ), N, WR, ZERO, X, 2, SCALE, XNORM, IERR );

                     // Scale if necessary

                     IF( SCALE.NE.ONE ) CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+1+N ) = X( 2, 1 )

                     VMAX = MAX( ABS( WORK( J+N ) ), ABS( WORK( J+1+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX

                  }
  170          CONTINUE

               // Copy the vector x or Q*x to VL and normalize.

               if ( .NOT.OVER ) {
                  dcopy(N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 );

                  II = IDAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
                  REMAX = ONE / ABS( VL( II, IS ) )
                  dscal(N-KI+1, REMAX, VL( KI, IS ), 1 );

                  DO 180 K = 1, KI - 1
                     VL( K, IS ) = ZERO
  180             CONTINUE

               } else {

                  IF( KI.LT.N ) CALL DGEMV( 'N', N, N-KI, ONE, VL( 1, KI+1 ), LDVL, WORK( KI+1+N ), 1, WORK( KI+N ), VL( 1, KI ), 1 )

                  II = IDAMAX( N, VL( 1, KI ), 1 )
                  REMAX = ONE / ABS( VL( II, KI ) )
                  dscal(N, REMAX, VL( 1, KI ), 1 );

               }

            } else {

               // Complex left eigenvector.

                // Initial solve:
                  // ((T(KI,KI)    T(KI,KI+1) )**T - (WR - I* WI))*X = 0.
                  // ((T(KI+1,KI) T(KI+1,KI+1))                )

               if ( ABS( T( KI, KI+1 ) ).GE.ABS( T( KI+1, KI ) ) ) {
                  WORK( KI+N ) = WI / T( KI, KI+1 )
                  WORK( KI+1+N2 ) = ONE
               } else {
                  WORK( KI+N ) = ONE
                  WORK( KI+1+N2 ) = -WI / T( KI+1, KI )
               }
               WORK( KI+1+N ) = ZERO
               WORK( KI+N2 ) = ZERO

               // Form right-hand side

               DO 190 K = KI + 2, N
                  WORK( K+N ) = -WORK( KI+N )*T( KI, K )
                  WORK( K+N2 ) = -WORK( KI+1+N2 )*T( KI+1, K )
  190          CONTINUE

               // Solve complex quasi-triangular system:
               // ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2

               VMAX = ONE
               VCRIT = BIGNUM

               JNXT = KI + 2
               DO 200 J = KI + 2, N
                  IF( J.LT.JNXT ) GO TO 200
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  if ( J.LT.N ) {
                     if ( T( J+1, J ).NE.ZERO ) {
                        J2 = J + 1
                        JNXT = J + 2
                     }
                  }

                  if ( J1.EQ.J2 ) {

                     // 1-by-1 diagonal block

                     // Scale if necessary to avoid overflow when
                     // forming the right-hand side elements.

                     if ( WORK( J ).GT.VCRIT ) {
                        REC = ONE / VMAX
                        dscal(N-KI+1, REC, WORK( KI+N ), 1 );
                        dscal(N-KI+1, REC, WORK( KI+N2 ), 1 );
                        VMAX = ONE
                        VCRIT = BIGNUM
                     }

                     WORK( J+N ) = WORK( J+N ) - DDOT( J-KI-2, T( KI+2, J ), 1, WORK( KI+2+N ), 1 )                      WORK( J+N2 ) = WORK( J+N2 ) - DDOT( J-KI-2, T( KI+2, J ), 1, WORK( KI+2+N2 ), 1 )

                     // Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2

                     dlaln2(.FALSE., 1, 2, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+N ), N, WR, -WI, X, 2, SCALE, XNORM, IERR );

                     // Scale if necessary

                     if ( SCALE.NE.ONE ) {
                        dscal(N-KI+1, SCALE, WORK( KI+N ), 1 );
                        dscal(N-KI+1, SCALE, WORK( KI+N2 ), 1 );
                     }
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     VMAX = MAX( ABS( WORK( J+N ) ), ABS( WORK( J+N2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX

                  } else {

                     // 2-by-2 diagonal block

                     // Scale if necessary to avoid overflow when forming
                     // the right-hand side elements.

                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     if ( BETA.GT.VCRIT ) {
                        REC = ONE / VMAX
                        dscal(N-KI+1, REC, WORK( KI+N ), 1 );
                        dscal(N-KI+1, REC, WORK( KI+N2 ), 1 );
                        VMAX = ONE
                        VCRIT = BIGNUM
                     }

                     WORK( J+N ) = WORK( J+N ) - DDOT( J-KI-2, T( KI+2, J ), 1, WORK( KI+2+N ), 1 )

                     WORK( J+N2 ) = WORK( J+N2 ) - DDOT( J-KI-2, T( KI+2, J ), 1, WORK( KI+2+N2 ), 1 )

                     WORK( J+1+N ) = WORK( J+1+N ) - DDOT( J-KI-2, T( KI+2, J+1 ), 1, WORK( KI+2+N ), 1 )

                     WORK( J+1+N2 ) = WORK( J+1+N2 ) - DDOT( J-KI-2, T( KI+2, J+1 ), 1, WORK( KI+2+N2 ), 1 )

                     // Solve 2-by-2 complex linear equation
                       // ([T(j,j)   T(j,j+1)  ]**T-(wr-i*wi)*I)*X = SCALE*B
                       // ([T(j+1,j) T(j+1,j+1)]               )

                     dlaln2(.TRUE., 2, 2, SMIN, ONE, T( J, J ), LDT, ONE, ONE, WORK( J+N ), N, WR, -WI, X, 2, SCALE, XNORM, IERR );

                     // Scale if necessary

                     if ( SCALE.NE.ONE ) {
                        dscal(N-KI+1, SCALE, WORK( KI+N ), 1 );
                        dscal(N-KI+1, SCALE, WORK( KI+N2 ), 1 );
                     }
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     WORK( J+1+N ) = X( 2, 1 )
                     WORK( J+1+N2 ) = X( 2, 2 )
                     VMAX = MAX( ABS( X( 1, 1 ) ), ABS( X( 1, 2 ) ), ABS( X( 2, 1 ) ), ABS( X( 2, 2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX

                  }
  200          CONTINUE

               // Copy the vector x or Q*x to VL and normalize.

               if ( .NOT.OVER ) {
                  dcopy(N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 );
                  dcopy(N-KI+1, WORK( KI+N2 ), 1, VL( KI, IS+1 ), 1 );

                  EMAX = ZERO
                  for (K = KI; K <= N; K++) { // 220
                     EMAX = MAX( EMAX, ABS( VL( K, IS ) )+ ABS( VL( K, IS+1 ) ) )
  220             CONTINUE
                  REMAX = ONE / EMAX
                  dscal(N-KI+1, REMAX, VL( KI, IS ), 1 );
                  dscal(N-KI+1, REMAX, VL( KI, IS+1 ), 1 );

                  DO 230 K = 1, KI - 1
                     VL( K, IS ) = ZERO
                     VL( K, IS+1 ) = ZERO
  230             CONTINUE
               } else {
                  if ( KI.LT.N-1 ) {
                     dgemv('N', N, N-KI-1, ONE, VL( 1, KI+2 ), LDVL, WORK( KI+2+N ), 1, WORK( KI+N ), VL( 1, KI ), 1 )                      CALL DGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ), LDVL, WORK( KI+2+N2 ), 1, WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 );
                  } else {
                     dscal(N, WORK( KI+N ), VL( 1, KI ), 1 );
                     dscal(N, WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 );
                  }

                  EMAX = ZERO
                  for (K = 1; K <= N; K++) { // 240
                     EMAX = MAX( EMAX, ABS( VL( K, KI ) )+ ABS( VL( K, KI+1 ) ) )
  240             CONTINUE
                  REMAX = ONE / EMAX
                  dscal(N, REMAX, VL( 1, KI ), 1 );
                  dscal(N, REMAX, VL( 1, KI+1 ), 1 );

               }

            }

            IS = IS + 1
            IF( IP.NE.0 ) IS = IS + 1
  250       CONTINUE
            IF( IP.EQ.-1 ) IP = 0             IF( IP.EQ.1 ) IP = -1

  260    CONTINUE

      }

      RETURN

      // End of DTREVC

      }
