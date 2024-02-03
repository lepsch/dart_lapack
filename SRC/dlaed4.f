      SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I, INFO, N;
      double             DLAM, RHO;
      // ..
      // .. Array Arguments ..
      double             D( * ), DELTA( * ), Z( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                MAXIT;
      const              MAXIT = 30 ;
      double             ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0, TEN = 10.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               ORGATI, SWTCH, SWTCH3;
      int                II, IIM1, IIP1, IP1, ITER, J, NITER;
      double             A, B, C, DEL, DLTLB, DLTUB, DPHI, DPSI, DW, EPS, ERRETM, ETA, MIDPT, PHI, PREW, PSI, RHOINV, TAU, TEMP, TEMP1, W;
      // ..
      // .. Local Arrays ..
      double             ZZ( 3 );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAED5, DLAED6
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Since this routine is called in an inner loop, we do no argument
      // checking.

      // Quick return for N=1 and 2.

      INFO = 0
      if ( N.EQ.1 ) {

          // Presumably, I=1 upon entry

         DLAM = D( 1 ) + RHO*Z( 1 )*Z( 1 )
         DELTA( 1 ) = ONE
         RETURN
      }
      if ( N.EQ.2 ) {
         dlaed5(I, D, Z, DELTA, RHO, DLAM );
         RETURN
      }

      // Compute machine epsilon

      EPS = DLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO

      // The case I = N

      if ( I.EQ.N ) {

         // Initialize some basic variables

         II = N - 1
         NITER = 1

         // Calculate initial guess

         MIDPT = RHO / TWO

         // If ||Z||_2 is not one, then TEMP should be set to
         // RHO * ||Z||_2^2 / TWO

         for (J = 1; J <= N; J++) { // 10
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
         } // 10

         PSI = ZERO
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
         } // 20

         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / DELTA( II ) + Z( N )*Z( N ) / DELTA( N )

         if ( W.LE.ZERO ) {
            TEMP = Z( N-1 )*Z( N-1 ) / ( D( N )-D( N-1 )+RHO ) + Z( N )*Z( N ) / RHO
            if ( C.LE.TEMP ) {
               TAU = RHO
            } else {
               DEL = D( N ) - D( N-1 )
               A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DEL
               if ( A.LT.ZERO ) {
                  TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               } else {
                  TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               }
            }

            // It can be proved that
                // D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO

            DLTLB = MIDPT
            DLTUB = RHO
         } else {
            DEL = D( N ) - D( N-1 )
            A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DEL
            if ( A.LT.ZERO ) {
               TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            } else {
               TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            }

            // It can be proved that
                // D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2

            DLTLB = ZERO
            DLTUB = MIDPT
         }

         for (J = 1; J <= N; J++) { // 30
            DELTA( J ) = ( D( J )-D( I ) ) - TAU
         } // 30

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         for (J = 1; J <= II; J++) { // 40
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
         } // 40
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + ABS( TAU )*( DPSI+DPHI )

         W = RHOINV + PHI + PSI

         // Test for convergence

         if ( ABS( W ).LE.EPS*ERRETM ) {
            DLAM = D( I ) + TAU
            GO TO 250
         }

         if ( W.LE.ZERO ) {
            DLTLB = MAX( DLTLB, TAU )
         } else {
            DLTUB = MIN( DLTUB, TAU )
         }

         // Calculate the new step

         NITER = NITER + 1
         C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
         A = ( DELTA( N-1 )+DELTA( N ) )*W - DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
         B = DELTA( N-1 )*DELTA( N )*W
         IF( C.LT.ZERO ) C = ABS( C )
         if ( C.EQ.ZERO ) {
            // ETA = B/A
            // ETA = RHO - TAU
            // ETA = DLTUB - TAU

            // Update proposed by Li, Ren-Cang:
            ETA = -W / ( DPSI+DPHI )
         } else if ( A.GE.ZERO ) {
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         } else {
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         }

         // Note, eta should be positive if w is negative, and
         // eta should be negative otherwise. However,
         // if for some reason caused by roundoff, eta*w > 0,
         // we simply use one Newton step instead. This way
         // will guarantee eta*w < 0.

         IF( W*ETA.GT.ZERO ) ETA = -W / ( DPSI+DPHI )
         TEMP = TAU + ETA
         if ( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) {
            if ( W.LT.ZERO ) {
               ETA = ( DLTUB-TAU ) / TWO
            } else {
               ETA = ( DLTLB-TAU ) / TWO
            }
         }
         for (J = 1; J <= N; J++) { // 50
            DELTA( J ) = DELTA( J ) - ETA
         } // 50

         TAU = TAU + ETA

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         for (J = 1; J <= II; J++) { // 60
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
         } // 60
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + ABS( TAU )*( DPSI+DPHI )

         W = RHOINV + PHI + PSI

         // Main loop to update the values of the array   DELTA

         ITER = NITER + 1

         for (NITER = ITER; NITER <= MAXIT; NITER++) { // 90

            // Test for convergence

            if ( ABS( W ).LE.EPS*ERRETM ) {
               DLAM = D( I ) + TAU
               GO TO 250
            }

            if ( W.LE.ZERO ) {
               DLTLB = MAX( DLTLB, TAU )
            } else {
               DLTUB = MIN( DLTUB, TAU )
            }

            // Calculate the new step

            C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
            A = ( DELTA( N-1 )+DELTA( N ) )*W - DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
            B = DELTA( N-1 )*DELTA( N )*W
            if ( A.GE.ZERO ) {
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            } else {
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            }

            // Note, eta should be positive if w is negative, and
            // eta should be negative otherwise. However,
            // if for some reason caused by roundoff, eta*w > 0,
            // we simply use one Newton step instead. This way
            // will guarantee eta*w < 0.

            IF( W*ETA.GT.ZERO ) ETA = -W / ( DPSI+DPHI )
            TEMP = TAU + ETA
            if ( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) {
               if ( W.LT.ZERO ) {
                  ETA = ( DLTUB-TAU ) / TWO
               } else {
                  ETA = ( DLTLB-TAU ) / TWO
               }
            }
            for (J = 1; J <= N; J++) { // 70
               DELTA( J ) = DELTA( J ) - ETA
            } // 70

            TAU = TAU + ETA

            // Evaluate PSI and the derivative DPSI

            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            for (J = 1; J <= II; J++) { // 80
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
            } // 80
            ERRETM = ABS( ERRETM )

            // Evaluate PHI and the derivative DPHI

            TEMP = Z( N ) / DELTA( N )
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + ABS( TAU )*( DPSI+DPHI )

            W = RHOINV + PHI + PSI
         } // 90

         // Return with INFO = 1, NITER = MAXIT and not converged

         INFO = 1
         DLAM = D( I ) + TAU
         GO TO 250

         // End for the case I = N

      } else {

         // The case for I < N

         NITER = 1
         IP1 = I + 1

         // Calculate initial guess

         DEL = D( IP1 ) - D( I )
         MIDPT = DEL / TWO
         for (J = 1; J <= N; J++) { // 100
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
         } // 100

         PSI = ZERO
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
         } // 110

         PHI = ZERO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / DELTA( J )
         } // 120
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / DELTA( I ) + Z( IP1 )*Z( IP1 ) / DELTA( IP1 )

         if ( W.GT.ZERO ) {

            // d(i)< the ith eigenvalue < (d(i)+d(i+1))/2

            // We choose d(i) as origin.

            ORGATI = .TRUE.
            A = C*DEL + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DEL
            if ( A.GT.ZERO ) {
               TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            } else {
               TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            }
            DLTLB = ZERO
            DLTUB = MIDPT
         } else {

            // (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)

            // We choose d(i+1) as origin.

            ORGATI = .FALSE.
            A = C*DEL - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DEL
            if ( A.LT.ZERO ) {
               TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            } else {
               TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            }
            DLTLB = -MIDPT
            DLTUB = ZERO
         }

         if ( ORGATI ) {
            for (J = 1; J <= N; J++) { // 130
               DELTA( J ) = ( D( J )-D( I ) ) - TAU
            } // 130
         } else {
            for (J = 1; J <= N; J++) { // 140
               DELTA( J ) = ( D( J )-D( IP1 ) ) - TAU
            } // 140
         }
         if ( ORGATI ) {
            II = I
         } else {
            II = I + 1
         }
         IIM1 = II - 1
         IIP1 = II + 1

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         for (J = 1; J <= IIM1; J++) { // 150
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
         } // 150
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
         } // 160

         W = RHOINV + PHI + PSI

         // W is the value of the secular function with
         // its ii-th element removed.

         SWTCH3 = .FALSE.
         if ( ORGATI ) {
            IF( W.LT.ZERO ) SWTCH3 = .TRUE.
         } else {
            IF( W.GT.ZERO ) SWTCH3 = .TRUE.
         }
         IF( II.EQ.1 .OR. II.EQ.N ) SWTCH3 = .FALSE.

         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP ) + ABS( TAU )*DW

         // Test for convergence

         if ( ABS( W ).LE.EPS*ERRETM ) {
            if ( ORGATI ) {
               DLAM = D( I ) + TAU
            } else {
               DLAM = D( IP1 ) + TAU
            }
            GO TO 250
         }

         if ( W.LE.ZERO ) {
            DLTLB = MAX( DLTLB, TAU )
         } else {
            DLTUB = MIN( DLTUB, TAU )
         }

         // Calculate the new step

         NITER = NITER + 1
         if ( .NOT.SWTCH3 ) {
            if ( ORGATI ) {
               C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )* ( Z( I ) / DELTA( I ) )**2
            } else {
               C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )* ( Z( IP1 ) / DELTA( IP1 ) )**2
            }
            A = ( DELTA( I )+DELTA( IP1 ) )*W - DELTA( I )*DELTA( IP1 )*DW
            B = DELTA( I )*DELTA( IP1 )*W
            if ( C.EQ.ZERO ) {
               if ( A.EQ.ZERO ) {
                  if ( ORGATI ) {
                     A = Z( I )*Z( I ) + DELTA( IP1 )*DELTA( IP1 )* ( DPSI+DPHI )
                  } else {
                     A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )* ( DPSI+DPHI )
                  }
               }
               ETA = B / A
            } else if ( A.LE.ZERO ) {
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            } else {
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            }
         } else {

            // Interpolation using THREE most relevant poles

            TEMP = RHOINV + PSI + PHI
            if ( ORGATI ) {
               TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) - ( D( IIM1 )-D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )* ( ( DPSI-TEMP1 )+DPHI )
            } else {
               TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) - ( D( IIP1 )-D( IIM1 ) )*TEMP1                ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )* ( DPSI+( DPHI-TEMP1 ) )
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            }
            ZZ( 2 ) = Z( II )*Z( II )
            CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA, INFO )             IF( INFO.NE.0 ) GO TO 250
         }

         // Note, eta should be positive if w is negative, and
         // eta should be negative otherwise. However,
         // if for some reason caused by roundoff, eta*w > 0,
         // we simply use one Newton step instead. This way
         // will guarantee eta*w < 0.

         IF( W*ETA.GE.ZERO ) ETA = -W / DW
         TEMP = TAU + ETA
         if ( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) {
            if ( W.LT.ZERO ) {
               ETA = ( DLTUB-TAU ) / TWO
            } else {
               ETA = ( DLTLB-TAU ) / TWO
            }
         }

         PREW = W

         for (J = 1; J <= N; J++) { // 180
            DELTA( J ) = DELTA( J ) - ETA
         } // 180

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         for (J = 1; J <= IIM1; J++) { // 190
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
         } // 190
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         DPHI = ZERO
         PHI = ZERO
         DO 200 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
         } // 200

         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP ) + ABS( TAU+ETA )*DW

         SWTCH = .FALSE.
         if ( ORGATI ) {
            IF( -W.GT.ABS( PREW ) / TEN ) SWTCH = .TRUE.
         } else {
            IF( W.GT.ABS( PREW ) / TEN ) SWTCH = .TRUE.
         }

         TAU = TAU + ETA

         // Main loop to update the values of the array   DELTA

         ITER = NITER + 1

         for (NITER = ITER; NITER <= MAXIT; NITER++) { // 240

            // Test for convergence

            if ( ABS( W ).LE.EPS*ERRETM ) {
               if ( ORGATI ) {
                  DLAM = D( I ) + TAU
               } else {
                  DLAM = D( IP1 ) + TAU
               }
               GO TO 250
            }

            if ( W.LE.ZERO ) {
               DLTLB = MAX( DLTLB, TAU )
            } else {
               DLTUB = MIN( DLTUB, TAU )
            }

            // Calculate the new step

            if ( .NOT.SWTCH3 ) {
               if ( .NOT.SWTCH ) {
                  if ( ORGATI ) {
                     C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )*( Z( I ) / DELTA( I ) )**2
                  } else {
                     C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )* ( Z( IP1 ) / DELTA( IP1 ) )**2
                  }
               } else {
                  TEMP = Z( II ) / DELTA( II )
                  if ( ORGATI ) {
                     DPSI = DPSI + TEMP*TEMP
                  } else {
                     DPHI = DPHI + TEMP*TEMP
                  }
                  C = W - DELTA( I )*DPSI - DELTA( IP1 )*DPHI
               }
               A = ( DELTA( I )+DELTA( IP1 ) )*W - DELTA( I )*DELTA( IP1 )*DW
               B = DELTA( I )*DELTA( IP1 )*W
               if ( C.EQ.ZERO ) {
                  if ( A.EQ.ZERO ) {
                     if ( .NOT.SWTCH ) {
                        if ( ORGATI ) {
                           A = Z( I )*Z( I ) + DELTA( IP1 )* DELTA( IP1 )*( DPSI+DPHI )
                        } else {
                           A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )*( DPSI+DPHI )
                        }
                     } else {
                        A = DELTA( I )*DELTA( I )*DPSI + DELTA( IP1 )*DELTA( IP1 )*DPHI
                     }
                  }
                  ETA = B / A
               } else if ( A.LE.ZERO ) {
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               } else {
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               }
            } else {

               // Interpolation using THREE most relevant poles

               TEMP = RHOINV + PSI + PHI
               if ( SWTCH ) {
                  C = TEMP - DELTA( IIM1 )*DPSI - DELTA( IIP1 )*DPHI
                  ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*DPSI
                  ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*DPHI
               } else {
                  if ( ORGATI ) {
                     TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) - ( D( IIM1 )-D( IIP1 ) )*TEMP1
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )* ( ( DPSI-TEMP1 )+DPHI )
                  } else {
                     TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) - ( D( IIP1 )-D( IIM1 ) )*TEMP1                      ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )* ( DPSI+( DPHI-TEMP1 ) )
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  }
               }
               CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA, INFO )                IF( INFO.NE.0 ) GO TO 250
            }

            // Note, eta should be positive if w is negative, and
            // eta should be negative otherwise. However,
            // if for some reason caused by roundoff, eta*w > 0,
            // we simply use one Newton step instead. This way
            // will guarantee eta*w < 0.

            IF( W*ETA.GE.ZERO ) ETA = -W / DW
            TEMP = TAU + ETA
            if ( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) {
               if ( W.LT.ZERO ) {
                  ETA = ( DLTUB-TAU ) / TWO
               } else {
                  ETA = ( DLTLB-TAU ) / TWO
               }
            }

            for (J = 1; J <= N; J++) { // 210
               DELTA( J ) = DELTA( J ) - ETA
            } // 210

            TAU = TAU + ETA
            PREW = W

            // Evaluate PSI and the derivative DPSI

            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            for (J = 1; J <= IIM1; J++) { // 220
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
            } // 220
            ERRETM = ABS( ERRETM )

            // Evaluate PHI and the derivative DPHI

            DPHI = ZERO
            PHI = ZERO
            DO 230 J = N, IIP1, -1
               TEMP = Z( J ) / DELTA( J )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
            } // 230

            TEMP = Z( II ) / DELTA( II )
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP ) + ABS( TAU )*DW             IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN ) SWTCH = .NOT.SWTCH

         } // 240

         // Return with INFO = 1, NITER = MAXIT and not converged

         INFO = 1
         if ( ORGATI ) {
            DLAM = D( I ) + TAU
         } else {
            DLAM = D( IP1 ) + TAU
         }

      }

      } // 250

      RETURN

      // End of DLAED4

      }
