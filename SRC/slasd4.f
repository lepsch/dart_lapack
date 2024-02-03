      SUBROUTINE SLASD4( N, I, D, Z, DELTA, RHO, SIGMA, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                I, INFO, N;
      REAL   RHO, SIGMA
      // ..
      // .. Array Arguments ..
      REAL   D( * ), DELTA( * ), WORK( * ), Z( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                MAXIT;
      const              MAXIT = 400 ;
      REAL               ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0, THREE = 3.0E+0, FOUR = 4.0E+0, EIGHT = 8.0E+0, TEN = 10.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ORGATI, SWTCH, SWTCH3, GEOMAVG;
      int                II, IIM1, IIP1, IP1, ITER, J, NITER;
      REAL               A, B, C, DELSQ, DELSQ2, SQ2, DPHI, DPSI, DTIIM, DTIIP, DTIPSQ, DTISQ, DTNSQ, DTNSQ1, DW, EPS, ERRETM, ETA, PHI, PREW, PSI, RHOINV, SGLB, SGUB, TAU, TAU2, TEMP, TEMP1, TEMP2, W
      // ..
      // .. Local Arrays ..
      REAL               DD( 3 ), ZZ( 3 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAED6, SLASD5
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Since this routine is called in an inner loop, we do no argument
      // checking.

      // Quick return for N=1 and 2.

      INFO = 0
      if ( N == 1 ) {

         // Presumably, I=1 upon entry

         SIGMA = SQRT( D( 1 )*D( 1 )+RHO*Z( 1 )*Z( 1 ) )
         DELTA( 1 ) = ONE
         WORK( 1 ) = ONE
         RETURN
      }
      if ( N == 2 ) {
         slasd5(I, D, Z, DELTA, RHO, SIGMA, WORK );
         RETURN
      }

      // Compute machine epsilon

      EPS = SLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
      TAU2= ZERO

      // The case I = N

      if ( I == N ) {

         // Initialize some basic variables

         II = N - 1
         NITER = 1

         // Calculate initial guess

         TEMP = RHO / TWO

         // If ||Z||_2 is not one, then TEMP should be set to
         // RHO * ||Z||_2^2 / TWO

         TEMP1 = TEMP / ( D( N )+SQRT( D( N )*D( N )+TEMP ) )
         for (J = 1; J <= N; J++) { // 10
            WORK( J ) = D( J ) + D( N ) + TEMP1
            DELTA( J ) = ( D( J )-D( N ) ) - TEMP1
         } // 10

         PSI = ZERO
         for (J = 1; J <= N - 2; J++) { // 20
            PSI = PSI + Z( J )*Z( J ) / ( DELTA( J )*WORK( J ) )
         } // 20

         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / ( DELTA( II )*WORK( II ) ) + Z( N )*Z( N ) / ( DELTA( N )*WORK( N ) )

         if ( W.LE.ZERO ) {
            TEMP1 = SQRT( D( N )*D( N )+RHO )
            TEMP = Z( N-1 )*Z( N-1 ) / ( ( D( N-1 )+TEMP1 )* ( D( N )-D( N-1 )+RHO / ( D( N )+TEMP1 ) ) ) + Z( N )*Z( N ) / RHO

            // The following TAU2 is to approximate
            // SIGMA_n^2 - D( N )*D( N )

            if ( C.LE.TEMP ) {
               TAU = RHO
            } else {
               DELSQ = ( D( N )-D( N-1 ) )*( D( N )+D( N-1 ) )
               A = -C*DELSQ + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DELSQ
               if ( A < ZERO ) {
                  TAU2 = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               } else {
                  TAU2 = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               }
               TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )
            }

            // It can be proved that
                // D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU2 <= D(N)^2+RHO

         } else {
            DELSQ = ( D( N )-D( N-1 ) )*( D( N )+D( N-1 ) )
            A = -C*DELSQ + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DELSQ

            // The following TAU2 is to approximate
            // SIGMA_n^2 - D( N )*D( N )

            if ( A < ZERO ) {
               TAU2 = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            } else {
               TAU2 = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            }
            TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )


            // It can be proved that
            // D(N)^2 < D(N)^2+TAU2 < SIGMA(N)^2 < D(N)^2+RHO/2

         }

         // The following TAU is to approximate SIGMA_n - D( N )

          // TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )

         SIGMA = D( N ) + TAU
         for (J = 1; J <= N; J++) { // 30
            DELTA( J ) = ( D( J )-D( N ) ) - TAU
            WORK( J ) = D( J ) + D( N ) + TAU
         } // 30

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         for (J = 1; J <= II; J++) { // 40
            TEMP = Z( J ) / ( DELTA( J )*WORK( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
         } // 40
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         TEMP = Z( N ) / ( DELTA( N )*WORK( N ) )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV
*    $          + ABS( TAU2 )*( DPSI+DPHI )

         W = RHOINV + PHI + PSI

         // Test for convergence

         if ( ABS( W ).LE.EPS*ERRETM ) {
            GO TO 240
         }

         // Calculate the new step

         NITER = NITER + 1
         DTNSQ1 = WORK( N-1 )*DELTA( N-1 )
         DTNSQ = WORK( N )*DELTA( N )
         C = W - DTNSQ1*DPSI - DTNSQ*DPHI
         A = ( DTNSQ+DTNSQ1 )*W - DTNSQ*DTNSQ1*( DPSI+DPHI )
         B = DTNSQ*DTNSQ1*W
         if (C < ZERO) C = ABS( C );
         if ( C == ZERO ) {
            ETA = RHO - SIGMA*SIGMA
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

         if (W*ETA > ZERO) ETA = -W / ( DPSI+DPHI );
         TEMP = ETA - DTNSQ
         if (TEMP > RHO) ETA = RHO + DTNSQ;

         ETA = ETA / ( SIGMA+SQRT( ETA+SIGMA*SIGMA ) )
         TAU = TAU + ETA
         SIGMA = SIGMA + ETA

         for (J = 1; J <= N; J++) { // 50
            DELTA( J ) = DELTA( J ) - ETA
            WORK( J ) = WORK( J ) + ETA
         } // 50

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         for (J = 1; J <= II; J++) { // 60
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
         } // 60
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         TAU2 = WORK( N )*DELTA( N )
         TEMP = Z( N ) / TAU2
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV
*    $          + ABS( TAU2 )*( DPSI+DPHI )

         W = RHOINV + PHI + PSI

         // Main loop to update the values of the array   DELTA

         ITER = NITER + 1

         for (NITER = ITER; NITER <= MAXIT; NITER++) { // 90

            // Test for convergence

            if ( ABS( W ).LE.EPS*ERRETM ) {
               GO TO 240
            }

            // Calculate the new step

            DTNSQ1 = WORK( N-1 )*DELTA( N-1 )
            DTNSQ = WORK( N )*DELTA( N )
            C = W - DTNSQ1*DPSI - DTNSQ*DPHI
            A = ( DTNSQ+DTNSQ1 )*W - DTNSQ1*DTNSQ*( DPSI+DPHI )
            B = DTNSQ1*DTNSQ*W
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

            if (W*ETA > ZERO) ETA = -W / ( DPSI+DPHI );
            TEMP = ETA - DTNSQ
            if (TEMP.LE.ZERO) ETA = ETA / TWO;

            ETA = ETA / ( SIGMA+SQRT( ETA+SIGMA*SIGMA ) )
            TAU = TAU + ETA
            SIGMA = SIGMA + ETA

            for (J = 1; J <= N; J++) { // 70
               DELTA( J ) = DELTA( J ) - ETA
               WORK( J ) = WORK( J ) + ETA
            } // 70

            // Evaluate PSI and the derivative DPSI

            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            for (J = 1; J <= II; J++) { // 80
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
            } // 80
            ERRETM = ABS( ERRETM )

            // Evaluate PHI and the derivative DPHI

            TAU2 = WORK( N )*DELTA( N )
            TEMP = Z( N ) / TAU2
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV
*    $             + ABS( TAU2 )*( DPSI+DPHI )

            W = RHOINV + PHI + PSI
         } // 90

         // Return with INFO = 1, NITER = MAXIT and not converged

         INFO = 1
         GO TO 240

         // End for the case I = N

      } else {

         // The case for I < N

         NITER = 1
         IP1 = I + 1

         // Calculate initial guess

         DELSQ = ( D( IP1 )-D( I ) )*( D( IP1 )+D( I ) )
         DELSQ2 = DELSQ / TWO
         SQ2=SQRT( ( D( I )*D( I )+D( IP1 )*D( IP1 ) ) / TWO )
         TEMP = DELSQ2 / ( D( I )+SQ2 )
         for (J = 1; J <= N; J++) { // 100
            WORK( J ) = D( J ) + D( I ) + TEMP
            DELTA( J ) = ( D( J )-D( I ) ) - TEMP
         } // 100

         PSI = ZERO
         for (J = 1; J <= I - 1; J++) { // 110
            PSI = PSI + Z( J )*Z( J ) / ( WORK( J )*DELTA( J ) )
         } // 110

         PHI = ZERO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / ( WORK( J )*DELTA( J ) )
         } // 120
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / ( WORK( I )*DELTA( I ) ) + Z( IP1 )*Z( IP1 ) / ( WORK( IP1 )*DELTA( IP1 ) )

         GEOMAVG = false;
         if ( W > ZERO ) {

            // d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2

            // We choose d(i) as origin.

            ORGATI = true;
            II = I
            SGLB = ZERO
            SGUB = DELSQ2  / ( D( I )+SQ2 )
            A = C*DELSQ + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DELSQ
            if ( A > ZERO ) {
               TAU2 = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            } else {
               TAU2 = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            }

            // TAU2 now is an estimation of SIGMA^2 - D( I )^2. The
            // following, however, is the corresponding estimation of
            // SIGMA - D( I ).

            TAU = TAU2 / ( D( I )+SQRT( D( I )*D( I )+TAU2 ) )
            TEMP = SQRT(EPS)
            if ( (D(I).LE.TEMP*D(IP1)) && (ABS(Z(I)).LE.TEMP) && (D(I) > ZERO) ) {
               TAU = MIN( TEN*D(I), SGUB )
               GEOMAVG = true;
            }
         } else {

            // (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2

            // We choose d(i+1) as origin.

            ORGATI = false;
            II = IP1
            SGLB = -DELSQ2  / ( D( II )+SQ2 )
            SGUB = ZERO
            A = C*DELSQ - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DELSQ
            if ( A < ZERO ) {
               TAU2 = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            } else {
               TAU2 = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            }

            // TAU2 now is an estimation of SIGMA^2 - D( IP1 )^2. The
            // following, however, is the corresponding estimation of
            // SIGMA - D( IP1 ).

            TAU = TAU2 / ( D( IP1 )+SQRT( ABS( D( IP1 )*D( IP1 )+ TAU2 ) ) )
         }

         SIGMA = D( II ) + TAU
         for (J = 1; J <= N; J++) { // 130
            WORK( J ) = D( J ) + D( II ) + TAU
            DELTA( J ) = ( D( J )-D( II ) ) - TAU
         } // 130
         IIM1 = II - 1
         IIP1 = II + 1

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         for (J = 1; J <= IIM1; J++) { // 150
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
         } // 150
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
         } // 160

         W = RHOINV + PHI + PSI

         // W is the value of the secular function with
         // its ii-th element removed.

         SWTCH3 = false;
         if ( ORGATI ) {
            if (W < ZERO) SWTCH3 = true ;
         } else {
            if (W > ZERO) SWTCH3 = true ;
         }
         if (II == 1 || II == N) SWTCH3 = false ;

         TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP )
*    $          + ABS( TAU2 )*DW

         // Test for convergence

         if ( ABS( W ).LE.EPS*ERRETM ) {
            GO TO 240
         }

         if ( W.LE.ZERO ) {
            SGLB = MAX( SGLB, TAU )
         } else {
            SGUB = MIN( SGUB, TAU )
         }

         // Calculate the new step

         NITER = NITER + 1
         if ( .NOT.SWTCH3 ) {
            DTIPSQ = WORK( IP1 )*DELTA( IP1 )
            DTISQ = WORK( I )*DELTA( I )
            if ( ORGATI ) {
               C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
            } else {
               C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
            }
            A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
            B = DTIPSQ*DTISQ*W
            if ( C == ZERO ) {
               if ( A == ZERO ) {
                  if ( ORGATI ) {
                     A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*( DPSI+DPHI )
                  } else {
                     A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI )
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

            DTIIM = WORK( IIM1 )*DELTA( IIM1 )
            DTIIP = WORK( IIP1 )*DELTA( IIP1 )
            TEMP = RHOINV + PSI + PHI
            if ( ORGATI ) {
               TEMP1 = Z( IIM1 ) / DTIIM
               TEMP1 = TEMP1*TEMP1
               C = ( TEMP - DTIIP*( DPSI+DPHI ) ) - ( D( IIM1 )-D( IIP1 ) )*( D( IIM1 )+D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               if ( DPSI < TEMP1 ) {
                  ZZ( 3 ) = DTIIP*DTIIP*DPHI
               } else {
                  ZZ( 3 ) = DTIIP*DTIIP*( ( DPSI-TEMP1 )+DPHI )
               }
            } else {
               TEMP1 = Z( IIP1 ) / DTIIP
               TEMP1 = TEMP1*TEMP1
               C = ( TEMP - DTIIM*( DPSI+DPHI ) ) - ( D( IIP1 )-D( IIM1 ) )*( D( IIM1 )+D( IIP1 ) )*TEMP1
               if ( DPHI < TEMP1 ) {
                  ZZ( 1 ) = DTIIM*DTIIM*DPSI
               } else {
                  ZZ( 1 ) = DTIIM*DTIIM*( DPSI+( DPHI-TEMP1 ) )
               }
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            }
            ZZ( 2 ) = Z( II )*Z( II )
            DD( 1 ) = DTIIM
            DD( 2 ) = DELTA( II )*WORK( II )
            DD( 3 ) = DTIIP
            slaed6(NITER, ORGATI, C, DD, ZZ, W, ETA, INFO );

            if ( INFO != 0 ) {

               // If INFO is not 0, i.e., SLAED6 failed, switch back
               // to 2 pole interpolation.

               SWTCH3 = false;
               INFO = 0
               DTIPSQ = WORK( IP1 )*DELTA( IP1 )
               DTISQ = WORK( I )*DELTA( I )
               if ( ORGATI ) {
                  C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
               } else {
                  C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
               }
               A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
               B = DTIPSQ*DTISQ*W
               if ( C == ZERO ) {
                  if ( A == ZERO ) {
                     if ( ORGATI ) {
                        A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*( DPSI+DPHI )
                     } else {
                        A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI)
                     }
                  }
                  ETA = B / A
               } else if ( A.LE.ZERO ) {
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               } else {
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               }
            }
         }

         // Note, eta should be positive if w is negative, and
         // eta should be negative otherwise. However,
         // if for some reason caused by roundoff, eta*w > 0,
         // we simply use one Newton step instead. This way
         // will guarantee eta*w < 0.

         if (W*ETA.GE.ZERO) ETA = -W / DW;

         ETA = ETA / ( SIGMA+SQRT( SIGMA*SIGMA+ETA ) )
         TEMP = TAU + ETA
         if ( TEMP > SGUB || TEMP < SGLB ) {
            if ( W < ZERO ) {
               ETA = ( SGUB-TAU ) / TWO
            } else {
               ETA = ( SGLB-TAU ) / TWO
            }
            if ( GEOMAVG ) {
               if ( W < ZERO ) {
                  if ( TAU > ZERO ) {
                     ETA = SQRT(SGUB*TAU)-TAU
                  }
               } else {
                  if ( SGLB > ZERO ) {
                     ETA = SQRT(SGLB*TAU)-TAU
                  }
               }
            }
         }

         PREW = W

         TAU = TAU + ETA
         SIGMA = SIGMA + ETA

         for (J = 1; J <= N; J++) { // 170
            WORK( J ) = WORK( J ) + ETA
            DELTA( J ) = DELTA( J ) - ETA
         } // 170

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         for (J = 1; J <= IIM1; J++) { // 180
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
         } // 180
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         DPHI = ZERO
         PHI = ZERO
         DO 190 J = N, IIP1, -1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
         } // 190

         TAU2 = WORK( II )*DELTA( II )
         TEMP = Z( II ) / TAU2
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP )
*    $          + ABS( TAU2 )*DW

         SWTCH = false;
         if ( ORGATI ) {
            IF( -W > ABS( PREW ) / TEN ) SWTCH = true;
         } else {
            IF( W > ABS( PREW ) / TEN ) SWTCH = true;
         }

         // Main loop to update the values of the array   DELTA and WORK

         ITER = NITER + 1

         for (NITER = ITER; NITER <= MAXIT; NITER++) { // 230

            // Test for convergence

            if ( ABS( W ).LE.EPS*ERRETM ) {
      // $ || (SGUB-SGLB).LE.EIGHT*ABS(SGUB+SGLB) ) THEN
               GO TO 240
            }

            if ( W.LE.ZERO ) {
               SGLB = MAX( SGLB, TAU )
            } else {
               SGUB = MIN( SGUB, TAU )
            }

            // Calculate the new step

            if ( .NOT.SWTCH3 ) {
               DTIPSQ = WORK( IP1 )*DELTA( IP1 )
               DTISQ = WORK( I )*DELTA( I )
               if ( .NOT.SWTCH ) {
                  if ( ORGATI ) {
                     C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
                  } else {
                     C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
                  }
               } else {
                  TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
                  if ( ORGATI ) {
                     DPSI = DPSI + TEMP*TEMP
                  } else {
                     DPHI = DPHI + TEMP*TEMP
                  }
                  C = W - DTISQ*DPSI - DTIPSQ*DPHI
               }
               A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
               B = DTIPSQ*DTISQ*W
               if ( C == ZERO ) {
                  if ( A == ZERO ) {
                     if ( .NOT.SWTCH ) {
                        if ( ORGATI ) {
                           A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ* ( DPSI+DPHI )
                        } else {
                           A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI )
                        }
                     } else {
                        A = DTISQ*DTISQ*DPSI + DTIPSQ*DTIPSQ*DPHI
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

               DTIIM = WORK( IIM1 )*DELTA( IIM1 )
               DTIIP = WORK( IIP1 )*DELTA( IIP1 )
               TEMP = RHOINV + PSI + PHI
               if ( SWTCH ) {
                  C = TEMP - DTIIM*DPSI - DTIIP*DPHI
                  ZZ( 1 ) = DTIIM*DTIIM*DPSI
                  ZZ( 3 ) = DTIIP*DTIIP*DPHI
               } else {
                  if ( ORGATI ) {
                     TEMP1 = Z( IIM1 ) / DTIIM
                     TEMP1 = TEMP1*TEMP1
                     TEMP2 = ( D( IIM1 )-D( IIP1 ) )* ( D( IIM1 )+D( IIP1 ) )*TEMP1
                     C = TEMP - DTIIP*( DPSI+DPHI ) - TEMP2
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     if ( DPSI < TEMP1 ) {
                        ZZ( 3 ) = DTIIP*DTIIP*DPHI
                     } else {
                        ZZ( 3 ) = DTIIP*DTIIP*( ( DPSI-TEMP1 )+DPHI )
                     }
                  } else {
                     TEMP1 = Z( IIP1 ) / DTIIP
                     TEMP1 = TEMP1*TEMP1
                     TEMP2 = ( D( IIP1 )-D( IIM1 ) )* ( D( IIM1 )+D( IIP1 ) )*TEMP1
                     C = TEMP - DTIIM*( DPSI+DPHI ) - TEMP2
                     if ( DPHI < TEMP1 ) {
                        ZZ( 1 ) = DTIIM*DTIIM*DPSI
                     } else {
                        ZZ( 1 ) = DTIIM*DTIIM*( DPSI+( DPHI-TEMP1 ) )
                     }
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  }
               }
               DD( 1 ) = DTIIM
               DD( 2 ) = DELTA( II )*WORK( II )
               DD( 3 ) = DTIIP
               slaed6(NITER, ORGATI, C, DD, ZZ, W, ETA, INFO );

               if ( INFO != 0 ) {

                  // If INFO is not 0, i.e., SLAED6 failed, switch
                  // back to two pole interpolation

                  SWTCH3 = false;
                  INFO = 0
                  DTIPSQ = WORK( IP1 )*DELTA( IP1 )
                  DTISQ = WORK( I )*DELTA( I )
                  if ( .NOT.SWTCH ) {
                     if ( ORGATI ) {
                        C = W - DTIPSQ*DW + DELSQ*( Z( I )/DTISQ )**2
                     } else {
                        C = W - DTISQ*DW - DELSQ*( Z( IP1 )/DTIPSQ )**2
                     }
                  } else {
                     TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
                     if ( ORGATI ) {
                        DPSI = DPSI + TEMP*TEMP
                     } else {
                        DPHI = DPHI + TEMP*TEMP
                     }
                     C = W - DTISQ*DPSI - DTIPSQ*DPHI
                  }
                  A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
                  B = DTIPSQ*DTISQ*W
                  if ( C == ZERO ) {
                     if ( A == ZERO ) {
                        if ( .NOT.SWTCH ) {
                           if ( ORGATI ) {
                              A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ* ( DPSI+DPHI )
                           } else {
                              A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI )
                           }
                        } else {
                           A = DTISQ*DTISQ*DPSI + DTIPSQ*DTIPSQ*DPHI
                        }
                     }
                     ETA = B / A
                  } else if ( A.LE.ZERO ) {
                     ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
                  } else {
                     ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
                  }
               }
            }

            // Note, eta should be positive if w is negative, and
            // eta should be negative otherwise. However,
            // if for some reason caused by roundoff, eta*w > 0,
            // we simply use one Newton step instead. This way
            // will guarantee eta*w < 0.

            if (W*ETA.GE.ZERO) ETA = -W / DW;

            ETA = ETA / ( SIGMA+SQRT( SIGMA*SIGMA+ETA ) )
            TEMP=TAU+ETA
            if ( TEMP > SGUB || TEMP < SGLB ) {
               if ( W < ZERO ) {
                  ETA = ( SGUB-TAU ) / TWO
               } else {
                  ETA = ( SGLB-TAU ) / TWO
               }
               if ( GEOMAVG ) {
                  if ( W < ZERO ) {
                     if ( TAU > ZERO ) {
                        ETA = SQRT(SGUB*TAU)-TAU
                     }
                  } else {
                     if ( SGLB > ZERO ) {
                        ETA = SQRT(SGLB*TAU)-TAU
                     }
                  }
               }
            }

            PREW = W

            TAU = TAU + ETA
            SIGMA = SIGMA + ETA

            for (J = 1; J <= N; J++) { // 200
               WORK( J ) = WORK( J ) + ETA
               DELTA( J ) = DELTA( J ) - ETA
            } // 200

            // Evaluate PSI and the derivative DPSI

            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            for (J = 1; J <= IIM1; J++) { // 210
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
            } // 210
            ERRETM = ABS( ERRETM )

            // Evaluate PHI and the derivative DPHI

            DPHI = ZERO
            PHI = ZERO
            DO 220 J = N, IIP1, -1
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
            } // 220

            TAU2 = WORK( II )*DELTA( II )
            TEMP = Z( II ) / TAU2
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP )
*    $             + ABS( TAU2 )*DW

            IF( W*PREW > ZERO && ABS( W ) > ABS( PREW ) / TEN ) SWTCH = .NOT.SWTCH

         } // 230

         // Return with INFO = 1, NITER = MAXIT and not converged

         INFO = 1

      }

      } // 240
      RETURN

      // End of SLASD4

      }
