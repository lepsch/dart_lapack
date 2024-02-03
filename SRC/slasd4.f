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
      IF( N.EQ.1 ) THEN

         // Presumably, I=1 upon entry

         SIGMA = SQRT( D( 1 )*D( 1 )+RHO*Z( 1 )*Z( 1 ) )
         DELTA( 1 ) = ONE
         WORK( 1 ) = ONE
         RETURN
      END IF
      IF( N.EQ.2 ) THEN
         CALL SLASD5( I, D, Z, DELTA, RHO, SIGMA, WORK )
         RETURN
      END IF

      // Compute machine epsilon

      EPS = SLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
      TAU2= ZERO

      // The case I = N

      IF( I.EQ.N ) THEN

         // Initialize some basic variables

         II = N - 1
         NITER = 1

         // Calculate initial guess

         TEMP = RHO / TWO

         // If ||Z||_2 is not one, then TEMP should be set to
         // RHO * ||Z||_2^2 / TWO

         TEMP1 = TEMP / ( D( N )+SQRT( D( N )*D( N )+TEMP ) )
         DO 10 J = 1, N
            WORK( J ) = D( J ) + D( N ) + TEMP1
            DELTA( J ) = ( D( J )-D( N ) ) - TEMP1
   10    CONTINUE

         PSI = ZERO
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / ( DELTA( J )*WORK( J ) )
   20    CONTINUE

         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / ( DELTA( II )*WORK( II ) ) + Z( N )*Z( N ) / ( DELTA( N )*WORK( N ) )

         IF( W.LE.ZERO ) THEN
            TEMP1 = SQRT( D( N )*D( N )+RHO )
            TEMP = Z( N-1 )*Z( N-1 ) / ( ( D( N-1 )+TEMP1 )* ( D( N )-D( N-1 )+RHO / ( D( N )+TEMP1 ) ) ) + Z( N )*Z( N ) / RHO

            // The following TAU2 is to approximate
            // SIGMA_n^2 - D( N )*D( N )

            IF( C.LE.TEMP ) THEN
               TAU = RHO
            } else {
               DELSQ = ( D( N )-D( N-1 ) )*( D( N )+D( N-1 ) )
               A = -C*DELSQ + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DELSQ
               IF( A.LT.ZERO ) THEN
                  TAU2 = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               } else {
                  TAU2 = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               END IF
               TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )
            END IF

            // It can be proved that
                // D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU2 <= D(N)^2+RHO

         } else {
            DELSQ = ( D( N )-D( N-1 ) )*( D( N )+D( N-1 ) )
            A = -C*DELSQ + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DELSQ

            // The following TAU2 is to approximate
            // SIGMA_n^2 - D( N )*D( N )

            IF( A.LT.ZERO ) THEN
               TAU2 = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            } else {
               TAU2 = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
            TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )


            // It can be proved that
            // D(N)^2 < D(N)^2+TAU2 < SIGMA(N)^2 < D(N)^2+RHO/2

         END IF

         // The following TAU is to approximate SIGMA_n - D( N )

          // TAU = TAU2 / ( D( N )+SQRT( D( N )*D( N )+TAU2 ) )

         SIGMA = D( N ) + TAU
         DO 30 J = 1, N
            DELTA( J ) = ( D( J )-D( N ) ) - TAU
            WORK( J ) = D( J ) + D( N ) + TAU
   30    CONTINUE

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 40 J = 1, II
            TEMP = Z( J ) / ( DELTA( J )*WORK( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   40    CONTINUE
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         TEMP = Z( N ) / ( DELTA( N )*WORK( N ) )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV
*    $          + ABS( TAU2 )*( DPSI+DPHI )

         W = RHOINV + PHI + PSI

         // Test for convergence

         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            GO TO 240
         END IF

         // Calculate the new step

         NITER = NITER + 1
         DTNSQ1 = WORK( N-1 )*DELTA( N-1 )
         DTNSQ = WORK( N )*DELTA( N )
         C = W - DTNSQ1*DPSI - DTNSQ*DPHI
         A = ( DTNSQ+DTNSQ1 )*W - DTNSQ*DTNSQ1*( DPSI+DPHI )
         B = DTNSQ*DTNSQ1*W
         IF( C.LT.ZERO ) C = ABS( C )
         IF( C.EQ.ZERO ) THEN
            ETA = RHO - SIGMA*SIGMA
         ELSE IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         } else {
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF

         // Note, eta should be positive if w is negative, and
         // eta should be negative otherwise. However,
         // if for some reason caused by roundoff, eta*w > 0,
         // we simply use one Newton step instead. This way
         // will guarantee eta*w < 0.

         IF( W*ETA.GT.ZERO ) ETA = -W / ( DPSI+DPHI )
         TEMP = ETA - DTNSQ
         IF( TEMP.GT.RHO ) ETA = RHO + DTNSQ

         ETA = ETA / ( SIGMA+SQRT( ETA+SIGMA*SIGMA ) )
         TAU = TAU + ETA
         SIGMA = SIGMA + ETA

         DO 50 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
            WORK( J ) = WORK( J ) + ETA
   50    CONTINUE

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 60 J = 1, II
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   60    CONTINUE
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

         DO 90 NITER = ITER, MAXIT

            // Test for convergence

            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               GO TO 240
            END IF

            // Calculate the new step

            DTNSQ1 = WORK( N-1 )*DELTA( N-1 )
            DTNSQ = WORK( N )*DELTA( N )
            C = W - DTNSQ1*DPSI - DTNSQ*DPHI
            A = ( DTNSQ+DTNSQ1 )*W - DTNSQ1*DTNSQ*( DPSI+DPHI )
            B = DTNSQ1*DTNSQ*W
            IF( A.GE.ZERO ) THEN
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            } else {
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF

            // Note, eta should be positive if w is negative, and
            // eta should be negative otherwise. However,
            // if for some reason caused by roundoff, eta*w > 0,
            // we simply use one Newton step instead. This way
            // will guarantee eta*w < 0.

            IF( W*ETA.GT.ZERO ) ETA = -W / ( DPSI+DPHI )
            TEMP = ETA - DTNSQ
            IF( TEMP.LE.ZERO ) ETA = ETA / TWO

            ETA = ETA / ( SIGMA+SQRT( ETA+SIGMA*SIGMA ) )
            TAU = TAU + ETA
            SIGMA = SIGMA + ETA

            DO 70 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
               WORK( J ) = WORK( J ) + ETA
   70       CONTINUE

            // Evaluate PSI and the derivative DPSI

            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 80 J = 1, II
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
   80       CONTINUE
            ERRETM = ABS( ERRETM )

            // Evaluate PHI and the derivative DPHI

            TAU2 = WORK( N )*DELTA( N )
            TEMP = Z( N ) / TAU2
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV
*    $             + ABS( TAU2 )*( DPSI+DPHI )

            W = RHOINV + PHI + PSI
   90    CONTINUE

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
         DO 100 J = 1, N
            WORK( J ) = D( J ) + D( I ) + TEMP
            DELTA( J ) = ( D( J )-D( I ) ) - TEMP
  100    CONTINUE

         PSI = ZERO
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / ( WORK( J )*DELTA( J ) )
  110    CONTINUE

         PHI = ZERO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / ( WORK( J )*DELTA( J ) )
  120    CONTINUE
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / ( WORK( I )*DELTA( I ) ) + Z( IP1 )*Z( IP1 ) / ( WORK( IP1 )*DELTA( IP1 ) )

         GEOMAVG = .FALSE.
         IF( W.GT.ZERO ) THEN

            // d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2

            // We choose d(i) as origin.

            ORGATI = .TRUE.
            II = I
            SGLB = ZERO
            SGUB = DELSQ2  / ( D( I )+SQ2 )
            A = C*DELSQ + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DELSQ
            IF( A.GT.ZERO ) THEN
               TAU2 = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            } else {
               TAU2 = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            END IF

            // TAU2 now is an estimation of SIGMA^2 - D( I )^2. The
            // following, however, is the corresponding estimation of
            // SIGMA - D( I ).

            TAU = TAU2 / ( D( I )+SQRT( D( I )*D( I )+TAU2 ) )
            TEMP = SQRT(EPS)
            IF( (D(I).LE.TEMP*D(IP1)).AND.(ABS(Z(I)).LE.TEMP) .AND.(D(I).GT.ZERO) ) THEN
               TAU = MIN( TEN*D(I), SGUB )
               GEOMAVG = .TRUE.
            END IF
         } else {

            // (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2

            // We choose d(i+1) as origin.

            ORGATI = .FALSE.
            II = IP1
            SGLB = -DELSQ2  / ( D( II )+SQ2 )
            SGUB = ZERO
            A = C*DELSQ - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DELSQ
            IF( A.LT.ZERO ) THEN
               TAU2 = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            } else {
               TAU2 = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            END IF

            // TAU2 now is an estimation of SIGMA^2 - D( IP1 )^2. The
            // following, however, is the corresponding estimation of
            // SIGMA - D( IP1 ).

            TAU = TAU2 / ( D( IP1 )+SQRT( ABS( D( IP1 )*D( IP1 )+ TAU2 ) ) )
         END IF

         SIGMA = D( II ) + TAU
         DO 130 J = 1, N
            WORK( J ) = D( J ) + D( II ) + TAU
            DELTA( J ) = ( D( J )-D( II ) ) - TAU
  130    CONTINUE
         IIM1 = II - 1
         IIP1 = II + 1

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 150 J = 1, IIM1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  150    CONTINUE
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  160    CONTINUE

         W = RHOINV + PHI + PSI

         // W is the value of the secular function with
         // its ii-th element removed.

         SWTCH3 = .FALSE.
         IF( ORGATI ) THEN
            IF( W.LT.ZERO ) SWTCH3 = .TRUE.
         } else {
            IF( W.GT.ZERO ) SWTCH3 = .TRUE.
         END IF
         IF( II.EQ.1 .OR. II.EQ.N ) SWTCH3 = .FALSE.

         TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP )
*    $          + ABS( TAU2 )*DW

         // Test for convergence

         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            GO TO 240
         END IF

         IF( W.LE.ZERO ) THEN
            SGLB = MAX( SGLB, TAU )
         } else {
            SGUB = MIN( SGUB, TAU )
         END IF

         // Calculate the new step

         NITER = NITER + 1
         IF( .NOT.SWTCH3 ) THEN
            DTIPSQ = WORK( IP1 )*DELTA( IP1 )
            DTISQ = WORK( I )*DELTA( I )
            IF( ORGATI ) THEN
               C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
            } else {
               C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
            END IF
            A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
            B = DTIPSQ*DTISQ*W
            IF( C.EQ.ZERO ) THEN
               IF( A.EQ.ZERO ) THEN
                  IF( ORGATI ) THEN
                     A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*( DPSI+DPHI )
                  } else {
                     A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI )
                  END IF
               END IF
               ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            } else {
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
         } else {

            // Interpolation using THREE most relevant poles

            DTIIM = WORK( IIM1 )*DELTA( IIM1 )
            DTIIP = WORK( IIP1 )*DELTA( IIP1 )
            TEMP = RHOINV + PSI + PHI
            IF( ORGATI ) THEN
               TEMP1 = Z( IIM1 ) / DTIIM
               TEMP1 = TEMP1*TEMP1
               C = ( TEMP - DTIIP*( DPSI+DPHI ) ) - ( D( IIM1 )-D( IIP1 ) )*( D( IIM1 )+D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               IF( DPSI.LT.TEMP1 ) THEN
                  ZZ( 3 ) = DTIIP*DTIIP*DPHI
               } else {
                  ZZ( 3 ) = DTIIP*DTIIP*( ( DPSI-TEMP1 )+DPHI )
               END IF
            } else {
               TEMP1 = Z( IIP1 ) / DTIIP
               TEMP1 = TEMP1*TEMP1
               C = ( TEMP - DTIIM*( DPSI+DPHI ) ) - ( D( IIP1 )-D( IIM1 ) )*( D( IIM1 )+D( IIP1 ) )*TEMP1
               IF( DPHI.LT.TEMP1 ) THEN
                  ZZ( 1 ) = DTIIM*DTIIM*DPSI
               } else {
                  ZZ( 1 ) = DTIIM*DTIIM*( DPSI+( DPHI-TEMP1 ) )
               END IF
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            END IF
            ZZ( 2 ) = Z( II )*Z( II )
            DD( 1 ) = DTIIM
            DD( 2 ) = DELTA( II )*WORK( II )
            DD( 3 ) = DTIIP
            CALL SLAED6( NITER, ORGATI, C, DD, ZZ, W, ETA, INFO )

            IF( INFO.NE.0 ) THEN

               // If INFO is not 0, i.e., SLAED6 failed, switch back
              t // o 2 pole interpolation.

               SWTCH3 = .FALSE.
               INFO = 0
               DTIPSQ = WORK( IP1 )*DELTA( IP1 )
               DTISQ = WORK( I )*DELTA( I )
               IF( ORGATI ) THEN
                  C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
               } else {
                  C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
               END IF
               A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
               B = DTIPSQ*DTISQ*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( ORGATI ) THEN
                        A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ*( DPSI+DPHI )
                     } else {
                        A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI)
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               } else {
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            END IF
         END IF

         // Note, eta should be positive if w is negative, and
         // eta should be negative otherwise. However,
         // if for some reason caused by roundoff, eta*w > 0,
         // we simply use one Newton step instead. This way
         // will guarantee eta*w < 0.

         IF( W*ETA.GE.ZERO ) ETA = -W / DW

         ETA = ETA / ( SIGMA+SQRT( SIGMA*SIGMA+ETA ) )
         TEMP = TAU + ETA
         IF( TEMP.GT.SGUB .OR. TEMP.LT.SGLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( SGUB-TAU ) / TWO
            } else {
               ETA = ( SGLB-TAU ) / TWO
            END IF
            IF( GEOMAVG ) THEN
               IF( W .LT. ZERO ) THEN
                  IF( TAU .GT. ZERO ) THEN
                     ETA = SQRT(SGUB*TAU)-TAU
                  END IF
               } else {
                  IF( SGLB .GT. ZERO ) THEN
                     ETA = SQRT(SGLB*TAU)-TAU
                  END IF
               END IF
            END IF
         END IF

         PREW = W

         TAU = TAU + ETA
         SIGMA = SIGMA + ETA

         DO 170 J = 1, N
            WORK( J ) = WORK( J ) + ETA
            DELTA( J ) = DELTA( J ) - ETA
  170    CONTINUE

         // Evaluate PSI and the derivative DPSI

         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 180 J = 1, IIM1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  180    CONTINUE
         ERRETM = ABS( ERRETM )

         // Evaluate PHI and the derivative DPHI

         DPHI = ZERO
         PHI = ZERO
         DO 190 J = N, IIP1, -1
            TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  190    CONTINUE

         TAU2 = WORK( II )*DELTA( II )
         TEMP = Z( II ) / TAU2
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP )
*    $          + ABS( TAU2 )*DW

         SWTCH = .FALSE.
         IF( ORGATI ) THEN
            IF( -W.GT.ABS( PREW ) / TEN ) SWTCH = .TRUE.
         } else {
            IF( W.GT.ABS( PREW ) / TEN ) SWTCH = .TRUE.
         END IF

         // Main loop to update the values of the array   DELTA and WORK

         ITER = NITER + 1

         DO 230 NITER = ITER, MAXIT

            // Test for convergence

            IF( ABS( W ).LE.EPS*ERRETM ) THEN
      // $          .OR. (SGUB-SGLB).LE.EIGHT*ABS(SGUB+SGLB) ) THEN
               GO TO 240
            END IF

            IF( W.LE.ZERO ) THEN
               SGLB = MAX( SGLB, TAU )
            } else {
               SGUB = MIN( SGUB, TAU )
            END IF

            // Calculate the new step

            IF( .NOT.SWTCH3 ) THEN
               DTIPSQ = WORK( IP1 )*DELTA( IP1 )
               DTISQ = WORK( I )*DELTA( I )
               IF( .NOT.SWTCH ) THEN
                  IF( ORGATI ) THEN
                     C = W - DTIPSQ*DW + DELSQ*( Z( I ) / DTISQ )**2
                  } else {
                     C = W - DTISQ*DW - DELSQ*( Z( IP1 ) / DTIPSQ )**2
                  END IF
               } else {
                  TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
                  IF( ORGATI ) THEN
                     DPSI = DPSI + TEMP*TEMP
                  } else {
                     DPHI = DPHI + TEMP*TEMP
                  END IF
                  C = W - DTISQ*DPSI - DTIPSQ*DPHI
               END IF
               A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
               B = DTIPSQ*DTISQ*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( .NOT.SWTCH ) THEN
                        IF( ORGATI ) THEN
                           A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ* ( DPSI+DPHI )
                        } else {
                           A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI )
                        END IF
                     } else {
                        A = DTISQ*DTISQ*DPSI + DTIPSQ*DTIPSQ*DPHI
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               } else {
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            } else {

               // Interpolation using THREE most relevant poles

               DTIIM = WORK( IIM1 )*DELTA( IIM1 )
               DTIIP = WORK( IIP1 )*DELTA( IIP1 )
               TEMP = RHOINV + PSI + PHI
               IF( SWTCH ) THEN
                  C = TEMP - DTIIM*DPSI - DTIIP*DPHI
                  ZZ( 1 ) = DTIIM*DTIIM*DPSI
                  ZZ( 3 ) = DTIIP*DTIIP*DPHI
               } else {
                  IF( ORGATI ) THEN
                     TEMP1 = Z( IIM1 ) / DTIIM
                     TEMP1 = TEMP1*TEMP1
                     TEMP2 = ( D( IIM1 )-D( IIP1 ) )* ( D( IIM1 )+D( IIP1 ) )*TEMP1
                     C = TEMP - DTIIP*( DPSI+DPHI ) - TEMP2
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     IF( DPSI.LT.TEMP1 ) THEN
                        ZZ( 3 ) = DTIIP*DTIIP*DPHI
                     } else {
                        ZZ( 3 ) = DTIIP*DTIIP*( ( DPSI-TEMP1 )+DPHI )
                     END IF
                  } else {
                     TEMP1 = Z( IIP1 ) / DTIIP
                     TEMP1 = TEMP1*TEMP1
                     TEMP2 = ( D( IIP1 )-D( IIM1 ) )* ( D( IIM1 )+D( IIP1 ) )*TEMP1
                     C = TEMP - DTIIM*( DPSI+DPHI ) - TEMP2
                     IF( DPHI.LT.TEMP1 ) THEN
                        ZZ( 1 ) = DTIIM*DTIIM*DPSI
                     } else {
                        ZZ( 1 ) = DTIIM*DTIIM*( DPSI+( DPHI-TEMP1 ) )
                     END IF
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  END IF
               END IF
               DD( 1 ) = DTIIM
               DD( 2 ) = DELTA( II )*WORK( II )
               DD( 3 ) = DTIIP
               CALL SLAED6( NITER, ORGATI, C, DD, ZZ, W, ETA, INFO )

               IF( INFO.NE.0 ) THEN

                  // If INFO is not 0, i.e., SLAED6 failed, switch
                  // back to two pole interpolation

                  SWTCH3 = .FALSE.
                  INFO = 0
                  DTIPSQ = WORK( IP1 )*DELTA( IP1 )
                  DTISQ = WORK( I )*DELTA( I )
                  IF( .NOT.SWTCH ) THEN
                     IF( ORGATI ) THEN
                        C = W - DTIPSQ*DW + DELSQ*( Z( I )/DTISQ )**2
                     } else {
                        C = W - DTISQ*DW - DELSQ*( Z( IP1 )/DTIPSQ )**2
                     END IF
                  } else {
                     TEMP = Z( II ) / ( WORK( II )*DELTA( II ) )
                     IF( ORGATI ) THEN
                        DPSI = DPSI + TEMP*TEMP
                     } else {
                        DPHI = DPHI + TEMP*TEMP
                     END IF
                     C = W - DTISQ*DPSI - DTIPSQ*DPHI
                  END IF
                  A = ( DTIPSQ+DTISQ )*W - DTIPSQ*DTISQ*DW
                  B = DTIPSQ*DTISQ*W
                  IF( C.EQ.ZERO ) THEN
                     IF( A.EQ.ZERO ) THEN
                        IF( .NOT.SWTCH ) THEN
                           IF( ORGATI ) THEN
                              A = Z( I )*Z( I ) + DTIPSQ*DTIPSQ* ( DPSI+DPHI )
                           } else {
                              A = Z( IP1 )*Z( IP1 ) + DTISQ*DTISQ*( DPSI+DPHI )
                           END IF
                        } else {
                           A = DTISQ*DTISQ*DPSI + DTIPSQ*DTIPSQ*DPHI
                        END IF
                     END IF
                     ETA = B / A
                  ELSE IF( A.LE.ZERO ) THEN
                     ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
                  } else {
                     ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
                  END IF
               END IF
            END IF

            // Note, eta should be positive if w is negative, and
            // eta should be negative otherwise. However,
            // if for some reason caused by roundoff, eta*w > 0,
            // we simply use one Newton step instead. This way
            // will guarantee eta*w < 0.

            IF( W*ETA.GE.ZERO ) ETA = -W / DW

            ETA = ETA / ( SIGMA+SQRT( SIGMA*SIGMA+ETA ) )
            TEMP=TAU+ETA
            IF( TEMP.GT.SGUB .OR. TEMP.LT.SGLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( SGUB-TAU ) / TWO
               } else {
                  ETA = ( SGLB-TAU ) / TWO
               END IF
               IF( GEOMAVG ) THEN
                  IF( W .LT. ZERO ) THEN
                     IF( TAU .GT. ZERO ) THEN
                        ETA = SQRT(SGUB*TAU)-TAU
                     END IF
                  } else {
                     IF( SGLB .GT. ZERO ) THEN
                        ETA = SQRT(SGLB*TAU)-TAU
                     END IF
                  END IF
               END IF
            END IF

            PREW = W

            TAU = TAU + ETA
            SIGMA = SIGMA + ETA

            DO 200 J = 1, N
               WORK( J ) = WORK( J ) + ETA
               DELTA( J ) = DELTA( J ) - ETA
  200       CONTINUE

            // Evaluate PSI and the derivative DPSI

            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 210 J = 1, IIM1
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
  210       CONTINUE
            ERRETM = ABS( ERRETM )

            // Evaluate PHI and the derivative DPHI

            DPHI = ZERO
            PHI = ZERO
            DO 220 J = N, IIP1, -1
               TEMP = Z( J ) / ( WORK( J )*DELTA( J ) )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
  220       CONTINUE

            TAU2 = WORK( II )*DELTA( II )
            TEMP = Z( II ) / TAU2
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP )
*    $             + ABS( TAU2 )*DW

            IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN ) SWTCH = .NOT.SWTCH

  230    CONTINUE

         // Return with INFO = 1, NITER = MAXIT and not converged

         INFO = 1

      END IF

  240 CONTINUE
      RETURN

      // End of SLASD4

      }
