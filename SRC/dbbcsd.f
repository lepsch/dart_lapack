      SUBROUTINE DBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E, B22D, B22E, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS;
      int                INFO, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      double             B11D( * ), B11E( * ), B12D( * ), B12E( * ), B21D( * ), B21E( * ), B22D( * ), B22E( * ), PHI( * ), THETA( * ), WORK( * );
      double             U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), V2T( LDV2T, * );
      // ..

*  ===================================================================

      // .. Parameters ..
      int                MAXITR;
      const              MAXITR = 6 ;
      double             HUNDRED, MEIGHTH, ONE, TEN, ZERO;
      const              HUNDRED = 100.0D0, MEIGHTH = -0.125D0, ONE = 1.0D0, TEN = 10.0D0, ZERO = 0.0D0 ;
      double             NEGONE;
      const              NEGONE = -1.0D0 ;
      double             PIOVER2;
      const     PIOVER2 = 1.57079632679489661923132169163975144210D0 ;
      // ..
      // .. Local Scalars ..

      bool               COLMAJOR, LQUERY, RESTART11, RESTART12, RESTART21, RESTART22, WANTU1, WANTU2, WANTV1T, WANTV2T;
      int                I, IMIN, IMAX, ITER, IU1CS, IU1SN, IU2CS, IU2SN, IV1TCS, IV1TSN, IV2TCS, IV2TSN, J, LWORKMIN, LWORKOPT, MAXIT, MINI;
      double             B11BULGE, B12BULGE, B21BULGE, B22BULGE, DUMMY, EPS, MU, NU, R, SIGMA11, SIGMA21, TEMP, THETAMAX, THETAMIN, THRESH, TOL, TOLMUL, UNFL, X1, X2, Y1, Y2;

      // .. External Subroutines ..
      // EXTERNAL DLASR, DSCAL, DSWAP, DLARTGP, DLARTGS, DLAS2, XERBLA
      // ..
      // .. External Functions ..
      double             DLAMCH;
      bool               LSAME;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, ATAN2, COS, MAX, MIN, SIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      LQUERY = LWORK == -1
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTV1T = LSAME( JOBV1T, 'Y' )
      WANTV2T = LSAME( JOBV2T, 'Y' )
      COLMAJOR = .NOT. LSAME( TRANS, 'T' )

      if ( M .LT. 0 ) {
         INFO = -6
      } else if ( P .LT. 0 || P .GT. M ) {
         INFO = -7
      } else if ( Q .LT. 0 || Q .GT. M ) {
         INFO = -8
      } else if ( Q .GT. P || Q .GT. M-P || Q .GT. M-Q ) {
         INFO = -8
      } else if ( WANTU1 && LDU1 .LT. P ) {
         INFO = -12
      } else if ( WANTU2 && LDU2 .LT. M-P ) {
         INFO = -14
      } else if ( WANTV1T && LDV1T .LT. Q ) {
         INFO = -16
      } else if ( WANTV2T && LDV2T .LT. M-Q ) {
         INFO = -18
      }

      // Quick return if Q = 0

      if ( INFO == 0 && Q == 0 ) {
         LWORKMIN = 1
         WORK(1) = LWORKMIN
         RETURN
      }

      // Compute workspace

      if ( INFO == 0 ) {
         IU1CS = 1
         IU1SN = IU1CS + Q
         IU2CS = IU1SN + Q
         IU2SN = IU2CS + Q
         IV1TCS = IU2SN + Q
         IV1TSN = IV1TCS + Q
         IV2TCS = IV1TSN + Q
         IV2TSN = IV2TCS + Q
         LWORKOPT = IV2TSN + Q - 1
         LWORKMIN = LWORKOPT
         WORK(1) = LWORKOPT
         if ( LWORK .LT. LWORKMIN && .NOT. LQUERY ) {
            INFO = -28
         }
      }

      if ( INFO != 0 ) {
         xerbla('DBBCSD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Get machine constants

      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      TOLMUL = MAX( TEN, MIN( HUNDRED, EPS**MEIGHTH ) )
      TOL = TOLMUL*EPS
      THRESH = MAX( TOL, MAXITR*Q*Q*UNFL )

      // Test for negligible sines or cosines

      for (I = 1; I <= Q; I++) {
         if ( THETA(I) .LT. THRESH ) {
            THETA(I) = ZERO
         } else if ( THETA(I) .GT. PIOVER2-THRESH ) {
            THETA(I) = PIOVER2
         }
      }
      for (I = 1; I <= Q-1; I++) {
         if ( PHI(I) .LT. THRESH ) {
            PHI(I) = ZERO
         } else if ( PHI(I) .GT. PIOVER2-THRESH ) {
            PHI(I) = PIOVER2
         }
      }

      // Initial deflation

      IMAX = Q
      DO WHILE( IMAX .GT. 1 )
         if ( PHI(IMAX-1) != ZERO ) {
            EXIT
         }
         IMAX = IMAX - 1
      }
      IMIN = IMAX - 1
      if ( IMIN .GT. 1 ) {
         DO WHILE( PHI(IMIN-1) != ZERO )
            IMIN = IMIN - 1
            if (IMIN .LE. 1) EXIT;
         }
      }

      // Initialize iteration counter

      MAXIT = MAXITR*Q*Q
      ITER = 0

      // Begin main iteration loop

      DO WHILE( IMAX .GT. 1 )

         // Compute the matrix entries

         B11D(IMIN) = COS( THETA(IMIN) )
         B21D(IMIN) = -SIN( THETA(IMIN) )
         for (I = IMIN; I <= IMAX - 1; I++) {
            B11E(I) = -SIN( THETA(I) ) * SIN( PHI(I) )
            B11D(I+1) = COS( THETA(I+1) ) * COS( PHI(I) )
            B12D(I) = SIN( THETA(I) ) * COS( PHI(I) )
            B12E(I) = COS( THETA(I+1) ) * SIN( PHI(I) )
            B21E(I) = -COS( THETA(I) ) * SIN( PHI(I) )
            B21D(I+1) = -SIN( THETA(I+1) ) * COS( PHI(I) )
            B22D(I) = COS( THETA(I) ) * COS( PHI(I) )
            B22E(I) = -SIN( THETA(I+1) ) * SIN( PHI(I) )
         }
         B12D(IMAX) = SIN( THETA(IMAX) )
         B22D(IMAX) = COS( THETA(IMAX) )

         // Abort if not converging; otherwise, increment ITER

         if ( ITER .GT. MAXIT ) {
            INFO = 0
            for (I = 1; I <= Q; I++) {
               IF( PHI(I) != ZERO ) INFO = INFO + 1
            }
            RETURN
         }

         ITER = ITER + IMAX - IMIN

         // Compute shifts

         THETAMAX = THETA(IMIN)
         THETAMIN = THETA(IMIN)
         for (I = IMIN+1; I <= IMAX; I++) {
            IF( THETA(I) > THETAMAX ) THETAMAX = THETA(I)             IF( THETA(I) < THETAMIN ) THETAMIN = THETA(I)
         }

         if ( THETAMAX .GT. PIOVER2 - THRESH ) {

            // Zero on diagonals of B11 and B22; induce deflation with a
            // zero shift

            MU = ZERO
            NU = ONE

         } else if ( THETAMIN .LT. THRESH ) {

            // Zero on diagonals of B12 and B22; induce deflation with a
            // zero shift

            MU = ONE
            NU = ZERO

         } else {

            // Compute shifts for B11 and B21 and use the lesser

            dlas2(B11D(IMAX-1), B11E(IMAX-1), B11D(IMAX), SIGMA11, DUMMY );
            dlas2(B21D(IMAX-1), B21E(IMAX-1), B21D(IMAX), SIGMA21, DUMMY );

            if ( SIGMA11 .LE. SIGMA21 ) {
               MU = SIGMA11
               NU = SQRT( ONE - MU**2 )
               if ( MU .LT. THRESH ) {
                  MU = ZERO
                  NU = ONE
               }
            } else {
               NU = SIGMA21
               MU = SQRT( 1.0 - NU**2 )
               if ( NU .LT. THRESH ) {
                  MU = ONE
                  NU = ZERO
               }
            }
         }

         // Rotate to produce bulges in B11 and B21

         if ( MU .LE. NU ) {
            dlartgs(B11D(IMIN), B11E(IMIN), MU, WORK(IV1TCS+IMIN-1), WORK(IV1TSN+IMIN-1) );
         } else {
            dlartgs(B21D(IMIN), B21E(IMIN), NU, WORK(IV1TCS+IMIN-1), WORK(IV1TSN+IMIN-1) );
         }

         TEMP = WORK(IV1TCS+IMIN-1)*B11D(IMIN) + WORK(IV1TSN+IMIN-1)*B11E(IMIN)          B11E(IMIN) = WORK(IV1TCS+IMIN-1)*B11E(IMIN) - WORK(IV1TSN+IMIN-1)*B11D(IMIN)
         B11D(IMIN) = TEMP
         B11BULGE = WORK(IV1TSN+IMIN-1)*B11D(IMIN+1)
         B11D(IMIN+1) = WORK(IV1TCS+IMIN-1)*B11D(IMIN+1)
         TEMP = WORK(IV1TCS+IMIN-1)*B21D(IMIN) + WORK(IV1TSN+IMIN-1)*B21E(IMIN)          B21E(IMIN) = WORK(IV1TCS+IMIN-1)*B21E(IMIN) - WORK(IV1TSN+IMIN-1)*B21D(IMIN)
         B21D(IMIN) = TEMP
         B21BULGE = WORK(IV1TSN+IMIN-1)*B21D(IMIN+1)
         B21D(IMIN+1) = WORK(IV1TCS+IMIN-1)*B21D(IMIN+1)

         // Compute THETA(IMIN)

         THETA( IMIN ) = ATAN2( SQRT( B21D(IMIN)**2+B21BULGE**2 ), SQRT( B11D(IMIN)**2+B11BULGE**2 ) )

         // Chase the bulges in B11(IMIN+1,IMIN) and B21(IMIN+1,IMIN)

         if ( B11D(IMIN)**2+B11BULGE**2 .GT. THRESH**2 ) {
            dlartgp(B11BULGE, B11D(IMIN), WORK(IU1SN+IMIN-1), WORK(IU1CS+IMIN-1), R );
         } else if ( MU .LE. NU ) {
            dlartgs(B11E( IMIN ), B11D( IMIN + 1 ), MU, WORK(IU1CS+IMIN-1), WORK(IU1SN+IMIN-1) );
         } else {
            dlartgs(B12D( IMIN ), B12E( IMIN ), NU, WORK(IU1CS+IMIN-1), WORK(IU1SN+IMIN-1) );
         }
         if ( B21D(IMIN)**2+B21BULGE**2 .GT. THRESH**2 ) {
            dlartgp(B21BULGE, B21D(IMIN), WORK(IU2SN+IMIN-1), WORK(IU2CS+IMIN-1), R );
         } else if ( NU .LT. MU ) {
            dlartgs(B21E( IMIN ), B21D( IMIN + 1 ), NU, WORK(IU2CS+IMIN-1), WORK(IU2SN+IMIN-1) );
         } else {
            dlartgs(B22D(IMIN), B22E(IMIN), MU, WORK(IU2CS+IMIN-1), WORK(IU2SN+IMIN-1) );
         }
         WORK(IU2CS+IMIN-1) = -WORK(IU2CS+IMIN-1)
         WORK(IU2SN+IMIN-1) = -WORK(IU2SN+IMIN-1)

         TEMP = WORK(IU1CS+IMIN-1)*B11E(IMIN) + WORK(IU1SN+IMIN-1)*B11D(IMIN+1)          B11D(IMIN+1) = WORK(IU1CS+IMIN-1)*B11D(IMIN+1) - WORK(IU1SN+IMIN-1)*B11E(IMIN)
         B11E(IMIN) = TEMP
         if ( IMAX .GT. IMIN+1 ) {
            B11BULGE = WORK(IU1SN+IMIN-1)*B11E(IMIN+1)
            B11E(IMIN+1) = WORK(IU1CS+IMIN-1)*B11E(IMIN+1)
         }
         TEMP = WORK(IU1CS+IMIN-1)*B12D(IMIN) + WORK(IU1SN+IMIN-1)*B12E(IMIN)          B12E(IMIN) = WORK(IU1CS+IMIN-1)*B12E(IMIN) - WORK(IU1SN+IMIN-1)*B12D(IMIN)
         B12D(IMIN) = TEMP
         B12BULGE = WORK(IU1SN+IMIN-1)*B12D(IMIN+1)
         B12D(IMIN+1) = WORK(IU1CS+IMIN-1)*B12D(IMIN+1)
         TEMP = WORK(IU2CS+IMIN-1)*B21E(IMIN) + WORK(IU2SN+IMIN-1)*B21D(IMIN+1)          B21D(IMIN+1) = WORK(IU2CS+IMIN-1)*B21D(IMIN+1) - WORK(IU2SN+IMIN-1)*B21E(IMIN)
         B21E(IMIN) = TEMP
         if ( IMAX .GT. IMIN+1 ) {
            B21BULGE = WORK(IU2SN+IMIN-1)*B21E(IMIN+1)
            B21E(IMIN+1) = WORK(IU2CS+IMIN-1)*B21E(IMIN+1)
         }
         TEMP = WORK(IU2CS+IMIN-1)*B22D(IMIN) + WORK(IU2SN+IMIN-1)*B22E(IMIN)          B22E(IMIN) = WORK(IU2CS+IMIN-1)*B22E(IMIN) - WORK(IU2SN+IMIN-1)*B22D(IMIN)
         B22D(IMIN) = TEMP
         B22BULGE = WORK(IU2SN+IMIN-1)*B22D(IMIN+1)
         B22D(IMIN+1) = WORK(IU2CS+IMIN-1)*B22D(IMIN+1)

         // Inner loop: chase bulges from B11(IMIN,IMIN+2),
         // B12(IMIN,IMIN+1), B21(IMIN,IMIN+2), and B22(IMIN,IMIN+1) to
         // bottom-right

         for (I = IMIN+1; I <= IMAX-1; I++) {

            // Compute PHI(I-1)

            X1 = SIN(THETA(I-1))*B11E(I-1) + COS(THETA(I-1))*B21E(I-1)
            X2 = SIN(THETA(I-1))*B11BULGE + COS(THETA(I-1))*B21BULGE
            Y1 = SIN(THETA(I-1))*B12D(I-1) + COS(THETA(I-1))*B22D(I-1)
            Y2 = SIN(THETA(I-1))*B12BULGE + COS(THETA(I-1))*B22BULGE

            PHI(I-1) = ATAN2( SQRT(X1**2+X2**2), SQRT(Y1**2+Y2**2) )

            // Determine if there are bulges to chase or if a new direct
            // summand has been reached

            RESTART11 = B11E(I-1)**2 + B11BULGE**2 .LE. THRESH**2
            RESTART21 = B21E(I-1)**2 + B21BULGE**2 .LE. THRESH**2
            RESTART12 = B12D(I-1)**2 + B12BULGE**2 .LE. THRESH**2
            RESTART22 = B22D(I-1)**2 + B22BULGE**2 .LE. THRESH**2

            // If possible, chase bulges from B11(I-1,I+1), B12(I-1,I),
            // B21(I-1,I+1), and B22(I-1,I). If necessary, restart bulge-
            // chasing by applying the original shift again.

            if ( .NOT. RESTART11 && .NOT. RESTART21 ) {
               dlartgp(X2, X1, WORK(IV1TSN+I-1), WORK(IV1TCS+I-1), R );
            } else if ( .NOT. RESTART11 && RESTART21 ) {
               dlartgp(B11BULGE, B11E(I-1), WORK(IV1TSN+I-1), WORK(IV1TCS+I-1), R );
            } else if ( RESTART11 && .NOT. RESTART21 ) {
               dlartgp(B21BULGE, B21E(I-1), WORK(IV1TSN+I-1), WORK(IV1TCS+I-1), R );
            } else if ( MU .LE. NU ) {
               dlartgs(B11D(I), B11E(I), MU, WORK(IV1TCS+I-1), WORK(IV1TSN+I-1) );
            } else {
               dlartgs(B21D(I), B21E(I), NU, WORK(IV1TCS+I-1), WORK(IV1TSN+I-1) );
            }
            WORK(IV1TCS+I-1) = -WORK(IV1TCS+I-1)
            WORK(IV1TSN+I-1) = -WORK(IV1TSN+I-1)
            if ( .NOT. RESTART12 && .NOT. RESTART22 ) {
               dlartgp(Y2, Y1, WORK(IV2TSN+I-1-1), WORK(IV2TCS+I-1-1), R );
            } else if ( .NOT. RESTART12 && RESTART22 ) {
               dlartgp(B12BULGE, B12D(I-1), WORK(IV2TSN+I-1-1), WORK(IV2TCS+I-1-1), R );
            } else if ( RESTART12 && .NOT. RESTART22 ) {
               dlartgp(B22BULGE, B22D(I-1), WORK(IV2TSN+I-1-1), WORK(IV2TCS+I-1-1), R );
            } else if ( NU .LT. MU ) {
               dlartgs(B12E(I-1), B12D(I), NU, WORK(IV2TCS+I-1-1), WORK(IV2TSN+I-1-1) );
            } else {
               dlartgs(B22E(I-1), B22D(I), MU, WORK(IV2TCS+I-1-1), WORK(IV2TSN+I-1-1) );
            }

            TEMP = WORK(IV1TCS+I-1)*B11D(I) + WORK(IV1TSN+I-1)*B11E(I)
            B11E(I) = WORK(IV1TCS+I-1)*B11E(I) - WORK(IV1TSN+I-1)*B11D(I)
            B11D(I) = TEMP
            B11BULGE = WORK(IV1TSN+I-1)*B11D(I+1)
            B11D(I+1) = WORK(IV1TCS+I-1)*B11D(I+1)
            TEMP = WORK(IV1TCS+I-1)*B21D(I) + WORK(IV1TSN+I-1)*B21E(I)
            B21E(I) = WORK(IV1TCS+I-1)*B21E(I) - WORK(IV1TSN+I-1)*B21D(I)
            B21D(I) = TEMP
            B21BULGE = WORK(IV1TSN+I-1)*B21D(I+1)
            B21D(I+1) = WORK(IV1TCS+I-1)*B21D(I+1)
            TEMP = WORK(IV2TCS+I-1-1)*B12E(I-1) + WORK(IV2TSN+I-1-1)*B12D(I)             B12D(I) = WORK(IV2TCS+I-1-1)*B12D(I) - WORK(IV2TSN+I-1-1)*B12E(I-1)
            B12E(I-1) = TEMP
            B12BULGE = WORK(IV2TSN+I-1-1)*B12E(I)
            B12E(I) = WORK(IV2TCS+I-1-1)*B12E(I)
            TEMP = WORK(IV2TCS+I-1-1)*B22E(I-1) + WORK(IV2TSN+I-1-1)*B22D(I)             B22D(I) = WORK(IV2TCS+I-1-1)*B22D(I) - WORK(IV2TSN+I-1-1)*B22E(I-1)
            B22E(I-1) = TEMP
            B22BULGE = WORK(IV2TSN+I-1-1)*B22E(I)
            B22E(I) = WORK(IV2TCS+I-1-1)*B22E(I)

            // Compute THETA(I)

            X1 = COS(PHI(I-1))*B11D(I) + SIN(PHI(I-1))*B12E(I-1)
            X2 = COS(PHI(I-1))*B11BULGE + SIN(PHI(I-1))*B12BULGE
            Y1 = COS(PHI(I-1))*B21D(I) + SIN(PHI(I-1))*B22E(I-1)
            Y2 = COS(PHI(I-1))*B21BULGE + SIN(PHI(I-1))*B22BULGE

            THETA(I) = ATAN2( SQRT(Y1**2+Y2**2), SQRT(X1**2+X2**2) )

            // Determine if there are bulges to chase or if a new direct
            // summand has been reached

            RESTART11 =   B11D(I)**2 + B11BULGE**2 .LE. THRESH**2
            RESTART12 = B12E(I-1)**2 + B12BULGE**2 .LE. THRESH**2
            RESTART21 =   B21D(I)**2 + B21BULGE**2 .LE. THRESH**2
            RESTART22 = B22E(I-1)**2 + B22BULGE**2 .LE. THRESH**2

            // If possible, chase bulges from B11(I+1,I), B12(I+1,I-1),
            // B21(I+1,I), and B22(I+1,I-1). If necessary, restart bulge-
            // chasing by applying the original shift again.

            if ( .NOT. RESTART11 && .NOT. RESTART12 ) {
               dlartgp(X2, X1, WORK(IU1SN+I-1), WORK(IU1CS+I-1), R );
            } else if ( .NOT. RESTART11 && RESTART12 ) {
               dlartgp(B11BULGE, B11D(I), WORK(IU1SN+I-1), WORK(IU1CS+I-1), R );
            } else if ( RESTART11 && .NOT. RESTART12 ) {
               dlartgp(B12BULGE, B12E(I-1), WORK(IU1SN+I-1), WORK(IU1CS+I-1), R );
            } else if ( MU .LE. NU ) {
               dlartgs(B11E(I), B11D(I+1), MU, WORK(IU1CS+I-1), WORK(IU1SN+I-1) );
            } else {
               dlartgs(B12D(I), B12E(I), NU, WORK(IU1CS+I-1), WORK(IU1SN+I-1) );
            }
            if ( .NOT. RESTART21 && .NOT. RESTART22 ) {
               dlartgp(Y2, Y1, WORK(IU2SN+I-1), WORK(IU2CS+I-1), R );
            } else if ( .NOT. RESTART21 && RESTART22 ) {
               dlartgp(B21BULGE, B21D(I), WORK(IU2SN+I-1), WORK(IU2CS+I-1), R );
            } else if ( RESTART21 && .NOT. RESTART22 ) {
               dlartgp(B22BULGE, B22E(I-1), WORK(IU2SN+I-1), WORK(IU2CS+I-1), R );
            } else if ( NU .LT. MU ) {
               dlartgs(B21E(I), B21E(I+1), NU, WORK(IU2CS+I-1), WORK(IU2SN+I-1) );
            } else {
               dlartgs(B22D(I), B22E(I), MU, WORK(IU2CS+I-1), WORK(IU2SN+I-1) );
            }
            WORK(IU2CS+I-1) = -WORK(IU2CS+I-1)
            WORK(IU2SN+I-1) = -WORK(IU2SN+I-1)

            TEMP = WORK(IU1CS+I-1)*B11E(I) + WORK(IU1SN+I-1)*B11D(I+1)
            B11D(I+1) = WORK(IU1CS+I-1)*B11D(I+1) - WORK(IU1SN+I-1)*B11E(I)
            B11E(I) = TEMP
            if ( I .LT. IMAX - 1 ) {
               B11BULGE = WORK(IU1SN+I-1)*B11E(I+1)
               B11E(I+1) = WORK(IU1CS+I-1)*B11E(I+1)
            }
            TEMP = WORK(IU2CS+I-1)*B21E(I) + WORK(IU2SN+I-1)*B21D(I+1)
            B21D(I+1) = WORK(IU2CS+I-1)*B21D(I+1) - WORK(IU2SN+I-1)*B21E(I)
            B21E(I) = TEMP
            if ( I .LT. IMAX - 1 ) {
               B21BULGE = WORK(IU2SN+I-1)*B21E(I+1)
               B21E(I+1) = WORK(IU2CS+I-1)*B21E(I+1)
            }
            TEMP = WORK(IU1CS+I-1)*B12D(I) + WORK(IU1SN+I-1)*B12E(I)
            B12E(I) = WORK(IU1CS+I-1)*B12E(I) - WORK(IU1SN+I-1)*B12D(I)
            B12D(I) = TEMP
            B12BULGE = WORK(IU1SN+I-1)*B12D(I+1)
            B12D(I+1) = WORK(IU1CS+I-1)*B12D(I+1)
            TEMP = WORK(IU2CS+I-1)*B22D(I) + WORK(IU2SN+I-1)*B22E(I)
            B22E(I) = WORK(IU2CS+I-1)*B22E(I) - WORK(IU2SN+I-1)*B22D(I)
            B22D(I) = TEMP
            B22BULGE = WORK(IU2SN+I-1)*B22D(I+1)
            B22D(I+1) = WORK(IU2CS+I-1)*B22D(I+1)

         }

         // Compute PHI(IMAX-1)

         X1 = SIN(THETA(IMAX-1))*B11E(IMAX-1) + COS(THETA(IMAX-1))*B21E(IMAX-1)          Y1 = SIN(THETA(IMAX-1))*B12D(IMAX-1) + COS(THETA(IMAX-1))*B22D(IMAX-1)
         Y2 = SIN(THETA(IMAX-1))*B12BULGE + COS(THETA(IMAX-1))*B22BULGE

         PHI(IMAX-1) = ATAN2( ABS(X1), SQRT(Y1**2+Y2**2) )

         // Chase bulges from B12(IMAX-1,IMAX) and B22(IMAX-1,IMAX)

         RESTART12 = B12D(IMAX-1)**2 + B12BULGE**2 .LE. THRESH**2
         RESTART22 = B22D(IMAX-1)**2 + B22BULGE**2 .LE. THRESH**2

         if ( .NOT. RESTART12 && .NOT. RESTART22 ) {
            dlartgp(Y2, Y1, WORK(IV2TSN+IMAX-1-1), WORK(IV2TCS+IMAX-1-1), R );
         } else if ( .NOT. RESTART12 && RESTART22 ) {
            dlartgp(B12BULGE, B12D(IMAX-1), WORK(IV2TSN+IMAX-1-1), WORK(IV2TCS+IMAX-1-1), R );
         } else if ( RESTART12 && .NOT. RESTART22 ) {
            dlartgp(B22BULGE, B22D(IMAX-1), WORK(IV2TSN+IMAX-1-1), WORK(IV2TCS+IMAX-1-1), R );
         } else if ( NU .LT. MU ) {
            dlartgs(B12E(IMAX-1), B12D(IMAX), NU, WORK(IV2TCS+IMAX-1-1), WORK(IV2TSN+IMAX-1-1) );
         } else {
            dlartgs(B22E(IMAX-1), B22D(IMAX), MU, WORK(IV2TCS+IMAX-1-1), WORK(IV2TSN+IMAX-1-1) );
         }

         TEMP = WORK(IV2TCS+IMAX-1-1)*B12E(IMAX-1) + WORK(IV2TSN+IMAX-1-1)*B12D(IMAX)          B12D(IMAX) = WORK(IV2TCS+IMAX-1-1)*B12D(IMAX) - WORK(IV2TSN+IMAX-1-1)*B12E(IMAX-1)
         B12E(IMAX-1) = TEMP
         TEMP = WORK(IV2TCS+IMAX-1-1)*B22E(IMAX-1) + WORK(IV2TSN+IMAX-1-1)*B22D(IMAX)          B22D(IMAX) = WORK(IV2TCS+IMAX-1-1)*B22D(IMAX) - WORK(IV2TSN+IMAX-1-1)*B22E(IMAX-1)
         B22E(IMAX-1) = TEMP

         // Update singular vectors

         if ( WANTU1 ) {
            if ( COLMAJOR ) {
               dlasr('R', 'V', 'F', P, IMAX-IMIN+1, WORK(IU1CS+IMIN-1), WORK(IU1SN+IMIN-1), U1(1,IMIN), LDU1 );
            } else {
               dlasr('L', 'V', 'F', IMAX-IMIN+1, P, WORK(IU1CS+IMIN-1), WORK(IU1SN+IMIN-1), U1(IMIN,1), LDU1 );
            }
         }
         if ( WANTU2 ) {
            if ( COLMAJOR ) {
               dlasr('R', 'V', 'F', M-P, IMAX-IMIN+1, WORK(IU2CS+IMIN-1), WORK(IU2SN+IMIN-1), U2(1,IMIN), LDU2 );
            } else {
               dlasr('L', 'V', 'F', IMAX-IMIN+1, M-P, WORK(IU2CS+IMIN-1), WORK(IU2SN+IMIN-1), U2(IMIN,1), LDU2 );
            }
         }
         if ( WANTV1T ) {
            if ( COLMAJOR ) {
               dlasr('L', 'V', 'F', IMAX-IMIN+1, Q, WORK(IV1TCS+IMIN-1), WORK(IV1TSN+IMIN-1), V1T(IMIN,1), LDV1T );
            } else {
               dlasr('R', 'V', 'F', Q, IMAX-IMIN+1, WORK(IV1TCS+IMIN-1), WORK(IV1TSN+IMIN-1), V1T(1,IMIN), LDV1T );
            }
         }
         if ( WANTV2T ) {
            if ( COLMAJOR ) {
               dlasr('L', 'V', 'F', IMAX-IMIN+1, M-Q, WORK(IV2TCS+IMIN-1), WORK(IV2TSN+IMIN-1), V2T(IMIN,1), LDV2T );
            } else {
               dlasr('R', 'V', 'F', M-Q, IMAX-IMIN+1, WORK(IV2TCS+IMIN-1), WORK(IV2TSN+IMIN-1), V2T(1,IMIN), LDV2T );
            }
         }

         // Fix signs on B11(IMAX-1,IMAX) and B21(IMAX-1,IMAX)

         if ( B11E(IMAX-1)+B21E(IMAX-1) .GT. 0 ) {
            B11D(IMAX) = -B11D(IMAX)
            B21D(IMAX) = -B21D(IMAX)
            if ( WANTV1T ) {
               if ( COLMAJOR ) {
                  dscal(Q, NEGONE, V1T(IMAX,1), LDV1T );
               } else {
                  dscal(Q, NEGONE, V1T(1,IMAX), 1 );
               }
            }
         }

         // Compute THETA(IMAX)

         X1 = COS(PHI(IMAX-1))*B11D(IMAX) + SIN(PHI(IMAX-1))*B12E(IMAX-1)          Y1 = COS(PHI(IMAX-1))*B21D(IMAX) + SIN(PHI(IMAX-1))*B22E(IMAX-1)

         THETA(IMAX) = ATAN2( ABS(Y1), ABS(X1) )

         // Fix signs on B11(IMAX,IMAX), B12(IMAX,IMAX-1), B21(IMAX,IMAX),
         // and B22(IMAX,IMAX-1)

         if ( B11D(IMAX)+B12E(IMAX-1) .LT. 0 ) {
            B12D(IMAX) = -B12D(IMAX)
            if ( WANTU1 ) {
               if ( COLMAJOR ) {
                  dscal(P, NEGONE, U1(1,IMAX), 1 );
               } else {
                  dscal(P, NEGONE, U1(IMAX,1), LDU1 );
               }
            }
         }
         if ( B21D(IMAX)+B22E(IMAX-1) .GT. 0 ) {
            B22D(IMAX) = -B22D(IMAX)
            if ( WANTU2 ) {
               if ( COLMAJOR ) {
                  dscal(M-P, NEGONE, U2(1,IMAX), 1 );
               } else {
                  dscal(M-P, NEGONE, U2(IMAX,1), LDU2 );
               }
            }
         }

         // Fix signs on B12(IMAX,IMAX) and B22(IMAX,IMAX)

         if ( B12D(IMAX)+B22D(IMAX) .LT. 0 ) {
            if ( WANTV2T ) {
               if ( COLMAJOR ) {
                  dscal(M-Q, NEGONE, V2T(IMAX,1), LDV2T );
               } else {
                  dscal(M-Q, NEGONE, V2T(1,IMAX), 1 );
               }
            }
         }

         // Test for negligible sines or cosines

         for (I = IMIN; I <= IMAX; I++) {
            if ( THETA(I) .LT. THRESH ) {
               THETA(I) = ZERO
            } else if ( THETA(I) .GT. PIOVER2-THRESH ) {
               THETA(I) = PIOVER2
            }
         }
         for (I = IMIN; I <= IMAX-1; I++) {
            if ( PHI(I) .LT. THRESH ) {
               PHI(I) = ZERO
            } else if ( PHI(I) .GT. PIOVER2-THRESH ) {
               PHI(I) = PIOVER2
            }
         }

         // Deflate

         if (IMAX .GT. 1) {
            DO WHILE( PHI(IMAX-1) == ZERO )
               IMAX = IMAX - 1
               if (IMAX .LE. 1) EXIT;
            }
         }
         if (IMIN .GT. IMAX - 1) IMIN = IMAX - 1;
         if (IMIN .GT. 1) {
            DO WHILE (PHI(IMIN-1) != ZERO)
                IMIN = IMIN - 1
                if (IMIN .LE. 1) EXIT;
            }
         }

         // Repeat main iteration loop

      }

      // Postprocessing: order THETA from least to greatest

      for (I = 1; I <= Q; I++) {

         MINI = I
         THETAMIN = THETA(I)
         for (J = I+1; J <= Q; J++) {
            if ( THETA(J) .LT. THETAMIN ) {
               MINI = J
               THETAMIN = THETA(J)
            }
         }

         if ( MINI != I ) {
            THETA(MINI) = THETA(I)
            THETA(I) = THETAMIN
            if ( COLMAJOR ) {
               if (WANTU1) CALL DSWAP( P, U1(1,I), 1, U1(1,MINI), 1 )                IF( WANTU2 ) CALL DSWAP( M-P, U2(1,I), 1, U2(1,MINI), 1 )                IF( WANTV1T ) CALL DSWAP( Q, V1T(I,1), LDV1T, V1T(MINI,1), LDV1T )                IF( WANTV2T ) CALL DSWAP( M-Q, V2T(I,1), LDV2T, V2T(MINI,1), LDV2T );
            } else {
               if (WANTU1) CALL DSWAP( P, U1(I,1), LDU1, U1(MINI,1), LDU1 )                IF( WANTU2 ) CALL DSWAP( M-P, U2(I,1), LDU2, U2(MINI,1), LDU2 )                IF( WANTV1T ) CALL DSWAP( Q, V1T(1,I), 1, V1T(1,MINI), 1 )                IF( WANTV2T ) CALL DSWAP( M-Q, V2T(1,I), 1, V2T(1,MINI), 1 );
            }
         }

      }

      RETURN

      // End of DBBCSD

      }
