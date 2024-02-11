      void cbbcsd(final int JOBU1, final int JOBU2, final int JOBV1T, final int JOBV2T, final int TRANS, final int M, final int P, final int Q, final int THETA, final int PHI, final Matrix<double> U1, final int LDU1, final Matrix<double> U2, final int LDU2, final Matrix<double> V1T, final int LDV1T, final Matrix<double> V2T, final int LDV2T, final int B11D, final int B11E, final int B12D, final int B12E, final int B21D, final int B21E, final int B22D, final int B22E, final Array<int> RWORK, final int LRWORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS;
      int                INFO, LDU1, LDU2, LDV1T, LDV2T, LRWORK, M, P, Q;
      double               B11D( * ), B11E( * ), B12D( * ), B12E( * ), B21D( * ), B21E( * ), B22D( * ), B22E( * ), PHI( * ), THETA( * ), RWORK( * );
      Complex            U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), V2T( LDV2T, * );
      // ..

// ===================================================================

      // .. Parameters ..
      int                MAXITR;
      const              MAXITR = 6 ;
      double               HUNDRED, MEIGHTH, ONE, TEN, ZERO;
      const              HUNDRED = 100.0, MEIGHTH = -0.125, ONE = 1.0, TEN = 10.0, ZERO = 0.0 ;
      Complex            NEGONECOMPLEX;
      const              NEGONECOMPLEX = (-1.0,0.0) ;
      double               PIOVER2;
      const     PIOVER2 = 1.57079632679489661923132169163975144210 ;
      bool               COLMAJOR, LQUERY, RESTART11, RESTART12, RESTART21, RESTART22, WANTU1, WANTU2, WANTV1T, WANTV2T;
      int                I, IMIN, IMAX, ITER, IU1CS, IU1SN, IU2CS, IU2SN, IV1TCS, IV1TSN, IV2TCS, IV2TSN, J, LRWORKMIN, LRWORKOPT, MAXIT, MINI;
      double               B11BULGE, B12BULGE, B21BULGE, B22BULGE, DUMMY, EPS, MU, NU, R, SIGMA11, SIGMA21, TEMP, THETAMAX, THETAMIN, THRESH, TOL, TOLMUL, UNFL, X1, X2, Y1, Y2;

      // .. External Subroutines ..
      // EXTERNAL CLASR, CSCAL, CSWAP, SLARTGP, SLARTGS, SLAS2, XERBLA
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      //- bool               lsame;
      // EXTERNAL lsame, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, ATAN2, COS, MAX, MIN, SIN, SQRT

      // Test input arguments

      INFO = 0;
      LQUERY = LRWORK == -1;
      WANTU1 = lsame( JOBU1, 'Y' );
      WANTU2 = lsame( JOBU2, 'Y' );
      WANTV1T = lsame( JOBV1T, 'Y' );
      WANTV2T = lsame( JOBV2T, 'Y' );
      COLMAJOR = !lsame( TRANS, 'T' );

      if ( M < 0 ) {
         INFO = -6;
      } else if ( P < 0 || P > M ) {
         INFO = -7;
      } else if ( Q < 0 || Q > M ) {
         INFO = -8;
      } else if ( Q > P || Q > M-P || Q > M-Q ) {
         INFO = -8;
      } else if ( WANTU1 && LDU1 < P ) {
         INFO = -12;
      } else if ( WANTU2 && LDU2 < M-P ) {
         INFO = -14;
      } else if ( WANTV1T && LDV1T < Q ) {
         INFO = -16;
      } else if ( WANTV2T && LDV2T < M-Q ) {
         INFO = -18;
      }

      // Quick return if Q = 0

      if ( INFO == 0 && Q == 0 ) {
         LRWORKMIN = 1;
         RWORK[1] = LRWORKMIN;
         return;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         IU1CS = 1;
         IU1SN = IU1CS + Q;
         IU2CS = IU1SN + Q;
         IU2SN = IU2CS + Q;
         IV1TCS = IU2SN + Q;
         IV1TSN = IV1TCS + Q;
         IV2TCS = IV1TSN + Q;
         IV2TSN = IV2TCS + Q;
         LRWORKOPT = IV2TSN + Q - 1;
         LRWORKMIN = LRWORKOPT;
         RWORK[1] = LRWORKOPT;
         if ( LRWORK < LRWORKMIN && !LQUERY ) {
            INFO = -28;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CBBCSD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Get machine constants

      EPS = SLAMCH( 'Epsilon' );
      UNFL = SLAMCH( 'Safe minimum' );
      TOLMUL = max( TEN, min( HUNDRED, EPS**MEIGHTH ) );
      TOL = TOLMUL*EPS;
      THRESH = max( TOL, MAXITR*Q*Q*UNFL );

      // Test for negligible sines or cosines

      for (I = 1; I <= Q; I++) {
         if ( THETA(I) < THRESH ) {
            THETA[I] = ZERO;
         } else if ( THETA(I) > PIOVER2-THRESH ) {
            THETA[I] = PIOVER2;
         }
      }
      for (I = 1; I <= Q-1; I++) {
         if ( PHI(I) < THRESH ) {
            PHI[I] = ZERO;
         } else if ( PHI(I) > PIOVER2-THRESH ) {
            PHI[I] = PIOVER2;
         }
      }

      // Initial deflation

      IMAX = Q;
      while (IMAX > 1) {
         if ( PHI(IMAX-1) != ZERO ) {
            break;
         }
         IMAX = IMAX - 1;
      }
      IMIN = IMAX - 1;
      if ( IMIN > 1 ) {
         while (PHI(IMIN-1) != ZERO) {
            IMIN = IMIN - 1;
            if (IMIN <= 1) break;
         }
      }

      // Initialize iteration counter

      MAXIT = MAXITR*Q*Q;
      ITER = 0;

      // Begin main iteration loop

      while (IMAX > 1) {

         // Compute the matrix entries

         B11D[IMIN] = COS( THETA(IMIN) );
         B21D[IMIN] = -SIN( THETA(IMIN) );
         for (I = IMIN; I <= IMAX - 1; I++) {
            B11E[I] = -SIN( THETA(I) ) * SIN( PHI(I) );
            B11D[I+1] = COS( THETA(I+1) ) * COS( PHI(I) );
            B12D[I] = SIN( THETA(I) ) * COS( PHI(I) );
            B12E[I] = COS( THETA(I+1) ) * SIN( PHI(I) );
            B21E[I] = -COS( THETA(I) ) * SIN( PHI(I) );
            B21D[I+1] = -SIN( THETA(I+1) ) * COS( PHI(I) );
            B22D[I] = COS( THETA(I) ) * COS( PHI(I) );
            B22E[I] = -SIN( THETA(I+1) ) * SIN( PHI(I) );
         }
         B12D[IMAX] = SIN( THETA(IMAX) );
         B22D[IMAX] = COS( THETA(IMAX) );

         // Abort if not converging; otherwise, increment ITER

         if ( ITER > MAXIT ) {
            INFO = 0;
            for (I = 1; I <= Q; I++) {
               if( PHI(I) != ZERO ) INFO = INFO + 1;
            }
            return;
         }

         ITER = ITER + IMAX - IMIN;

         // Compute shifts

         THETAMAX = THETA(IMIN);
         THETAMIN = THETA(IMIN);
         for (I = IMIN+1; I <= IMAX; I++) {
            if( THETA(I) > THETAMAX ) THETAMAX = THETA(I);
            IF( THETA(I) < THETAMIN ) THETAMIN = THETA(I);
         }

         if ( THETAMAX > PIOVER2 - THRESH ) {

            // Zero on diagonals of B11 and B22; induce deflation with a
            // zero shift

            MU = ZERO;
            NU = ONE;

         } else if ( THETAMIN < THRESH ) {

            // Zero on diagonals of B12 and B22; induce deflation with a
            // zero shift

            MU = ONE;
            NU = ZERO;

         } else {

            // Compute shifts for B11 and B21 and use the lesser

            slas2(B11D(IMAX-1), B11E(IMAX-1), B11D(IMAX), SIGMA11, DUMMY );
            slas2(B21D(IMAX-1), B21E(IMAX-1), B21D(IMAX), SIGMA21, DUMMY );

            if ( SIGMA11 <= SIGMA21 ) {
               MU = SIGMA11;
               NU = sqrt( ONE - MU**2 );
               if ( MU < THRESH ) {
                  MU = ZERO;
                  NU = ONE;
               }
            } else {
               NU = SIGMA21;
               MU = sqrt( 1.0 - NU**2 );
               if ( NU < THRESH ) {
                  MU = ONE;
                  NU = ZERO;
               }
            }
         }

         // Rotate to produce bulges in B11 and B21

         if ( MU <= NU ) {
            slartgs(B11D(IMIN), B11E(IMIN), MU, RWORK(IV1TCS+IMIN-1), RWORK(IV1TSN+IMIN-1) );
         } else {
            slartgs(B21D(IMIN), B21E(IMIN), NU, RWORK(IV1TCS+IMIN-1), RWORK(IV1TSN+IMIN-1) );
         }

         TEMP = RWORK(IV1TCS+IMIN-1)*B11D(IMIN) + RWORK(IV1TSN+IMIN-1)*B11E(IMIN)          B11E(IMIN) = RWORK(IV1TCS+IMIN-1)*B11E(IMIN) - RWORK(IV1TSN+IMIN-1)*B11D(IMIN);
         B11D[IMIN] = TEMP;
         B11BULGE = RWORK(IV1TSN+IMIN-1)*B11D(IMIN+1);
         B11D[IMIN+1] = RWORK(IV1TCS+IMIN-1)*B11D(IMIN+1);
         TEMP = RWORK(IV1TCS+IMIN-1)*B21D(IMIN) + RWORK(IV1TSN+IMIN-1)*B21E(IMIN)          B21E(IMIN) = RWORK(IV1TCS+IMIN-1)*B21E(IMIN) - RWORK(IV1TSN+IMIN-1)*B21D(IMIN);
         B21D[IMIN] = TEMP;
         B21BULGE = RWORK(IV1TSN+IMIN-1)*B21D(IMIN+1);
         B21D[IMIN+1] = RWORK(IV1TCS+IMIN-1)*B21D(IMIN+1);

         // Compute THETA(IMIN)

         THETA[IMIN] = ATAN2( sqrt( B21D(IMIN)**2+B21BULGE**2 ), sqrt( B11D(IMIN)**2+B11BULGE**2 ) );

         // Chase the bulges in B11(IMIN+1,IMIN) and B21(IMIN+1,IMIN)

         if ( B11D(IMIN)**2+B11BULGE**2 > THRESH**2 ) {
            slartgp(B11BULGE, B11D(IMIN), RWORK(IU1SN+IMIN-1), RWORK(IU1CS+IMIN-1), R );
         } else if ( MU <= NU ) {
            slartgs(B11E( IMIN ), B11D( IMIN + 1 ), MU, RWORK(IU1CS+IMIN-1), RWORK(IU1SN+IMIN-1) );
         } else {
            slartgs(B12D( IMIN ), B12E( IMIN ), NU, RWORK(IU1CS+IMIN-1), RWORK(IU1SN+IMIN-1) );
         }
         if ( B21D(IMIN)**2+B21BULGE**2 > THRESH**2 ) {
            slartgp(B21BULGE, B21D(IMIN), RWORK(IU2SN+IMIN-1), RWORK(IU2CS+IMIN-1), R );
         } else if ( NU < MU ) {
            slartgs(B21E( IMIN ), B21D( IMIN + 1 ), NU, RWORK(IU2CS+IMIN-1), RWORK(IU2SN+IMIN-1) );
         } else {
            slartgs(B22D(IMIN), B22E(IMIN), MU, RWORK(IU2CS+IMIN-1), RWORK(IU2SN+IMIN-1) );
         }
         RWORK[IU2CS+IMIN-1] = -RWORK(IU2CS+IMIN-1);
         RWORK[IU2SN+IMIN-1] = -RWORK(IU2SN+IMIN-1);

         TEMP = RWORK(IU1CS+IMIN-1)*B11E(IMIN) + RWORK(IU1SN+IMIN-1)*B11D(IMIN+1)          B11D(IMIN+1) = RWORK(IU1CS+IMIN-1)*B11D(IMIN+1) - RWORK(IU1SN+IMIN-1)*B11E(IMIN);
         B11E[IMIN] = TEMP;
         if ( IMAX > IMIN+1 ) {
            B11BULGE = RWORK(IU1SN+IMIN-1)*B11E(IMIN+1);
            B11E[IMIN+1] = RWORK(IU1CS+IMIN-1)*B11E(IMIN+1);
         }
         TEMP = RWORK(IU1CS+IMIN-1)*B12D(IMIN) + RWORK(IU1SN+IMIN-1)*B12E(IMIN)          B12E(IMIN) = RWORK(IU1CS+IMIN-1)*B12E(IMIN) - RWORK(IU1SN+IMIN-1)*B12D(IMIN);
         B12D[IMIN] = TEMP;
         B12BULGE = RWORK(IU1SN+IMIN-1)*B12D(IMIN+1);
         B12D[IMIN+1] = RWORK(IU1CS+IMIN-1)*B12D(IMIN+1);
         TEMP = RWORK(IU2CS+IMIN-1)*B21E(IMIN) + RWORK(IU2SN+IMIN-1)*B21D(IMIN+1)          B21D(IMIN+1) = RWORK(IU2CS+IMIN-1)*B21D(IMIN+1) - RWORK(IU2SN+IMIN-1)*B21E(IMIN);
         B21E[IMIN] = TEMP;
         if ( IMAX > IMIN+1 ) {
            B21BULGE = RWORK(IU2SN+IMIN-1)*B21E(IMIN+1);
            B21E[IMIN+1] = RWORK(IU2CS+IMIN-1)*B21E(IMIN+1);
         }
         TEMP = RWORK(IU2CS+IMIN-1)*B22D(IMIN) + RWORK(IU2SN+IMIN-1)*B22E(IMIN)          B22E(IMIN) = RWORK(IU2CS+IMIN-1)*B22E(IMIN) - RWORK(IU2SN+IMIN-1)*B22D(IMIN);
         B22D[IMIN] = TEMP;
         B22BULGE = RWORK(IU2SN+IMIN-1)*B22D(IMIN+1);
         B22D[IMIN+1] = RWORK(IU2CS+IMIN-1)*B22D(IMIN+1);

         // Inner loop: chase bulges from B11(IMIN,IMIN+2),
         // B12(IMIN,IMIN+1), B21(IMIN,IMIN+2), and B22(IMIN,IMIN+1) to
         // bottom-right

         for (I = IMIN+1; I <= IMAX-1; I++) {

            // Compute PHI(I-1)

            X1 = SIN(THETA(I-1))*B11E(I-1) + COS(THETA(I-1))*B21E(I-1);
            X2 = SIN(THETA(I-1))*B11BULGE + COS(THETA(I-1))*B21BULGE;
            Y1 = SIN(THETA(I-1))*B12D(I-1) + COS(THETA(I-1))*B22D(I-1);
            Y2 = SIN(THETA(I-1))*B12BULGE + COS(THETA(I-1))*B22BULGE;

            PHI[I-1] = ATAN2( sqrt(X1**2+X2**2), sqrt(Y1**2+Y2**2) );

            // Determine if there are bulges to chase or if a new direct
            // summand has been reached

            RESTART11 = B11E(I-1)**2 + B11BULGE**2 <= THRESH**2;
            RESTART21 = B21E(I-1)**2 + B21BULGE**2 <= THRESH**2;
            RESTART12 = B12D(I-1)**2 + B12BULGE**2 <= THRESH**2;
            RESTART22 = B22D(I-1)**2 + B22BULGE**2 <= THRESH**2;

            // If possible, chase bulges from B11(I-1,I+1), B12(I-1,I),
            // B21(I-1,I+1), and B22(I-1,I). If necessary, restart bulge-
            // chasing by applying the original shift again.

            if ( !RESTART11 && !RESTART21 ) {
               slartgp(X2, X1, RWORK(IV1TSN+I-1), RWORK(IV1TCS+I-1), R );
            } else if ( !RESTART11 && RESTART21 ) {
               slartgp(B11BULGE, B11E(I-1), RWORK(IV1TSN+I-1), RWORK(IV1TCS+I-1), R );
            } else if ( RESTART11 && !RESTART21 ) {
               slartgp(B21BULGE, B21E(I-1), RWORK(IV1TSN+I-1), RWORK(IV1TCS+I-1), R );
            } else if ( MU <= NU ) {
               slartgs(B11D(I), B11E(I), MU, RWORK(IV1TCS+I-1), RWORK(IV1TSN+I-1) );
            } else {
               slartgs(B21D(I), B21E(I), NU, RWORK(IV1TCS+I-1), RWORK(IV1TSN+I-1) );
            }
            RWORK[IV1TCS+I-1] = -RWORK(IV1TCS+I-1);
            RWORK[IV1TSN+I-1] = -RWORK(IV1TSN+I-1);
            if ( !RESTART12 && !RESTART22 ) {
               slartgp(Y2, Y1, RWORK(IV2TSN+I-1-1), RWORK(IV2TCS+I-1-1), R );
            } else if ( !RESTART12 && RESTART22 ) {
               slartgp(B12BULGE, B12D(I-1), RWORK(IV2TSN+I-1-1), RWORK(IV2TCS+I-1-1), R );
            } else if ( RESTART12 && !RESTART22 ) {
               slartgp(B22BULGE, B22D(I-1), RWORK(IV2TSN+I-1-1), RWORK(IV2TCS+I-1-1), R );
            } else if ( NU < MU ) {
               slartgs(B12E(I-1), B12D(I), NU, RWORK(IV2TCS+I-1-1), RWORK(IV2TSN+I-1-1) );
            } else {
               slartgs(B22E(I-1), B22D(I), MU, RWORK(IV2TCS+I-1-1), RWORK(IV2TSN+I-1-1) );
            }

            TEMP = RWORK(IV1TCS+I-1)*B11D(I) + RWORK(IV1TSN+I-1)*B11E(I);
            B11E[I] = RWORK(IV1TCS+I-1)*B11E(I) - RWORK(IV1TSN+I-1)*B11D(I);
            B11D[I] = TEMP;
            B11BULGE = RWORK(IV1TSN+I-1)*B11D(I+1);
            B11D[I+1] = RWORK(IV1TCS+I-1)*B11D(I+1);
            TEMP = RWORK(IV1TCS+I-1)*B21D(I) + RWORK(IV1TSN+I-1)*B21E(I);
            B21E[I] = RWORK(IV1TCS+I-1)*B21E(I) - RWORK(IV1TSN+I-1)*B21D(I);
            B21D[I] = TEMP;
            B21BULGE = RWORK(IV1TSN+I-1)*B21D(I+1);
            B21D[I+1] = RWORK(IV1TCS+I-1)*B21D(I+1);
            TEMP = RWORK(IV2TCS+I-1-1)*B12E(I-1) + RWORK(IV2TSN+I-1-1)*B12D(I)             B12D(I) = RWORK(IV2TCS+I-1-1)*B12D(I) - RWORK(IV2TSN+I-1-1)*B12E(I-1);
            B12E[I-1] = TEMP;
            B12BULGE = RWORK(IV2TSN+I-1-1)*B12E(I);
            B12E[I] = RWORK(IV2TCS+I-1-1)*B12E(I);
            TEMP = RWORK(IV2TCS+I-1-1)*B22E(I-1) + RWORK(IV2TSN+I-1-1)*B22D(I)             B22D(I) = RWORK(IV2TCS+I-1-1)*B22D(I) - RWORK(IV2TSN+I-1-1)*B22E(I-1);
            B22E[I-1] = TEMP;
            B22BULGE = RWORK(IV2TSN+I-1-1)*B22E(I);
            B22E[I] = RWORK(IV2TCS+I-1-1)*B22E(I);

            // Compute THETA(I)

            X1 = COS(PHI(I-1))*B11D(I) + SIN(PHI(I-1))*B12E(I-1);
            X2 = COS(PHI(I-1))*B11BULGE + SIN(PHI(I-1))*B12BULGE;
            Y1 = COS(PHI(I-1))*B21D(I) + SIN(PHI(I-1))*B22E(I-1);
            Y2 = COS(PHI(I-1))*B21BULGE + SIN(PHI(I-1))*B22BULGE;

            THETA[I] = ATAN2( sqrt(Y1**2+Y2**2), sqrt(X1**2+X2**2) );

            // Determine if there are bulges to chase or if a new direct
            // summand has been reached

            RESTART11 =   B11D(I)**2 + B11BULGE**2 <= THRESH**2;
            RESTART12 = B12E(I-1)**2 + B12BULGE**2 <= THRESH**2;
            RESTART21 =   B21D(I)**2 + B21BULGE**2 <= THRESH**2;
            RESTART22 = B22E(I-1)**2 + B22BULGE**2 <= THRESH**2;

            // If possible, chase bulges from B11(I+1,I), B12(I+1,I-1),
            // B21(I+1,I), and B22(I+1,I-1). If necessary, restart bulge-
            // chasing by applying the original shift again.

            if ( !RESTART11 && !RESTART12 ) {
               slartgp(X2, X1, RWORK(IU1SN+I-1), RWORK(IU1CS+I-1), R );
            } else if ( !RESTART11 && RESTART12 ) {
               slartgp(B11BULGE, B11D(I), RWORK(IU1SN+I-1), RWORK(IU1CS+I-1), R );
            } else if ( RESTART11 && !RESTART12 ) {
               slartgp(B12BULGE, B12E(I-1), RWORK(IU1SN+I-1), RWORK(IU1CS+I-1), R );
            } else if ( MU <= NU ) {
               slartgs(B11E(I), B11D(I+1), MU, RWORK(IU1CS+I-1), RWORK(IU1SN+I-1) );
            } else {
               slartgs(B12D(I), B12E(I), NU, RWORK(IU1CS+I-1), RWORK(IU1SN+I-1) );
            }
            if ( !RESTART21 && !RESTART22 ) {
               slartgp(Y2, Y1, RWORK(IU2SN+I-1), RWORK(IU2CS+I-1), R );
            } else if ( !RESTART21 && RESTART22 ) {
               slartgp(B21BULGE, B21D(I), RWORK(IU2SN+I-1), RWORK(IU2CS+I-1), R );
            } else if ( RESTART21 && !RESTART22 ) {
               slartgp(B22BULGE, B22E(I-1), RWORK(IU2SN+I-1), RWORK(IU2CS+I-1), R );
            } else if ( NU < MU ) {
               slartgs(B21E(I), B21E(I+1), NU, RWORK(IU2CS+I-1), RWORK(IU2SN+I-1) );
            } else {
               slartgs(B22D(I), B22E(I), MU, RWORK(IU2CS+I-1), RWORK(IU2SN+I-1) );
            }
            RWORK[IU2CS+I-1] = -RWORK(IU2CS+I-1);
            RWORK[IU2SN+I-1] = -RWORK(IU2SN+I-1);

            TEMP = RWORK(IU1CS+I-1)*B11E(I) + RWORK(IU1SN+I-1)*B11D(I+1);
            B11D[I+1] = RWORK(IU1CS+I-1)*B11D(I+1) - RWORK(IU1SN+I-1)*B11E(I);
            B11E[I] = TEMP;
            if ( I < IMAX - 1 ) {
               B11BULGE = RWORK(IU1SN+I-1)*B11E(I+1);
               B11E[I+1] = RWORK(IU1CS+I-1)*B11E(I+1);
            }
            TEMP = RWORK(IU2CS+I-1)*B21E(I) + RWORK(IU2SN+I-1)*B21D(I+1);
            B21D[I+1] = RWORK(IU2CS+I-1)*B21D(I+1) - RWORK(IU2SN+I-1)*B21E(I);
            B21E[I] = TEMP;
            if ( I < IMAX - 1 ) {
               B21BULGE = RWORK(IU2SN+I-1)*B21E(I+1);
               B21E[I+1] = RWORK(IU2CS+I-1)*B21E(I+1);
            }
            TEMP = RWORK(IU1CS+I-1)*B12D(I) + RWORK(IU1SN+I-1)*B12E(I);
            B12E[I] = RWORK(IU1CS+I-1)*B12E(I) - RWORK(IU1SN+I-1)*B12D(I);
            B12D[I] = TEMP;
            B12BULGE = RWORK(IU1SN+I-1)*B12D(I+1);
            B12D[I+1] = RWORK(IU1CS+I-1)*B12D(I+1);
            TEMP = RWORK(IU2CS+I-1)*B22D(I) + RWORK(IU2SN+I-1)*B22E(I);
            B22E[I] = RWORK(IU2CS+I-1)*B22E(I) - RWORK(IU2SN+I-1)*B22D(I);
            B22D[I] = TEMP;
            B22BULGE = RWORK(IU2SN+I-1)*B22D(I+1);
            B22D[I+1] = RWORK(IU2CS+I-1)*B22D(I+1);

         }

         // Compute PHI(IMAX-1)

         X1 = SIN(THETA(IMAX-1))*B11E(IMAX-1) + COS(THETA(IMAX-1))*B21E(IMAX-1)          Y1 = SIN(THETA(IMAX-1))*B12D(IMAX-1) + COS(THETA(IMAX-1))*B22D(IMAX-1);
         Y2 = SIN(THETA(IMAX-1))*B12BULGE + COS(THETA(IMAX-1))*B22BULGE;

         PHI[IMAX-1] = ATAN2( (X1).abs(), sqrt(Y1**2+Y2**2) );

         // Chase bulges from B12(IMAX-1,IMAX) and B22(IMAX-1,IMAX)

         RESTART12 = B12D(IMAX-1)**2 + B12BULGE**2 <= THRESH**2;
         RESTART22 = B22D(IMAX-1)**2 + B22BULGE**2 <= THRESH**2;

         if ( !RESTART12 && !RESTART22 ) {
            slartgp(Y2, Y1, RWORK(IV2TSN+IMAX-1-1), RWORK(IV2TCS+IMAX-1-1), R );
         } else if ( !RESTART12 && RESTART22 ) {
            slartgp(B12BULGE, B12D(IMAX-1), RWORK(IV2TSN+IMAX-1-1), RWORK(IV2TCS+IMAX-1-1), R );
         } else if ( RESTART12 && !RESTART22 ) {
            slartgp(B22BULGE, B22D(IMAX-1), RWORK(IV2TSN+IMAX-1-1), RWORK(IV2TCS+IMAX-1-1), R );
         } else if ( NU < MU ) {
            slartgs(B12E(IMAX-1), B12D(IMAX), NU, RWORK(IV2TCS+IMAX-1-1), RWORK(IV2TSN+IMAX-1-1) );
         } else {
            slartgs(B22E(IMAX-1), B22D(IMAX), MU, RWORK(IV2TCS+IMAX-1-1), RWORK(IV2TSN+IMAX-1-1) );
         }

         TEMP = RWORK(IV2TCS+IMAX-1-1)*B12E(IMAX-1) + RWORK(IV2TSN+IMAX-1-1)*B12D(IMAX)          B12D(IMAX) = RWORK(IV2TCS+IMAX-1-1)*B12D(IMAX) - RWORK(IV2TSN+IMAX-1-1)*B12E(IMAX-1);
         B12E[IMAX-1] = TEMP;
         TEMP = RWORK(IV2TCS+IMAX-1-1)*B22E(IMAX-1) + RWORK(IV2TSN+IMAX-1-1)*B22D(IMAX)          B22D(IMAX) = RWORK(IV2TCS+IMAX-1-1)*B22D(IMAX) - RWORK(IV2TSN+IMAX-1-1)*B22E(IMAX-1);
         B22E[IMAX-1] = TEMP;

         // Update singular vectors

         if ( WANTU1 ) {
            if ( COLMAJOR ) {
               clasr('R', 'V', 'F', P, IMAX-IMIN+1, RWORK(IU1CS+IMIN-1), RWORK(IU1SN+IMIN-1), U1(1,IMIN), LDU1 );
            } else {
               clasr('L', 'V', 'F', IMAX-IMIN+1, P, RWORK(IU1CS+IMIN-1), RWORK(IU1SN+IMIN-1), U1(IMIN,1), LDU1 );
            }
         }
         if ( WANTU2 ) {
            if ( COLMAJOR ) {
               clasr('R', 'V', 'F', M-P, IMAX-IMIN+1, RWORK(IU2CS+IMIN-1), RWORK(IU2SN+IMIN-1), U2(1,IMIN), LDU2 );
            } else {
               clasr('L', 'V', 'F', IMAX-IMIN+1, M-P, RWORK(IU2CS+IMIN-1), RWORK(IU2SN+IMIN-1), U2(IMIN,1), LDU2 );
            }
         }
         if ( WANTV1T ) {
            if ( COLMAJOR ) {
               clasr('L', 'V', 'F', IMAX-IMIN+1, Q, RWORK(IV1TCS+IMIN-1), RWORK(IV1TSN+IMIN-1), V1T(IMIN,1), LDV1T );
            } else {
               clasr('R', 'V', 'F', Q, IMAX-IMIN+1, RWORK(IV1TCS+IMIN-1), RWORK(IV1TSN+IMIN-1), V1T(1,IMIN), LDV1T );
            }
         }
         if ( WANTV2T ) {
            if ( COLMAJOR ) {
               clasr('L', 'V', 'F', IMAX-IMIN+1, M-Q, RWORK(IV2TCS+IMIN-1), RWORK(IV2TSN+IMIN-1), V2T(IMIN,1), LDV2T );
            } else {
               clasr('R', 'V', 'F', M-Q, IMAX-IMIN+1, RWORK(IV2TCS+IMIN-1), RWORK(IV2TSN+IMIN-1), V2T(1,IMIN), LDV2T );
            }
         }

         // Fix signs on B11(IMAX-1,IMAX) and B21(IMAX-1,IMAX)

         if ( B11E(IMAX-1)+B21E(IMAX-1) > 0 ) {
            B11D[IMAX] = -B11D(IMAX);
            B21D[IMAX] = -B21D(IMAX);
            if ( WANTV1T ) {
               if ( COLMAJOR ) {
                  cscal(Q, NEGONECOMPLEX, V1T(IMAX,1), LDV1T );
               } else {
                  cscal(Q, NEGONECOMPLEX, V1T(1,IMAX), 1 );
               }
            }
         }

         // Compute THETA(IMAX)

         X1 = COS(PHI(IMAX-1))*B11D(IMAX) + SIN(PHI(IMAX-1))*B12E(IMAX-1)          Y1 = COS(PHI(IMAX-1))*B21D(IMAX) + SIN(PHI(IMAX-1))*B22E(IMAX-1);

         THETA[IMAX] = ATAN2( (Y1).abs(), (X1).abs() );

         // Fix signs on B11(IMAX,IMAX), B12(IMAX,IMAX-1), B21(IMAX,IMAX),
         // and B22(IMAX,IMAX-1)

         if ( B11D(IMAX)+B12E(IMAX-1) < 0 ) {
            B12D[IMAX] = -B12D(IMAX);
            if ( WANTU1 ) {
               if ( COLMAJOR ) {
                  cscal(P, NEGONECOMPLEX, U1(1,IMAX), 1 );
               } else {
                  cscal(P, NEGONECOMPLEX, U1(IMAX,1), LDU1 );
               }
            }
         }
         if ( B21D(IMAX)+B22E(IMAX-1) > 0 ) {
            B22D[IMAX] = -B22D(IMAX);
            if ( WANTU2 ) {
               if ( COLMAJOR ) {
                  cscal(M-P, NEGONECOMPLEX, U2(1,IMAX), 1 );
               } else {
                  cscal(M-P, NEGONECOMPLEX, U2(IMAX,1), LDU2 );
               }
            }
         }

         // Fix signs on B12(IMAX,IMAX) and B22(IMAX,IMAX)

         if ( B12D(IMAX)+B22D(IMAX) < 0 ) {
            if ( WANTV2T ) {
               if ( COLMAJOR ) {
                  cscal(M-Q, NEGONECOMPLEX, V2T(IMAX,1), LDV2T );
               } else {
                  cscal(M-Q, NEGONECOMPLEX, V2T(1,IMAX), 1 );
               }
            }
         }

         // Test for negligible sines or cosines

         for (I = IMIN; I <= IMAX; I++) {
            if ( THETA(I) < THRESH ) {
               THETA[I] = ZERO;
            } else if ( THETA(I) > PIOVER2-THRESH ) {
               THETA[I] = PIOVER2;
            }
         }
         for (I = IMIN; I <= IMAX-1; I++) {
            if ( PHI(I) < THRESH ) {
               PHI[I] = ZERO;
            } else if ( PHI(I) > PIOVER2-THRESH ) {
               PHI[I] = PIOVER2;
            }
         }

         // Deflate

         if (IMAX > 1) {
            while (PHI(IMAX-1) == ZERO) {
               IMAX = IMAX - 1;
               if (IMAX <= 1) break;
            }
         }
         if (IMIN > IMAX - 1) IMIN = IMAX - 1;
         if (IMIN > 1) {
            while (PHI(IMIN-1) != ZERO) {
                IMIN = IMIN - 1;
                if (IMIN <= 1) break;
            }
         }

         // Repeat main iteration loop

      }

      // Postprocessing: order THETA from least to greatest

      for (I = 1; I <= Q; I++) {

         MINI = I;
         THETAMIN = THETA(I);
         for (J = I+1; J <= Q; J++) {
            if ( THETA(J) < THETAMIN ) {
               MINI = J;
               THETAMIN = THETA(J);
            }
         }

         if ( MINI != I ) {
            THETA[MINI] = THETA(I);
            THETA[I] = THETAMIN;
            if ( COLMAJOR ) {
               if (WANTU1) cswap( P, U1(1,I), 1, U1(1,MINI), 1 );
               if( WANTU2 ) cswap( M-P, U2(1,I), 1, U2(1,MINI), 1 );
               if( WANTV1T ) cswap( Q, V1T(I,1), LDV1T, V1T(MINI,1), LDV1T );
               IF( WANTV2T ) cswap( M-Q, V2T(I,1), LDV2T, V2T(MINI,1), LDV2T );
            } else {
               if (WANTU1) cswap( P, U1(I,1), LDU1, U1(MINI,1), LDU1 );
               if( WANTU2 ) cswap( M-P, U2(I,1), LDU2, U2(MINI,1), LDU2 );
               if( WANTV1T ) cswap( Q, V1T(1,I), 1, V1T(1,MINI), 1 );
               IF( WANTV2T ) cswap( M-Q, V2T(1,I), 1, V2T(1,MINI), 1 );
            }
         }

      }

      }
