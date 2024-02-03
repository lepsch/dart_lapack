      SUBROUTINE DORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIGNS, TRANS;
      int                INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      double             PHI( * ), THETA( * );
      double             TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), WORK( * ), X11( LDX11, * ), X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, * );
      // ..

*  ====================================================================

      // .. Parameters ..
      double             REALONE;
      const              REALONE = 1.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               COLMAJOR, LQUERY;
      int                I, LWORKMIN, LWORKOPT;
      double             Z1, Z2, Z3, Z4;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLARF, DLARFGP, DSCAL, XERBLA
      // ..
      // .. External Functions ..
      double             DNRM2;
      bool               LSAME;
      // EXTERNAL DNRM2, LSAME
      // ..
      // .. Intrinsic Functions
      // INTRINSIC ATAN2, COS, MAX, SIN
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0;
      COLMAJOR = !LSAME( TRANS, 'T' );
      if ( !LSAME( SIGNS, 'O' ) ) {
         Z1 = REALONE;
         Z2 = REALONE;
         Z3 = REALONE;
         Z4 = REALONE;
      } else {
         Z1 = REALONE;
         Z2 = -REALONE;
         Z3 = REALONE;
         Z4 = -REALONE;
      }
      LQUERY = LWORK == -1;

      if ( M < 0 ) {
         INFO = -3;
      } else if ( P < 0 || P > M ) {
         INFO = -4;
      } else if ( Q < 0 || Q > P || Q > M-P || Q > M-Q ) {
         INFO = -5;
      } else if ( COLMAJOR && LDX11 < MAX( 1, P ) ) {
         INFO = -7;
      } else if ( !COLMAJOR && LDX11 < MAX( 1, Q ) ) {
         INFO = -7;
      } else if ( COLMAJOR && LDX12 < MAX( 1, P ) ) {
         INFO = -9;
      } else if ( !COLMAJOR && LDX12 < MAX( 1, M-Q ) ) {
         INFO = -9;
      } else if ( COLMAJOR && LDX21 < MAX( 1, M-P ) ) {
         INFO = -11;
      } else if ( !COLMAJOR && LDX21 < MAX( 1, Q ) ) {
         INFO = -11;
      } else if ( COLMAJOR && LDX22 < MAX( 1, M-P ) ) {
         INFO = -13;
      } else if ( !COLMAJOR && LDX22 < MAX( 1, M-Q ) ) {
         INFO = -13;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         LWORKOPT = M - Q;
         LWORKMIN = M - Q;
         WORK(1) = LWORKOPT;
         if ( LWORK < LWORKMIN && !LQUERY ) {
            INFO = -21;
         }
      }
      if ( INFO != 0 ) {
         xerbla('xORBDB', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Handle column-major and row-major separately

      if ( COLMAJOR ) {

         // Reduce columns 1, ..., Q of X11, X12, X21, and X22

         for (I = 1; I <= Q; I++) {

            if ( I == 1 ) {
               dscal(P-I+1, Z1, X11(I,I), 1 );
            } else {
               dscal(P-I+1, Z1*COS(PHI(I-1)), X11(I,I), 1 );
               daxpy(P-I+1, -Z1*Z3*Z4*SIN(PHI(I-1)), X12(I,I-1), 1, X11(I,I), 1 );
            }
            if ( I == 1 ) {
               dscal(M-P-I+1, Z2, X21(I,I), 1 );
            } else {
               dscal(M-P-I+1, Z2*COS(PHI(I-1)), X21(I,I), 1 );
               daxpy(M-P-I+1, -Z2*Z3*Z4*SIN(PHI(I-1)), X22(I,I-1), 1, X21(I,I), 1 );
            }

            THETA(I) = ATAN2( DNRM2( M-P-I+1, X21(I,I), 1 ), DNRM2( P-I+1, X11(I,I), 1 ) );

            if ( P > I ) {
               dlarfgp(P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) );
            } else if ( P == I ) {
               dlarfgp(P-I+1, X11(I,I), X11(I,I), 1, TAUP1(I) );
            }
            X11(I,I) = ONE;
            if ( M-P > I ) {
               dlarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
            } else if ( M-P == I ) {
               dlarfgp(M-P-I+1, X21(I,I), X21(I,I), 1, TAUP2(I) );
            }
            X21(I,I) = ONE;

            if ( Q > I ) {
               dlarf('L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK );
            }
            if ( M-Q+1 > I ) {
               dlarf('L', P-I+1, M-Q-I+1, X11(I,I), 1, TAUP1(I), X12(I,I), LDX12, WORK );
            }
            if ( Q > I ) {
               dlarf('L', M-P-I+1, Q-I, X21(I,I), 1, TAUP2(I), X21(I,I+1), LDX21, WORK );
            }
            if ( M-Q+1 > I ) {
               dlarf('L', M-P-I+1, M-Q-I+1, X21(I,I), 1, TAUP2(I), X22(I,I), LDX22, WORK );
            }

            if ( I < Q ) {
               dscal(Q-I, -Z1*Z3*SIN(THETA(I)), X11(I,I+1), LDX11 );
               daxpy(Q-I, Z2*Z3*COS(THETA(I)), X21(I,I+1), LDX21, X11(I,I+1), LDX11 );
            }
            dscal(M-Q-I+1, -Z1*Z4*SIN(THETA(I)), X12(I,I), LDX12 );
            daxpy(M-Q-I+1, Z2*Z4*COS(THETA(I)), X22(I,I), LDX22, X12(I,I), LDX12 );

            if (I < Q) PHI(I) = ATAN2( DNRM2( Q-I, X11(I,I+1), LDX11 ), DNRM2( M-Q-I+1, X12(I,I), LDX12 ) );

            if ( I < Q ) {
               if ( Q-I == 1 ) {
                  dlarfgp(Q-I, X11(I,I+1), X11(I,I+1), LDX11, TAUQ1(I) );
               } else {
                  dlarfgp(Q-I, X11(I,I+1), X11(I,I+2), LDX11, TAUQ1(I) );
               }
               X11(I,I+1) = ONE;
            }
            if ( Q+I-1 < M ) {
               if ( M-Q == I ) {
                  dlarfgp(M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) );
               } else {
                  dlarfgp(M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) );
               }
            }
            X12(I,I) = ONE;

            if ( I < Q ) {
               dlarf('R', P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X11(I+1,I+1), LDX11, WORK );
               dlarf('R', M-P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X21(I+1,I+1), LDX21, WORK );
            }
            if ( P > I ) {
               dlarf('R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK );
            }
            if ( M-P > I ) {
               dlarf('R', M-P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(I+1,I), LDX22, WORK );
            }

         }

         // Reduce columns Q + 1, ..., P of X12, X22

         for (I = Q + 1; I <= P; I++) {

            dscal(M-Q-I+1, -Z1*Z4, X12(I,I), LDX12 );
            if ( I >= M-Q ) {
               dlarfgp(M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) );
            } else {
               dlarfgp(M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) );
            }
            X12(I,I) = ONE;

            if ( P > I ) {
               dlarf('R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK );
            }
            if (M-P-Q >= 1) CALL DLARF( 'R', M-P-Q, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(Q+1,I), LDX22, WORK );

         }

         // Reduce columns P + 1, ..., M - Q of X12, X22

         for (I = 1; I <= M - P - Q; I++) {

            dscal(M-P-Q-I+1, Z2*Z4, X22(Q+I,P+I), LDX22 );
            if ( I == M-P-Q ) {
               dlarfgp(M-P-Q-I+1, X22(Q+I,P+I), X22(Q+I,P+I), LDX22, TAUQ2(P+I) );
            } else {
               dlarfgp(M-P-Q-I+1, X22(Q+I,P+I), X22(Q+I,P+I+1), LDX22, TAUQ2(P+I) );
            }
            X22(Q+I,P+I) = ONE;
            if ( I < M-P-Q ) {
               dlarf('R', M-P-Q-I, M-P-Q-I+1, X22(Q+I,P+I), LDX22, TAUQ2(P+I), X22(Q+I+1,P+I), LDX22, WORK );
            }

         }

      } else {

         // Reduce columns 1, ..., Q of X11, X12, X21, X22

         for (I = 1; I <= Q; I++) {

            if ( I == 1 ) {
               dscal(P-I+1, Z1, X11(I,I), LDX11 );
            } else {
               dscal(P-I+1, Z1*COS(PHI(I-1)), X11(I,I), LDX11 );
               daxpy(P-I+1, -Z1*Z3*Z4*SIN(PHI(I-1)), X12(I-1,I), LDX12, X11(I,I), LDX11 );
            }
            if ( I == 1 ) {
               dscal(M-P-I+1, Z2, X21(I,I), LDX21 );
            } else {
               dscal(M-P-I+1, Z2*COS(PHI(I-1)), X21(I,I), LDX21 );
               daxpy(M-P-I+1, -Z2*Z3*Z4*SIN(PHI(I-1)), X22(I-1,I), LDX22, X21(I,I), LDX21 );
            }

            THETA(I) = ATAN2( DNRM2( M-P-I+1, X21(I,I), LDX21 ), DNRM2( P-I+1, X11(I,I), LDX11 ) );

            dlarfgp(P-I+1, X11(I,I), X11(I,I+1), LDX11, TAUP1(I) );
            X11(I,I) = ONE;
            if ( I == M-P ) {
               dlarfgp(M-P-I+1, X21(I,I), X21(I,I), LDX21, TAUP2(I) );
            } else {
               dlarfgp(M-P-I+1, X21(I,I), X21(I,I+1), LDX21, TAUP2(I) );
            }
            X21(I,I) = ONE;

            if ( Q > I ) {
               dlarf('R', Q-I, P-I+1, X11(I,I), LDX11, TAUP1(I), X11(I+1,I), LDX11, WORK );
            }
            if ( M-Q+1 > I ) {
               dlarf('R', M-Q-I+1, P-I+1, X11(I,I), LDX11, TAUP1(I), X12(I,I), LDX12, WORK );
            }
            if ( Q > I ) {
               dlarf('R', Q-I, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X21(I+1,I), LDX21, WORK );
            }
            if ( M-Q+1 > I ) {
               dlarf('R', M-Q-I+1, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X22(I,I), LDX22, WORK );
            }

            if ( I < Q ) {
               dscal(Q-I, -Z1*Z3*SIN(THETA(I)), X11(I+1,I), 1 );
               daxpy(Q-I, Z2*Z3*COS(THETA(I)), X21(I+1,I), 1, X11(I+1,I), 1 );
            }
            dscal(M-Q-I+1, -Z1*Z4*SIN(THETA(I)), X12(I,I), 1 );
            daxpy(M-Q-I+1, Z2*Z4*COS(THETA(I)), X22(I,I), 1, X12(I,I), 1 );

            if (I < Q) PHI(I) = ATAN2( DNRM2( Q-I, X11(I+1,I), 1 ), DNRM2( M-Q-I+1, X12(I,I), 1 ) );

            if ( I < Q ) {
               if ( Q-I == 1) {
                  dlarfgp(Q-I, X11(I+1,I), X11(I+1,I), 1, TAUQ1(I) );
               } else {
                  dlarfgp(Q-I, X11(I+1,I), X11(I+2,I), 1, TAUQ1(I) );
               }
               X11(I+1,I) = ONE;
            }
            if ( M-Q > I ) {
               dlarfgp(M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) );
            } else {
               dlarfgp(M-Q-I+1, X12(I,I), X12(I,I), 1, TAUQ2(I) );
            }
            X12(I,I) = ONE;

            if ( I < Q ) {
               dlarf('L', Q-I, P-I, X11(I+1,I), 1, TAUQ1(I), X11(I+1,I+1), LDX11, WORK );
               dlarf('L', Q-I, M-P-I, X11(I+1,I), 1, TAUQ1(I), X21(I+1,I+1), LDX21, WORK );
            }
            dlarf('L', M-Q-I+1, P-I, X12(I,I), 1, TAUQ2(I), X12(I,I+1), LDX12, WORK );
            if ( M-P-I > 0 ) {
               dlarf('L', M-Q-I+1, M-P-I, X12(I,I), 1, TAUQ2(I), X22(I,I+1), LDX22, WORK );
            }

         }

         // Reduce columns Q + 1, ..., P of X12, X22

         for (I = Q + 1; I <= P; I++) {

            dscal(M-Q-I+1, -Z1*Z4, X12(I,I), 1 );
            dlarfgp(M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) );
            X12(I,I) = ONE;

            if ( P > I ) {
               dlarf('L', M-Q-I+1, P-I, X12(I,I), 1, TAUQ2(I), X12(I,I+1), LDX12, WORK );
            }
            if (M-P-Q >= 1) CALL DLARF( 'L', M-Q-I+1, M-P-Q, X12(I,I), 1, TAUQ2(I), X22(I,Q+1), LDX22, WORK );

         }

         // Reduce columns P + 1, ..., M - Q of X12, X22

         for (I = 1; I <= M - P - Q; I++) {

            dscal(M-P-Q-I+1, Z2*Z4, X22(P+I,Q+I), 1 );
            if ( M-P-Q == I ) {
               dlarfgp(M-P-Q-I+1, X22(P+I,Q+I), X22(P+I,Q+I), 1, TAUQ2(P+I) );
            } else {
               dlarfgp(M-P-Q-I+1, X22(P+I,Q+I), X22(P+I+1,Q+I), 1, TAUQ2(P+I) );
               dlarf('L', M-P-Q-I+1, M-P-Q-I, X22(P+I,Q+I), 1, TAUQ2(P+I), X22(P+I,Q+I+1), LDX22, WORK );
            }
            X22(P+I,Q+I) = ONE;

         }

      }

      return;

      // End of DORBDB

      }
