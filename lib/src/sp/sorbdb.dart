      void sorbdb(TRANS, SIGNS, M, P, Q, final Matrix<double> X11, final int LDX11, final Matrix<double> X12, final int LDX12, final Matrix<double> X21, final int LDX21, final Matrix<double> X22, final int LDX22, THETA, PHI, TAUP1, TAUP2, TAUQ1, TAUQ2, final Array<double> WORK, final int LWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIGNS, TRANS;
      int                INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, Q;
      double               PHI( * ), THETA( * );
      double               TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), WORK( * ), X11( LDX11, * ), X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, * );
      // ..

// ====================================================================

      // .. Parameters ..
      double               REALONE;
      const              REALONE = 1.0 ;
      double               ONE;
      const              ONE = 1.0 ;
      bool               COLMAJOR, LQUERY;
      int                I, LWORKMIN, LWORKOPT;
      double               Z1, Z2, Z3, Z4;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SLARF, SLARFGP, SSCAL, XERBLA
      // ..
      // .. External Functions ..
      //- REAL               SNRM2;
      //- bool               lsame;
      // EXTERNAL SNRM2, lsame
      // ..
      // .. Intrinsic Functions
      // INTRINSIC ATAN2, COS, MAX, SIN

      // Test input arguments

      INFO = 0;
      COLMAJOR = !lsame( TRANS, 'T' );
      if ( !lsame( SIGNS, 'O' ) ) {
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
      } else if ( COLMAJOR && LDX11 < max( 1, P ) ) {
         INFO = -7;
      } else if ( !COLMAJOR && LDX11 < max( 1, Q ) ) {
         INFO = -7;
      } else if ( COLMAJOR && LDX12 < max( 1, P ) ) {
         INFO = -9;
      } else if ( !COLMAJOR && LDX12 < max( 1, M-Q ) ) {
         INFO = -9;
      } else if ( COLMAJOR && LDX21 < max( 1, M-P ) ) {
         INFO = -11;
      } else if ( !COLMAJOR && LDX21 < max( 1, Q ) ) {
         INFO = -11;
      } else if ( COLMAJOR && LDX22 < max( 1, M-P ) ) {
         INFO = -13;
      } else if ( !COLMAJOR && LDX22 < max( 1, M-Q ) ) {
         INFO = -13;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         LWORKOPT = M - Q;
         LWORKMIN = M - Q;
         WORK[1] = LWORKOPT;
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
               sscal(P-I+1, Z1, X11(I,I), 1 );
            } else {
               sscal(P-I+1, Z1*COS(PHI(I-1)), X11(I,I), 1 );
               saxpy(P-I+1, -Z1*Z3*Z4*SIN(PHI(I-1)), X12(I,I-1), 1, X11(I,I), 1 );
            }
            if ( I == 1 ) {
               sscal(M-P-I+1, Z2, X21(I,I), 1 );
            } else {
               sscal(M-P-I+1, Z2*COS(PHI(I-1)), X21(I,I), 1 );
               saxpy(M-P-I+1, -Z2*Z3*Z4*SIN(PHI(I-1)), X22(I,I-1), 1, X21(I,I), 1 );
            }

            THETA[I] = ATAN2( SNRM2( M-P-I+1, X21(I,I), 1 ), SNRM2( P-I+1, X11(I,I), 1 ) );

            if ( P > I ) {
               slarfgp(P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) );
            } else if ( P == I ) {
               slarfgp(P-I+1, X11(I,I), X11(I,I), 1, TAUP1(I) );
            }
            X11[I][I] = ONE;
            if ( M-P > I ) {
               slarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
            } else if ( M-P == I ) {
               slarfgp(M-P-I+1, X21(I,I), X21(I,I), 1, TAUP2(I) );
            }
            X21[I][I] = ONE;

            if ( Q > I ) {
               slarf('L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK );
            }
            if ( M-Q+1 > I ) {
               slarf('L', P-I+1, M-Q-I+1, X11(I,I), 1, TAUP1(I), X12(I,I), LDX12, WORK );
            }
            if ( Q > I ) {
               slarf('L', M-P-I+1, Q-I, X21(I,I), 1, TAUP2(I), X21(I,I+1), LDX21, WORK );
            }
            if ( M-Q+1 > I ) {
               slarf('L', M-P-I+1, M-Q-I+1, X21(I,I), 1, TAUP2(I), X22(I,I), LDX22, WORK );
            }

            if ( I < Q ) {
               sscal(Q-I, -Z1*Z3*SIN(THETA(I)), X11(I,I+1), LDX11 );
               saxpy(Q-I, Z2*Z3*COS(THETA(I)), X21(I,I+1), LDX21, X11(I,I+1), LDX11 );
            }
            sscal(M-Q-I+1, -Z1*Z4*SIN(THETA(I)), X12(I,I), LDX12 );
            saxpy(M-Q-I+1, Z2*Z4*COS(THETA(I)), X22(I,I), LDX22, X12(I,I), LDX12 );

            if (I < Q) PHI(I) = ATAN2( SNRM2( Q-I, X11(I,I+1), LDX11 ), SNRM2( M-Q-I+1, X12(I,I), LDX12 ) );

            if ( I < Q ) {
               if ( Q-I == 1 ) {
                  slarfgp(Q-I, X11(I,I+1), X11(I,I+1), LDX11, TAUQ1(I) );
               } else {
                  slarfgp(Q-I, X11(I,I+1), X11(I,I+2), LDX11, TAUQ1(I) );
               }
               X11[I][I+1] = ONE;
            }
            if ( Q+I-1 < M ) {
               if ( M-Q == I ) {
                  slarfgp(M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) );
               } else {
                  slarfgp(M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) );
               }
            }
            X12[I][I] = ONE;

            if ( I < Q ) {
               slarf('R', P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X11(I+1,I+1), LDX11, WORK );
               slarf('R', M-P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X21(I+1,I+1), LDX21, WORK );
            }
            if ( P > I ) {
               slarf('R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK );
            }
            if ( M-P > I ) {
               slarf('R', M-P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(I+1,I), LDX22, WORK );
            }

         }

         // Reduce columns Q + 1, ..., P of X12, X22

         for (I = Q + 1; I <= P; I++) {

            sscal(M-Q-I+1, -Z1*Z4, X12(I,I), LDX12 );
            if ( I >= M-Q ) {
               slarfgp(M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) );
            } else {
               slarfgp(M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) );
            }
            X12[I][I] = ONE;

            if ( P > I ) {
               slarf('R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK );
            }
            if (M-P-Q >= 1) slarf( 'R', M-P-Q, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(Q+1,I), LDX22, WORK );

         }

         // Reduce columns P + 1, ..., M - Q of X12, X22

         for (I = 1; I <= M - P - Q; I++) {

            sscal(M-P-Q-I+1, Z2*Z4, X22(Q+I,P+I), LDX22 );
            if ( I == M-P-Q ) {
               slarfgp(M-P-Q-I+1, X22(Q+I,P+I), X22(Q+I,P+I), LDX22, TAUQ2(P+I) );
            } else {
               slarfgp(M-P-Q-I+1, X22(Q+I,P+I), X22(Q+I,P+I+1), LDX22, TAUQ2(P+I) );
            }
            X22[Q+I][P+I] = ONE;
            if ( I < M-P-Q ) {
               slarf('R', M-P-Q-I, M-P-Q-I+1, X22(Q+I,P+I), LDX22, TAUQ2(P+I), X22(Q+I+1,P+I), LDX22, WORK );
            }

         }

      } else {

         // Reduce columns 1, ..., Q of X11, X12, X21, X22

         for (I = 1; I <= Q; I++) {

            if ( I == 1 ) {
               sscal(P-I+1, Z1, X11(I,I), LDX11 );
            } else {
               sscal(P-I+1, Z1*COS(PHI(I-1)), X11(I,I), LDX11 );
               saxpy(P-I+1, -Z1*Z3*Z4*SIN(PHI(I-1)), X12(I-1,I), LDX12, X11(I,I), LDX11 );
            }
            if ( I == 1 ) {
               sscal(M-P-I+1, Z2, X21(I,I), LDX21 );
            } else {
               sscal(M-P-I+1, Z2*COS(PHI(I-1)), X21(I,I), LDX21 );
               saxpy(M-P-I+1, -Z2*Z3*Z4*SIN(PHI(I-1)), X22(I-1,I), LDX22, X21(I,I), LDX21 );
            }

            THETA[I] = ATAN2( SNRM2( M-P-I+1, X21(I,I), LDX21 ), SNRM2( P-I+1, X11(I,I), LDX11 ) );

            slarfgp(P-I+1, X11(I,I), X11(I,I+1), LDX11, TAUP1(I) );
            X11[I][I] = ONE;
            if ( I == M-P ) {
               slarfgp(M-P-I+1, X21(I,I), X21(I,I), LDX21, TAUP2(I) );
            } else {
               slarfgp(M-P-I+1, X21(I,I), X21(I,I+1), LDX21, TAUP2(I) );
            }
            X21[I][I] = ONE;

            if ( Q > I ) {
               slarf('R', Q-I, P-I+1, X11(I,I), LDX11, TAUP1(I), X11(I+1,I), LDX11, WORK );
            }
            if ( M-Q+1 > I ) {
               slarf('R', M-Q-I+1, P-I+1, X11(I,I), LDX11, TAUP1(I), X12(I,I), LDX12, WORK );
            }
            if ( Q > I ) {
               slarf('R', Q-I, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X21(I+1,I), LDX21, WORK );
            }
            if ( M-Q+1 > I ) {
               slarf('R', M-Q-I+1, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X22(I,I), LDX22, WORK );
            }

            if ( I < Q ) {
               sscal(Q-I, -Z1*Z3*SIN(THETA(I)), X11(I+1,I), 1 );
               saxpy(Q-I, Z2*Z3*COS(THETA(I)), X21(I+1,I), 1, X11(I+1,I), 1 );
            }
            sscal(M-Q-I+1, -Z1*Z4*SIN(THETA(I)), X12(I,I), 1 );
            saxpy(M-Q-I+1, Z2*Z4*COS(THETA(I)), X22(I,I), 1, X12(I,I), 1 );

            if (I < Q) PHI(I) = ATAN2( SNRM2( Q-I, X11(I+1,I), 1 ), SNRM2( M-Q-I+1, X12(I,I), 1 ) );

            if ( I < Q ) {
               if ( Q-I == 1) {
                  slarfgp(Q-I, X11(I+1,I), X11(I+1,I), 1, TAUQ1(I) );
               } else {
                  slarfgp(Q-I, X11(I+1,I), X11(I+2,I), 1, TAUQ1(I) );
               }
               X11[I+1][I] = ONE;
            }
            if ( M-Q > I ) {
               slarfgp(M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) );
            } else {
               slarfgp(M-Q-I+1, X12(I,I), X12(I,I), 1, TAUQ2(I) );
            }
            X12[I][I] = ONE;

            if ( I < Q ) {
               slarf('L', Q-I, P-I, X11(I+1,I), 1, TAUQ1(I), X11(I+1,I+1), LDX11, WORK );
               slarf('L', Q-I, M-P-I, X11(I+1,I), 1, TAUQ1(I), X21(I+1,I+1), LDX21, WORK );
            }
            slarf('L', M-Q-I+1, P-I, X12(I,I), 1, TAUQ2(I), X12(I,I+1), LDX12, WORK );
            if ( M-P-I > 0 ) {
               slarf('L', M-Q-I+1, M-P-I, X12(I,I), 1, TAUQ2(I), X22(I,I+1), LDX22, WORK );
            }

         }

         // Reduce columns Q + 1, ..., P of X12, X22

         for (I = Q + 1; I <= P; I++) {

            sscal(M-Q-I+1, -Z1*Z4, X12(I,I), 1 );
            slarfgp(M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) );
            X12[I][I] = ONE;

            if ( P > I ) {
               slarf('L', M-Q-I+1, P-I, X12(I,I), 1, TAUQ2(I), X12(I,I+1), LDX12, WORK );
            }
            if (M-P-Q >= 1) slarf( 'L', M-Q-I+1, M-P-Q, X12(I,I), 1, TAUQ2(I), X22(I,Q+1), LDX22, WORK );

         }

         // Reduce columns P + 1, ..., M - Q of X12, X22

         for (I = 1; I <= M - P - Q; I++) {

            sscal(M-P-Q-I+1, Z2*Z4, X22(P+I,Q+I), 1 );
            if ( M-P-Q == I ) {
               slarfgp(M-P-Q-I+1, X22(P+I,Q+I), X22(P+I,Q+I), 1, TAUQ2(P+I) );
               X22[P+I][Q+I] = ONE;
            } else {
               slarfgp(M-P-Q-I+1, X22(P+I,Q+I), X22(P+I+1,Q+I), 1, TAUQ2(P+I) );
               X22[P+I][Q+I] = ONE;
               slarf('L', M-P-Q-I+1, M-P-Q-I, X22(P+I,Q+I), 1, TAUQ2(P+I), X22(P+I,Q+I+1), LDX22, WORK );
            }


         }

      }

      }
