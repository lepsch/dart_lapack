      SUBROUTINE CUNBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIGNS, TRANS;
      int                INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      REAL               PHI( * ), THETA( * )
      COMPLEX            TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), WORK( * ), X11( LDX11, * ), X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, * )
      // ..

*  ====================================================================

      // .. Parameters ..
      REAL               REALONE
      const              REALONE = 1.0E0 ;
      COMPLEX            ONE
      const              ONE = (1.0E0,0.0E0) ;
      // ..
      // .. Local Scalars ..
      bool               COLMAJOR, LQUERY;
      int                I, LWORKMIN, LWORKOPT;
      REAL               Z1, Z2, Z3, Z4
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CLARF, CLARFGP, CSCAL, XERBLA
      // EXTERNAL CLACGV

      // ..
      // .. External Functions ..
      REAL               SCNRM2, SROUNDUP_LWORK
      bool               LSAME;
      // EXTERNAL SCNRM2, SROUNDUP_LWORK, LSAME
      // ..
      // .. Intrinsic Functions
      // INTRINSIC ATAN2, COS, MAX, MIN, SIN
      // INTRINSIC CMPLX, CONJG
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      COLMAJOR = .NOT. LSAME( TRANS, 'T' )
      if ( .NOT. LSAME( SIGNS, 'O' ) ) {
         Z1 = REALONE
         Z2 = REALONE
         Z3 = REALONE
         Z4 = REALONE
      } else {
         Z1 = REALONE
         Z2 = -REALONE
         Z3 = REALONE
         Z4 = -REALONE
      }
      LQUERY = LWORK == -1

      if ( M .LT. 0 ) {
         INFO = -3
      } else if ( P .LT. 0 .OR. P .GT. M ) {
         INFO = -4
      } else if ( Q .LT. 0 .OR. Q .GT. P .OR. Q .GT. M-P .OR. Q .GT. M-Q ) {
         INFO = -5
      } else if ( COLMAJOR && LDX11 .LT. MAX( 1, P ) ) {
         INFO = -7
      } else if ( .NOT.COLMAJOR && LDX11 .LT. MAX( 1, Q ) ) {
         INFO = -7
      } else if ( COLMAJOR && LDX12 .LT. MAX( 1, P ) ) {
         INFO = -9
      } else if ( .NOT.COLMAJOR && LDX12 .LT. MAX( 1, M-Q ) ) {
         INFO = -9
      } else if ( COLMAJOR && LDX21 .LT. MAX( 1, M-P ) ) {
         INFO = -11
      } else if ( .NOT.COLMAJOR && LDX21 .LT. MAX( 1, Q ) ) {
         INFO = -11
      } else if ( COLMAJOR && LDX22 .LT. MAX( 1, M-P ) ) {
         INFO = -13
      } else if ( .NOT.COLMAJOR && LDX22 .LT. MAX( 1, M-Q ) ) {
         INFO = -13
      }

      // Compute workspace

      if ( INFO == 0 ) {
         LWORKOPT = M - Q
         LWORKMIN = M - Q
         WORK(1) = SROUNDUP_LWORK(LWORKOPT)
         if ( LWORK .LT. LWORKMIN && .NOT. LQUERY ) {
            INFO = -21
         }
      }
      if ( INFO != 0 ) {
         xerbla('xORBDB', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Handle column-major and row-major separately

      if ( COLMAJOR ) {

         // Reduce columns 1, ..., Q of X11, X12, X21, and X22

         for (I = 1; I <= Q; I++) {

            if ( I == 1 ) {
               cscal(P-I+1, CMPLX( Z1, 0.0E0 ), X11(I,I), 1 );
            } else {
               cscal(P-I+1, CMPLX( Z1*COS(PHI(I-1)), 0.0E0 ), X11(I,I), 1 );
               caxpy(P-I+1, CMPLX( -Z1*Z3*Z4*SIN(PHI(I-1)), 0.0E0 ), X12(I,I-1), 1, X11(I,I), 1 );
            }
            if ( I == 1 ) {
               cscal(M-P-I+1, CMPLX( Z2, 0.0E0 ), X21(I,I), 1 );
            } else {
               cscal(M-P-I+1, CMPLX( Z2*COS(PHI(I-1)), 0.0E0 ), X21(I,I), 1 );
               caxpy(M-P-I+1, CMPLX( -Z2*Z3*Z4*SIN(PHI(I-1)), 0.0E0 ), X22(I,I-1), 1, X21(I,I), 1 );
            }

            THETA(I) = ATAN2( SCNRM2( M-P-I+1, X21(I,I), 1 ), SCNRM2( P-I+1, X11(I,I), 1 ) )

            if ( P .GT. I ) {
               clarfgp(P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) );
            } else if ( P == I ) {
               clarfgp(P-I+1, X11(I,I), X11(I,I), 1, TAUP1(I) );
            }
            X11(I,I) = ONE
            if ( M-P .GT. I ) {
               clarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
            } else if ( M-P == I ) {
               clarfgp(M-P-I+1, X21(I,I), X21(I,I), 1, TAUP2(I) );
            }
            X21(I,I) = ONE

            if ( Q .GT. I ) {
               clarf('L', P-I+1, Q-I, X11(I,I), 1, CONJG(TAUP1(I)), X11(I,I+1), LDX11, WORK );
               clarf('L', M-P-I+1, Q-I, X21(I,I), 1, CONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK );
            }
            if ( M-Q+1 .GT. I ) {
               clarf('L', P-I+1, M-Q-I+1, X11(I,I), 1, CONJG(TAUP1(I)), X12(I,I), LDX12, WORK );
               clarf('L', M-P-I+1, M-Q-I+1, X21(I,I), 1, CONJG(TAUP2(I)), X22(I,I), LDX22, WORK );
            }

            if ( I .LT. Q ) {
               cscal(Q-I, CMPLX( -Z1*Z3*SIN(THETA(I)), 0.0E0 ), X11(I,I+1), LDX11 );
               caxpy(Q-I, CMPLX( Z2*Z3*COS(THETA(I)), 0.0E0 ), X21(I,I+1), LDX21, X11(I,I+1), LDX11 );
            }
            cscal(M-Q-I+1, CMPLX( -Z1*Z4*SIN(THETA(I)), 0.0E0 ), X12(I,I), LDX12 );
            caxpy(M-Q-I+1, CMPLX( Z2*Z4*COS(THETA(I)), 0.0E0 ), X22(I,I), LDX22, X12(I,I), LDX12 );

            if (I .LT. Q) PHI(I) = ATAN2( SCNRM2( Q-I, X11(I,I+1), LDX11 ), SCNRM2( M-Q-I+1, X12(I,I), LDX12 ) );

            if ( I .LT. Q ) {
               clacgv(Q-I, X11(I,I+1), LDX11 );
               if ( I == Q-1 ) {
                  clarfgp(Q-I, X11(I,I+1), X11(I,I+1), LDX11, TAUQ1(I) );
               } else {
                  clarfgp(Q-I, X11(I,I+1), X11(I,I+2), LDX11, TAUQ1(I) );
               }
               X11(I,I+1) = ONE
            }
            if ( M-Q+1 .GT. I ) {
               clacgv(M-Q-I+1, X12(I,I), LDX12 );
               if ( M-Q == I ) {
                  clarfgp(M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) );
               } else {
                  clarfgp(M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) );
               }
            }
            X12(I,I) = ONE

            if ( I .LT. Q ) {
               clarf('R', P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X11(I+1,I+1), LDX11, WORK );
               clarf('R', M-P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X21(I+1,I+1), LDX21, WORK );
            }
            if ( P .GT. I ) {
               clarf('R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK );
            }
            if ( M-P .GT. I ) {
               clarf('R', M-P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(I+1,I), LDX22, WORK );
            }

            if (I .LT. Q) CALL CLACGV( Q-I, X11(I,I+1), LDX11 );
            clacgv(M-Q-I+1, X12(I,I), LDX12 );

         }

         // Reduce columns Q + 1, ..., P of X12, X22

         for (I = Q + 1; I <= P; I++) {

            cscal(M-Q-I+1, CMPLX( -Z1*Z4, 0.0E0 ), X12(I,I), LDX12 );
            clacgv(M-Q-I+1, X12(I,I), LDX12 );
            if ( I .GE. M-Q ) {
               clarfgp(M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) );
            } else {
               clarfgp(M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) );
            }
            X12(I,I) = ONE

            if ( P .GT. I ) {
               clarf('R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK );
            }
            if (M-P-Q .GE. 1) CALL CLARF( 'R', M-P-Q, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(Q+1,I), LDX22, WORK );

            clacgv(M-Q-I+1, X12(I,I), LDX12 );

         }

         // Reduce columns P + 1, ..., M - Q of X12, X22

         for (I = 1; I <= M - P - Q; I++) {

            cscal(M-P-Q-I+1, CMPLX( Z2*Z4, 0.0E0 ), X22(Q+I,P+I), LDX22 );
            clacgv(M-P-Q-I+1, X22(Q+I,P+I), LDX22 );
            clarfgp(M-P-Q-I+1, X22(Q+I,P+I), X22(Q+I,P+I+1), LDX22, TAUQ2(P+I) );
            X22(Q+I,P+I) = ONE
            clarf('R', M-P-Q-I, M-P-Q-I+1, X22(Q+I,P+I), LDX22, TAUQ2(P+I), X22(Q+I+1,P+I), LDX22, WORK );

            clacgv(M-P-Q-I+1, X22(Q+I,P+I), LDX22 );

         }

      } else {

         // Reduce columns 1, ..., Q of X11, X12, X21, X22

         for (I = 1; I <= Q; I++) {

            if ( I == 1 ) {
               cscal(P-I+1, CMPLX( Z1, 0.0E0 ), X11(I,I), LDX11 );
            } else {
               cscal(P-I+1, CMPLX( Z1*COS(PHI(I-1)), 0.0E0 ), X11(I,I), LDX11 );
               caxpy(P-I+1, CMPLX( -Z1*Z3*Z4*SIN(PHI(I-1)), 0.0E0 ), X12(I-1,I), LDX12, X11(I,I), LDX11 );
            }
            if ( I == 1 ) {
               cscal(M-P-I+1, CMPLX( Z2, 0.0E0 ), X21(I,I), LDX21 );
            } else {
               cscal(M-P-I+1, CMPLX( Z2*COS(PHI(I-1)), 0.0E0 ), X21(I,I), LDX21 );
               caxpy(M-P-I+1, CMPLX( -Z2*Z3*Z4*SIN(PHI(I-1)), 0.0E0 ), X22(I-1,I), LDX22, X21(I,I), LDX21 );
            }

            THETA(I) = ATAN2( SCNRM2( M-P-I+1, X21(I,I), LDX21 ), SCNRM2( P-I+1, X11(I,I), LDX11 ) )

            clacgv(P-I+1, X11(I,I), LDX11 );
            clacgv(M-P-I+1, X21(I,I), LDX21 );

            clarfgp(P-I+1, X11(I,I), X11(I,I+1), LDX11, TAUP1(I) );
            X11(I,I) = ONE
            if ( I == M-P ) {
               clarfgp(M-P-I+1, X21(I,I), X21(I,I), LDX21, TAUP2(I) );
            } else {
               clarfgp(M-P-I+1, X21(I,I), X21(I,I+1), LDX21, TAUP2(I) );
            }
            X21(I,I) = ONE

            clarf('R', Q-I, P-I+1, X11(I,I), LDX11, TAUP1(I), X11(I+1,I), LDX11, WORK );
            clarf('R', M-Q-I+1, P-I+1, X11(I,I), LDX11, TAUP1(I), X12(I,I), LDX12, WORK );
            clarf('R', Q-I, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X21(I+1,I), LDX21, WORK );
            clarf('R', M-Q-I+1, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X22(I,I), LDX22, WORK );

            clacgv(P-I+1, X11(I,I), LDX11 );
            clacgv(M-P-I+1, X21(I,I), LDX21 );

            if ( I .LT. Q ) {
               cscal(Q-I, CMPLX( -Z1*Z3*SIN(THETA(I)), 0.0E0 ), X11(I+1,I), 1 );
               caxpy(Q-I, CMPLX( Z2*Z3*COS(THETA(I)), 0.0E0 ), X21(I+1,I), 1, X11(I+1,I), 1 );
            }
            cscal(M-Q-I+1, CMPLX( -Z1*Z4*SIN(THETA(I)), 0.0E0 ), X12(I,I), 1 );
            caxpy(M-Q-I+1, CMPLX( Z2*Z4*COS(THETA(I)), 0.0E0 ), X22(I,I), 1, X12(I,I), 1 );

            if (I .LT. Q) PHI(I) = ATAN2( SCNRM2( Q-I, X11(I+1,I), 1 ), SCNRM2( M-Q-I+1, X12(I,I), 1 ) );

            if ( I .LT. Q ) {
               clarfgp(Q-I, X11(I+1,I), X11(I+2,I), 1, TAUQ1(I) );
               X11(I+1,I) = ONE
            }
            clarfgp(M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) );
            X12(I,I) = ONE

            if ( I .LT. Q ) {
               clarf('L', Q-I, P-I, X11(I+1,I), 1, CONJG(TAUQ1(I)), X11(I+1,I+1), LDX11, WORK );
               clarf('L', Q-I, M-P-I, X11(I+1,I), 1, CONJG(TAUQ1(I)), X21(I+1,I+1), LDX21, WORK );
            }
            clarf('L', M-Q-I+1, P-I, X12(I,I), 1, CONJG(TAUQ2(I)), X12(I,I+1), LDX12, WORK );

            if ( M-P .GT. I ) {
               clarf('L', M-Q-I+1, M-P-I, X12(I,I), 1, CONJG(TAUQ2(I)), X22(I,I+1), LDX22, WORK );
            }
         }

         // Reduce columns Q + 1, ..., P of X12, X22

         for (I = Q + 1; I <= P; I++) {

            cscal(M-Q-I+1, CMPLX( -Z1*Z4, 0.0E0 ), X12(I,I), 1 );
            clarfgp(M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) );
            X12(I,I) = ONE

            if ( P .GT. I ) {
               clarf('L', M-Q-I+1, P-I, X12(I,I), 1, CONJG(TAUQ2(I)), X12(I,I+1), LDX12, WORK );
            }
            if (M-P-Q .GE. 1) CALL CLARF( 'L', M-Q-I+1, M-P-Q, X12(I,I), 1, CONJG(TAUQ2(I)), X22(I,Q+1), LDX22, WORK );

         }

         // Reduce columns P + 1, ..., M - Q of X12, X22

         for (I = 1; I <= M - P - Q; I++) {

            cscal(M-P-Q-I+1, CMPLX( Z2*Z4, 0.0E0 ), X22(P+I,Q+I), 1 );
            clarfgp(M-P-Q-I+1, X22(P+I,Q+I), X22(P+I+1,Q+I), 1, TAUQ2(P+I) );
            X22(P+I,Q+I) = ONE
            if ( M-P-Q != I ) {
               clarf('L', M-P-Q-I+1, M-P-Q-I, X22(P+I,Q+I), 1, CONJG(TAUQ2(P+I)), X22(P+I,Q+I+1), LDX22, WORK );
            }
         }

      }

      RETURN

      // End of CUNBDB

      }
