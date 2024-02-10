      void cunbdb2(final int M, final int P, final int Q, final Matrix<double> X11, final int LDX11, final Matrix<double> X21, final int LDX21, final int THETA, final int PHI, final int TAUP1, final int TAUP2, final int TAUQ1, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      double               PHI(*), THETA(*);
      Complex            TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*);
      // ..

// ====================================================================

      // .. Parameters ..
      Complex            NEGONE, ONE;
      const              NEGONE = (-1.0,0.0), ONE = (1.0,0.0) ;
      double               C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, CLARFGP, CUNBDB5, CSROT, CSCAL, CLACGV, XERBLA
      // ..
      // .. External Functions ..
      //- REAL               SCNRM2, SROUNDUP_LWORK;
      // EXTERNAL SCNRM2, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT

      // Test input arguments

      INFO = 0;
      LQUERY = LWORK == -1;

      if ( M < 0 ) {
         INFO = -1;
      } else if ( P < 0 || P > M-P ) {
         INFO = -2;
      } else if ( Q < 0 || Q < P || M-Q < P ) {
         INFO = -3;
      } else if ( LDX11 < max( 1, P ) ) {
         INFO = -5;
      } else if ( LDX21 < max( 1, M-P ) ) {
         INFO = -7;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         ILARF = 2;
         LLARF = max( P-1, M-P, Q-1 );
         IORBDB5 = 2;
         LORBDB5 = Q-1;
         LWORKOPT = max( ILARF+LLARF-1, IORBDB5+LORBDB5-1 );
         LWORKMIN = LWORKOPT;
         WORK[1] = SROUNDUP_LWORK(LWORKOPT);
         if ( LWORK < LWORKMIN && !LQUERY ) {
           INFO = -14;
         }
      }
      if ( INFO != 0 ) {
         xerbla('CUNBDB2', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Reduce rows 1, ..., P of X11 and X21

      for (I = 1; I <= P; I++) {

         if ( I > 1 ) {
            csrot(Q-I+1, X11(I,I), LDX11, X21(I-1,I), LDX21, C, S );
         }
         clacgv(Q-I+1, X11(I,I), LDX11 );
         clarfgp(Q-I+1, X11(I,I), X11(I,I+1), LDX11, TAUQ1(I) );
         C = double( X11(I,I) );
         X11[I][I] = ONE;
         clarf('R', P-I, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) );
         clarf('R', M-P-I+1, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X21(I,I), LDX21, WORK(ILARF) );
         clacgv(Q-I+1, X11(I,I), LDX11 );
         S = sqrt( SCNRM2( P-I, X11(I+1,I), 1 )**2 + SCNRM2( M-P-I+1, X21(I,I), 1 )**2 );
         THETA[I] = ATAN2( S, C );

         cunbdb5(P-I, M-P-I+1, Q-I, X11(I+1,I), 1, X21(I,I), 1, X11(I+1,I+1), LDX11, X21(I,I+1), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
         cscal(P-I, NEGONE, X11(I+1,I), 1 );
         clarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
         if ( I < P ) {
            clarfgp(P-I, X11(I+1,I), X11(I+2,I), 1, TAUP1(I) );
            PHI[I] = ATAN2( REAL( X11(I+1,I) ), double( X21(I,I) ) );
            C = COS( PHI(I) );
            S = SIN( PHI(I) );
            X11[I+1][I] = ONE;
            clarf('L', P-I, Q-I, X11(I+1,I), 1, CONJG(TAUP1(I)), X11(I+1,I+1), LDX11, WORK(ILARF) );
         }
         X21[I][I] = ONE;
         clarf('L', M-P-I+1, Q-I, X21(I,I), 1, CONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK(ILARF) );

      }

      // Reduce the bottom-right portion of X21 to the identity matrix

      for (I = P + 1; I <= Q; I++) {
         clarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
         X21[I][I] = ONE;
         clarf('L', M-P-I+1, Q-I, X21(I,I), 1, CONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK(ILARF) );
      }

      }
