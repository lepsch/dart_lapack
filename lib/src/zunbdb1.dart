      void zunbdb1(final int M, final int P, final int Q, final Matrix<double> X11, final int LDX11, final Matrix<double> X21, final int LDX21, final int THETA, final int PHI, final int TAUP1, final int TAUP2, final int TAUQ1, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      double             PHI(*), THETA(*);
      Complex         TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*);
      // ..

// ====================================================================

      // .. Parameters ..
      Complex         ONE;
      const              ONE = (1.0,0.0) ;
      double             C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARF, ZLARFGP, ZUNBDB5, ZDROT, XERBLA
      // EXTERNAL ZLACGV
      // ..
      // .. External Functions ..
      //- double             DZNRM2;
      // EXTERNAL DZNRM2
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT

      // Test input arguments

      INFO = 0;
      LQUERY = LWORK == -1;

      if ( M < 0 ) {
         INFO = -1;
      } else if ( P < Q || M-P < Q ) {
         INFO = -2;
      } else if ( Q < 0 || M-Q < Q ) {
         INFO = -3;
      } else if ( LDX11 < max( 1, P ) ) {
         INFO = -5;
      } else if ( LDX21 < max( 1, M-P ) ) {
         INFO = -7;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         ILARF = 2;
         LLARF = max( P-1, M-P-1, Q-1 );
         IORBDB5 = 2;
         LORBDB5 = Q-2;
         LWORKOPT = max( ILARF+LLARF-1, IORBDB5+LORBDB5-1 );
         LWORKMIN = LWORKOPT;
         WORK[1] = LWORKOPT;
         if ( LWORK < LWORKMIN && !LQUERY ) {
           INFO = -14;
         }
      }
      if ( INFO != 0 ) {
         xerbla('ZUNBDB1', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Reduce columns 1, ..., Q of X11 and X21

      for (I = 1; I <= Q; I++) {

         zlarfgp(P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) );
         zlarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
         THETA[I] = ATAN2( (X21(I,I)).toDouble(), (X11(I,I)).toDouble() );
         C = COS( THETA(I) );
         S = SIN( THETA(I) );
         X11[I][I] = ONE;
         X21[I][I] = ONE;
         zlarf('L', P-I+1, Q-I, X11(I,I), 1, DCONJG(TAUP1(I)), X11(I,I+1), LDX11, WORK(ILARF) );
         zlarf('L', M-P-I+1, Q-I, X21(I,I), 1, DCONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK(ILARF) );

         if ( I < Q ) {
            zdrot(Q-I, X11(I,I+1), LDX11, X21(I,I+1), LDX21, C, S );
            zlacgv(Q-I, X21(I,I+1), LDX21 );
            zlarfgp(Q-I, X21(I,I+1), X21(I,I+2), LDX21, TAUQ1(I) );
            S = (X21(I,I+1)).toDouble();
            X21[I][I+1] = ONE;
            zlarf('R', P-I, Q-I, X21(I,I+1), LDX21, TAUQ1(I), X11(I+1,I+1), LDX11, WORK(ILARF) );
            zlarf('R', M-P-I, Q-I, X21(I,I+1), LDX21, TAUQ1(I), X21(I+1,I+1), LDX21, WORK(ILARF) );
            zlacgv(Q-I, X21(I,I+1), LDX21 );
            C = sqrt( DZNRM2( P-I, X11(I+1,I+1), 1 )**2 + DZNRM2( M-P-I, X21(I+1,I+1), 1 )**2 );
            PHI[I] = ATAN2( S, C );
            zunbdb5(P-I, M-P-I, Q-I-1, X11(I+1,I+1), 1, X21(I+1,I+1), 1, X11(I+1,I+2), LDX11, X21(I+1,I+2), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
         }

      }

      }
