      void zunbdb4(M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      // ..
      // .. Array Arguments ..
      double             PHI(*), THETA(*);
      Complex         PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*);
      // ..

// ====================================================================

      // .. Parameters ..
      Complex         NEGONE, ONE, ZERO;
      const              NEGONE = (-1.0,0.0), ONE = (1.0,0.0), ZERO = (0.0,0.0) ;
      // ..
      // .. Local Scalars ..
      double             C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, J, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARF, ZLARFGP, ZUNBDB5, ZDROT, ZSCAL, ZLACGV, XERBLA
      // ..
      // .. External Functions ..
      //- double             DZNRM2;
      // EXTERNAL DZNRM2
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0;
      LQUERY = LWORK == -1;

      if ( M < 0 ) {
         INFO = -1;
      } else if ( P < M-Q || M-P < M-Q ) {
         INFO = -2;
      } else if ( Q < M-Q || Q > M ) {
         INFO = -3;
      } else if ( LDX11 < max( 1, P ) ) {
         INFO = -5;
      } else if ( LDX21 < max( 1, M-P ) ) {
         INFO = -7;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         ILARF = 2;
         LLARF = max( Q-1, P-1, M-P-1 );
         IORBDB5 = 2;
         LORBDB5 = Q;
         LWORKOPT = ILARF + LLARF - 1;
         LWORKOPT = max( LWORKOPT, IORBDB5 + LORBDB5 - 1 );
         LWORKMIN = LWORKOPT;
         WORK[1] = LWORKOPT;
         if ( LWORK < LWORKMIN && !LQUERY ) {
           INFO = -14;
         }
      }
      if ( INFO != 0 ) {
         xerbla('ZUNBDB4', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Reduce columns 1, ..., M-Q of X11 and X21

      for (I = 1; I <= M-Q; I++) {

         if ( I == 1 ) {
            for (J = 1; J <= M; J++) {
               PHANTOM[J] = ZERO;
            }
            zunbdb5(P, M-P, Q, PHANTOM(1), 1, PHANTOM(P+1), 1, X11, LDX11, X21, LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
            zscal(P, NEGONE, PHANTOM(1), 1 );
            zlarfgp(P, PHANTOM(1), PHANTOM(2), 1, TAUP1(1) );
            zlarfgp(M-P, PHANTOM(P+1), PHANTOM(P+2), 1, TAUP2(1) );
            THETA[I] = ATAN2( (PHANTOM(1)).toDouble(), (PHANTOM(P+1)).toDouble() );
            C = COS( THETA(I) );
            S = SIN( THETA(I) );
            PHANTOM[1] = ONE;
            PHANTOM[P+1] = ONE;
            zlarf('L', P, Q, PHANTOM(1), 1, DCONJG(TAUP1(1)), X11, LDX11, WORK(ILARF) );
            zlarf('L', M-P, Q, PHANTOM(P+1), 1, DCONJG(TAUP2(1)), X21, LDX21, WORK(ILARF) );
         } else {
            zunbdb5(P-I+1, M-P-I+1, Q-I+1, X11(I,I-1), 1, X21(I,I-1), 1, X11(I,I), LDX11, X21(I,I), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
            zscal(P-I+1, NEGONE, X11(I,I-1), 1 );
            zlarfgp(P-I+1, X11(I,I-1), X11(I+1,I-1), 1, TAUP1(I) );
            zlarfgp(M-P-I+1, X21(I,I-1), X21(I+1,I-1), 1, TAUP2(I) );
            THETA[I] = ATAN2( (X11(I,I-1)).toDouble(), (X21(I,I-1)).toDouble() );
            C = COS( THETA(I) );
            S = SIN( THETA(I) );
            X11[I,I-1] = ONE;
            X21[I,I-1] = ONE;
            zlarf('L', P-I+1, Q-I+1, X11(I,I-1), 1, DCONJG(TAUP1(I)), X11(I,I), LDX11, WORK(ILARF) );
            zlarf('L', M-P-I+1, Q-I+1, X21(I,I-1), 1, DCONJG(TAUP2(I)), X21(I,I), LDX21, WORK(ILARF) );
         }

         zdrot(Q-I+1, X11(I,I), LDX11, X21(I,I), LDX21, S, -C );
         zlacgv(Q-I+1, X21(I,I), LDX21 );
         zlarfgp(Q-I+1, X21(I,I), X21(I,I+1), LDX21, TAUQ1(I) );
         C = (X21(I,I)).toDouble();
         X21[I,I] = ONE;
         zlarf('R', P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) );
         zlarf('R', M-P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X21(I+1,I), LDX21, WORK(ILARF) );
         zlacgv(Q-I+1, X21(I,I), LDX21 );
         if ( I < M-Q ) {
            S = sqrt( DZNRM2( P-I, X11(I+1,I), 1 )**2 + DZNRM2( M-P-I, X21(I+1,I), 1 )**2 );
            PHI[I] = ATAN2( S, C );
         }

      }

      // Reduce the bottom-right portion of X11 to [ I 0 ]

      for (I = M - Q + 1; I <= P; I++) {
         zlacgv(Q-I+1, X11(I,I), LDX11 );
         zlarfgp(Q-I+1, X11(I,I), X11(I,I+1), LDX11, TAUQ1(I) );
         X11[I,I] = ONE;
         zlarf('R', P-I, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) );
         zlarf('R', Q-P, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X21(M-Q+1,I), LDX21, WORK(ILARF) );
         zlacgv(Q-I+1, X11(I,I), LDX11 );
      }

      // Reduce the bottom-right portion of X21 to [ 0 I ]

      for (I = P + 1; I <= Q; I++) {
         zlacgv(Q-I+1, X21(M-Q+I-P,I), LDX21 );
         zlarfgp(Q-I+1, X21(M-Q+I-P,I), X21(M-Q+I-P,I+1), LDX21, TAUQ1(I) );
         X21[M-Q+I-P,I] = ONE;
         zlarf('R', Q-I, Q-I+1, X21(M-Q+I-P,I), LDX21, TAUQ1(I), X21(M-Q+I-P+1,I), LDX21, WORK(ILARF) );
         zlacgv(Q-I+1, X21(M-Q+I-P,I), LDX21 );
      }

      return;
      }