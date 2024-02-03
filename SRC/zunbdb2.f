      SUBROUTINE ZUNBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      // ..
      // .. Array Arguments ..
      double             PHI(*), THETA(*);
      COMPLEX*16         TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*);
      // ..

// ====================================================================

      // .. Parameters ..
      COMPLEX*16         NEGONE, ONE;
      const              NEGONE = (-1.0,0.0), ONE = (1.0,0.0) ;
      // ..
      // .. Local Scalars ..
      double             C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARF, ZLARFGP, ZUNBDB5, ZDROT, ZSCAL, ZLACGV, XERBLA
      // ..
      // .. External Functions ..
      double             DZNRM2;
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
         WORK(1) = LWORKOPT;
         if ( LWORK < LWORKMIN && !LQUERY ) {
           INFO = -14;
         }
      }
      if ( INFO != 0 ) {
         xerbla('ZUNBDB2', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Reduce rows 1, ..., P of X11 and X21

      for (I = 1; I <= P; I++) {

         if ( I > 1 ) {
            zdrot(Q-I+1, X11(I,I), LDX11, X21(I-1,I), LDX21, C, S );
         }
         zlacgv(Q-I+1, X11(I,I), LDX11 );
         zlarfgp(Q-I+1, X11(I,I), X11(I,I+1), LDX11, TAUQ1(I) );
         C = DBLE( X11(I,I) );
         X11(I,I) = ONE;
         zlarf('R', P-I, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) );
         zlarf('R', M-P-I+1, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X21(I,I), LDX21, WORK(ILARF) );
         zlacgv(Q-I+1, X11(I,I), LDX11 );
         S = sqrt( DZNRM2( P-I, X11(I+1,I), 1 )**2 + DZNRM2( M-P-I+1, X21(I,I), 1 )**2 );
         THETA(I) = ATAN2( S, C );

         zunbdb5(P-I, M-P-I+1, Q-I, X11(I+1,I), 1, X21(I,I), 1, X11(I+1,I+1), LDX11, X21(I,I+1), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
         zscal(P-I, NEGONE, X11(I+1,I), 1 );
         zlarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
         if ( I < P ) {
            zlarfgp(P-I, X11(I+1,I), X11(I+2,I), 1, TAUP1(I) );
            PHI(I) = ATAN2( DBLE( X11(I+1,I) ), DBLE( X21(I,I) ) );
            C = COS( PHI(I) );
            S = SIN( PHI(I) );
            X11(I+1,I) = ONE;
            zlarf('L', P-I, Q-I, X11(I+1,I), 1, DCONJG(TAUP1(I)), X11(I+1,I+1), LDX11, WORK(ILARF) );
         }
         X21(I,I) = ONE;
         zlarf('L', M-P-I+1, Q-I, X21(I,I), 1, DCONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK(ILARF) );

      }

      // Reduce the bottom-right portion of X21 to the identity matrix

      for (I = P + 1; I <= Q; I++) {
         zlarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
         X21(I,I) = ONE;
         zlarf('L', M-P-I+1, Q-I, X21(I,I), 1, DCONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK(ILARF) );
      }

      return;

      // End of ZUNBDB2

      }
