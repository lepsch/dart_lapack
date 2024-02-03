      SUBROUTINE SORBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      // ..
      // .. Array Arguments ..
      REAL               PHI(*), THETA(*)
      REAL               TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*)
      // ..

*  ====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      REAL               C, S
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, SLARFGP, SORBDB5, SROT, XERBLA
      // ..
      // .. External Functions ..
      REAL               SNRM2
      // EXTERNAL SNRM2
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      LQUERY = LWORK == -1

      if ( M < 0 ) {
         INFO = -1
      } else if ( P < Q || M-P < Q ) {
         INFO = -2
      } else if ( Q < 0 || M-Q < Q ) {
         INFO = -3
      } else if ( LDX11 < MAX( 1, P ) ) {
         INFO = -5
      } else if ( LDX21 < MAX( 1, M-P ) ) {
         INFO = -7
      }

      // Compute workspace

      if ( INFO == 0 ) {
         ILARF = 2
         LLARF = MAX( P-1, M-P-1, Q-1 )
         IORBDB5 = 2
         LORBDB5 = Q-2
         LWORKOPT = MAX( ILARF+LLARF-1, IORBDB5+LORBDB5-1 )
         LWORKMIN = LWORKOPT
         WORK(1) = LWORKOPT
         if ( LWORK < LWORKMIN && .NOT.LQUERY ) {
           INFO = -14
         }
      }
      if ( INFO != 0 ) {
         xerbla('SORBDB1', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Reduce columns 1, ..., Q of X11 and X21

      for (I = 1; I <= Q; I++) {

         slarfgp(P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) );
         slarfgp(M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) );
         THETA(I) = ATAN2( X21(I,I), X11(I,I) )
         C = COS( THETA(I) )
         S = SIN( THETA(I) )
         X11(I,I) = ONE
         X21(I,I) = ONE
         slarf('L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK(ILARF) );
         slarf('L', M-P-I+1, Q-I, X21(I,I), 1, TAUP2(I), X21(I,I+1), LDX21, WORK(ILARF) );

         if ( I < Q ) {
            srot(Q-I, X11(I,I+1), LDX11, X21(I,I+1), LDX21, C, S );
            slarfgp(Q-I, X21(I,I+1), X21(I,I+2), LDX21, TAUQ1(I) );
            S = X21(I,I+1)
            X21(I,I+1) = ONE
            slarf('R', P-I, Q-I, X21(I,I+1), LDX21, TAUQ1(I), X11(I+1,I+1), LDX11, WORK(ILARF) );
            slarf('R', M-P-I, Q-I, X21(I,I+1), LDX21, TAUQ1(I), X21(I+1,I+1), LDX21, WORK(ILARF) )             C = SQRT( SNRM2( P-I, X11(I+1,I+1), 1 )**2 + SNRM2( M-P-I, X21(I+1,I+1), 1 )**2 );
            PHI(I) = ATAN2( S, C )
            sorbdb5(P-I, M-P-I, Q-I-1, X11(I+1,I+1), 1, X21(I+1,I+1), 1, X11(I+1,I+2), LDX11, X21(I+1,I+2), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
         }

      }

      RETURN

      // End of SORBDB1

      }
