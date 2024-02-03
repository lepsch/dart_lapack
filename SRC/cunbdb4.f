*> \verbatim
*>
*>  The upper-bidiagonal blocks B11, B21 are represented implicitly by
*>  angles THETA(1), ..., THETA(Q) and PHI(1), ..., PHI(Q-1). Every entry
*>  in each bidiagonal band is a product of a sine or cosine of a THETA
*>  with a sine or cosine of a PHI. See [1] or CUNCSD for details.
*>
*>  P1, P2, and Q1 are represented as products of elementary reflectors.
*>  See CUNCSD2BY1 for details on generating P1, P2, and Q1 using CUNGQR
*>  and CUNGLQ.
*> \endverbatim

*> \par References:
*  ================
*>
*>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
*>      Algorithms, 50(1):33-65, 2009.
*>
*  =====================================================================
      SUBROUTINE CUNBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      // ..
      // .. Array Arguments ..
      REAL               PHI(*), THETA(*)
      COMPLEX            PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*)
      // ..

*  ====================================================================

      // .. Parameters ..
      COMPLEX            NEGONE, ONE, ZERO
      const              NEGONE = (-1.0E0,0.0E0), ONE = (1.0E0,0.0E0), ZERO = (0.0E0,0.0E0) ;
      // ..
      // .. Local Scalars ..
      REAL               C, S
      int                CHILDINFO, I, ILARF, IORBDB5, J, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, CLARFGP, CUNBDB5, CSROT, CSCAL, CLACGV, XERBLA
      // ..
      // .. External Functions ..
      REAL               SCNRM2, SROUNDUP_LWORK
      // EXTERNAL SCNRM2, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      LQUERY = LWORK == -1

      if ( M .LT. 0 ) {
         INFO = -1
      } else if ( P .LT. M-Q .OR. M-P .LT. M-Q ) {
         INFO = -2
      } else if ( Q .LT. M-Q .OR. Q .GT. M ) {
         INFO = -3
      } else if ( LDX11 .LT. MAX( 1, P ) ) {
         INFO = -5
      } else if ( LDX21 .LT. MAX( 1, M-P ) ) {
         INFO = -7
      }

      // Compute workspace

      if ( INFO == 0 ) {
         ILARF = 2
         LLARF = MAX( Q-1, P-1, M-P-1 )
         IORBDB5 = 2
         LORBDB5 = Q
         LWORKOPT = ILARF + LLARF - 1
         LWORKOPT = MAX( LWORKOPT, IORBDB5 + LORBDB5 - 1 )
         LWORKMIN = LWORKOPT
         WORK(1) = SROUNDUP_LWORK(LWORKOPT)
         if ( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) {
           INFO = -14
         }
      }
      if ( INFO != 0 ) {
         xerbla('CUNBDB4', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Reduce columns 1, ..., M-Q of X11 and X21

      for (I = 1; I <= M-Q; I++) {

         if ( I == 1 ) {
            for (J = 1; J <= M; J++) {
               PHANTOM(J) = ZERO
            }
            cunbdb5(P, M-P, Q, PHANTOM(1), 1, PHANTOM(P+1), 1, X11, LDX11, X21, LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
            cscal(P, NEGONE, PHANTOM(1), 1 );
            clarfgp(P, PHANTOM(1), PHANTOM(2), 1, TAUP1(1) );
            clarfgp(M-P, PHANTOM(P+1), PHANTOM(P+2), 1, TAUP2(1) );
            THETA(I) = ATAN2( REAL( PHANTOM(1) ), REAL( PHANTOM(P+1) ) )
            C = COS( THETA(I) )
            S = SIN( THETA(I) )
            PHANTOM(1) = ONE
            PHANTOM(P+1) = ONE
            clarf('L', P, Q, PHANTOM(1), 1, CONJG(TAUP1(1)), X11, LDX11, WORK(ILARF) );
            clarf('L', M-P, Q, PHANTOM(P+1), 1, CONJG(TAUP2(1)), X21, LDX21, WORK(ILARF) );
         } else {
            cunbdb5(P-I+1, M-P-I+1, Q-I+1, X11(I,I-1), 1, X21(I,I-1), 1, X11(I,I), LDX11, X21(I,I), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO );
            cscal(P-I+1, NEGONE, X11(I,I-1), 1 );
            clarfgp(P-I+1, X11(I,I-1), X11(I+1,I-1), 1, TAUP1(I) );
            clarfgp(M-P-I+1, X21(I,I-1), X21(I+1,I-1), 1, TAUP2(I) );
            THETA(I) = ATAN2( REAL( X11(I,I-1) ), REAL( X21(I,I-1) ) )
            C = COS( THETA(I) )
            S = SIN( THETA(I) )
            X11(I,I-1) = ONE
            X21(I,I-1) = ONE
            clarf('L', P-I+1, Q-I+1, X11(I,I-1), 1, CONJG(TAUP1(I)), X11(I,I), LDX11, WORK(ILARF) );
            clarf('L', M-P-I+1, Q-I+1, X21(I,I-1), 1, CONJG(TAUP2(I)), X21(I,I), LDX21, WORK(ILARF) );
         }

         csrot(Q-I+1, X11(I,I), LDX11, X21(I,I), LDX21, S, -C );
         clacgv(Q-I+1, X21(I,I), LDX21 );
         clarfgp(Q-I+1, X21(I,I), X21(I,I+1), LDX21, TAUQ1(I) );
         C = REAL( X21(I,I) )
         X21(I,I) = ONE
         clarf('R', P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) );
         clarf('R', M-P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X21(I+1,I), LDX21, WORK(ILARF) );
         clacgv(Q-I+1, X21(I,I), LDX21 );
         if ( I .LT. M-Q ) {
            S = SQRT( SCNRM2( P-I, X11(I+1,I), 1 )**2 + SCNRM2( M-P-I, X21(I+1,I), 1 )**2 )
            PHI(I) = ATAN2( S, C )
         }

      }

      // Reduce the bottom-right portion of X11 to [ I 0 ]

      for (I = M - Q + 1; I <= P; I++) {
         clacgv(Q-I+1, X11(I,I), LDX11 );
         clarfgp(Q-I+1, X11(I,I), X11(I,I+1), LDX11, TAUQ1(I) );
         X11(I,I) = ONE
         clarf('R', P-I, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) );
         clarf('R', Q-P, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X21(M-Q+1,I), LDX21, WORK(ILARF) );
         clacgv(Q-I+1, X11(I,I), LDX11 );
      }

      // Reduce the bottom-right portion of X21 to [ 0 I ]

      for (I = P + 1; I <= Q; I++) {
         clacgv(Q-I+1, X21(M-Q+I-P,I), LDX21 );
         clarfgp(Q-I+1, X21(M-Q+I-P,I), X21(M-Q+I-P,I+1), LDX21, TAUQ1(I) );
         X21(M-Q+I-P,I) = ONE
         clarf('R', Q-I, Q-I+1, X21(M-Q+I-P,I), LDX21, TAUQ1(I), X21(M-Q+I-P+1,I), LDX21, WORK(ILARF) );
         clacgv(Q-I+1, X21(M-Q+I-P,I), LDX21 );
      }

      RETURN

      // End of CUNBDB4

      }
