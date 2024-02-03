      SUBROUTINE DORBDB3( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      // ..
      // .. Array Arguments ..
      double             PHI(*), THETA(*);
      double             TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*);
      // ..

*  ====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      double             C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARF, DLARFGP, DORBDB5, DROT, XERBLA
      // ..
      // .. External Functions ..
      double             DNRM2;
      // EXTERNAL DNRM2
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      LQUERY = LWORK .EQ. -1

      if ( M .LT. 0 ) {
         INFO = -1
      } else if ( 2*P .LT. M .OR. P .GT. M ) {
         INFO = -2
      } else if ( Q .LT. M-P .OR. M-Q .LT. M-P ) {
         INFO = -3
      } else if ( LDX11 .LT. MAX( 1, P ) ) {
         INFO = -5
      } else if ( LDX21 .LT. MAX( 1, M-P ) ) {
         INFO = -7
      }

      // Compute workspace

      if ( INFO .EQ. 0 ) {
         ILARF = 2
         LLARF = MAX( P, M-P-1, Q-1 )
         IORBDB5 = 2
         LORBDB5 = Q-1
         LWORKOPT = MAX( ILARF+LLARF-1, IORBDB5+LORBDB5-1 )
         LWORKMIN = LWORKOPT
         WORK(1) = LWORKOPT
         if ( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) {
           INFO = -14
         }
      }
      if ( INFO .NE. 0 ) {
         CALL XERBLA( 'DORBDB3', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Reduce rows 1, ..., M-P of X11 and X21

      DO I = 1, M-P

         if ( I .GT. 1 ) {
            CALL DROT( Q-I+1, X11(I-1,I), LDX11, X21(I,I), LDX11, C, S )
         }

         CALL DLARFGP( Q-I+1, X21(I,I), X21(I,I+1), LDX21, TAUQ1(I) )
         S = X21(I,I)
         X21(I,I) = ONE
         CALL DLARF( 'R', P-I+1, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X11(I,I), LDX11, WORK(ILARF) )          CALL DLARF( 'R', M-P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X21(I+1,I), LDX21, WORK(ILARF) )          C = SQRT( DNRM2( P-I+1, X11(I,I), 1 )**2 + DNRM2( M-P-I, X21(I+1,I), 1 )**2 )
         THETA(I) = ATAN2( S, C )

         CALL DORBDB5( P-I+1, M-P-I, Q-I, X11(I,I), 1, X21(I+1,I), 1, X11(I,I+1), LDX11, X21(I+1,I+1), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO )
         CALL DLARFGP( P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) )
         if ( I .LT. M-P ) {
            CALL DLARFGP( M-P-I, X21(I+1,I), X21(I+2,I), 1, TAUP2(I) )
            PHI(I) = ATAN2( X21(I+1,I), X11(I,I) )
            C = COS( PHI(I) )
            S = SIN( PHI(I) )
            X21(I+1,I) = ONE
            CALL DLARF( 'L', M-P-I, Q-I, X21(I+1,I), 1, TAUP2(I), X21(I+1,I+1), LDX21, WORK(ILARF) )
         }
         X11(I,I) = ONE
         CALL DLARF( 'L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK(ILARF) )

      END DO

      // Reduce the bottom-right portion of X11 to the identity matrix

      DO I = M-P + 1, Q
         CALL DLARFGP( P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) )
         X11(I,I) = ONE
         CALL DLARF( 'L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK(ILARF) )
      END DO

      RETURN

      // End of DORBDB3

      }
