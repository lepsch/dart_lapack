      SUBROUTINE CUNBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      // ..
      // .. Array Arguments ..
      REAL               PHI(*), THETA(*)
      COMPLEX            TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*)
      // ..

*  ====================================================================

      // .. Parameters ..
      COMPLEX            NEGONE, ONE
      const              NEGONE = (-1.0E0,0.0E0), ONE = (1.0E0,0.0E0) ;
      // ..
      // .. Local Scalars ..
      REAL               C, S
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
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
      LQUERY = LWORK .EQ. -1

      if ( M .LT. 0 ) {
         INFO = -1
      } else if ( P .LT. 0 .OR. P .GT. M-P ) {
         INFO = -2
      } else if ( Q .LT. 0 .OR. Q .LT. P .OR. M-Q .LT. P ) {
         INFO = -3
      } else if ( LDX11 .LT. MAX( 1, P ) ) {
         INFO = -5
      } else if ( LDX21 .LT. MAX( 1, M-P ) ) {
         INFO = -7
      }

      // Compute workspace

      if ( INFO .EQ. 0 ) {
         ILARF = 2
         LLARF = MAX( P-1, M-P, Q-1 )
         IORBDB5 = 2
         LORBDB5 = Q-1
         LWORKOPT = MAX( ILARF+LLARF-1, IORBDB5+LORBDB5-1 )
         LWORKMIN = LWORKOPT
         WORK(1) = SROUNDUP_LWORK(LWORKOPT)
         if ( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) {
           INFO = -14
         }
      }
      if ( INFO .NE. 0 ) {
         CALL XERBLA( 'CUNBDB2', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Reduce rows 1, ..., P of X11 and X21

      DO I = 1, P

         if ( I .GT. 1 ) {
            CALL CSROT( Q-I+1, X11(I,I), LDX11, X21(I-1,I), LDX21, C, S )
         }
         CALL CLACGV( Q-I+1, X11(I,I), LDX11 )
         CALL CLARFGP( Q-I+1, X11(I,I), X11(I,I+1), LDX11, TAUQ1(I) )
         C = REAL( X11(I,I) )
         X11(I,I) = ONE
         CALL CLARF( 'R', P-I, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) )          CALL CLARF( 'R', M-P-I+1, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X21(I,I), LDX21, WORK(ILARF) )
         CALL CLACGV( Q-I+1, X11(I,I), LDX11 )
         S = SQRT( SCNRM2( P-I, X11(I+1,I), 1 )**2 + SCNRM2( M-P-I+1, X21(I,I), 1 )**2 )
         THETA(I) = ATAN2( S, C )

         CALL CUNBDB5( P-I, M-P-I+1, Q-I, X11(I+1,I), 1, X21(I,I), 1, X11(I+1,I+1), LDX11, X21(I,I+1), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO )
         CALL CSCAL( P-I, NEGONE, X11(I+1,I), 1 )
         CALL CLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
         if ( I .LT. P ) {
            CALL CLARFGP( P-I, X11(I+1,I), X11(I+2,I), 1, TAUP1(I) )
            PHI(I) = ATAN2( REAL( X11(I+1,I) ), REAL( X21(I,I) ) )
            C = COS( PHI(I) )
            S = SIN( PHI(I) )
            X11(I+1,I) = ONE
            CALL CLARF( 'L', P-I, Q-I, X11(I+1,I), 1, CONJG(TAUP1(I)), X11(I+1,I+1), LDX11, WORK(ILARF) )
         }
         X21(I,I) = ONE
         CALL CLARF( 'L', M-P-I+1, Q-I, X21(I,I), 1, CONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK(ILARF) )

      END DO

      // Reduce the bottom-right portion of X21 to the identity matrix

      DO I = P + 1, Q
         CALL CLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
         X21(I,I) = ONE
         CALL CLARF( 'L', M-P-I+1, Q-I, X21(I,I), 1, CONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK(ILARF) )
      END DO

      RETURN

      // End of CUNBDB2

      }
