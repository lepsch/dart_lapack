      SUBROUTINE SORCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11, X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, WORK, LWORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU1, JOBU2, JOBV1T;
      int                INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21, M, P, Q;
      // ..
      // .. Array Arguments ..
      REAL               THETA(*)
      REAL               U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*), X11(LDX11,*), X21(LDX21,*)
      int                IWORK(*);
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E0, ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                CHILDINFO, I, IB11D, IB11E, IB12D, IB12E, IB21D, IB21E, IB22D, IB22E, IBBCSD, IORBDB, IORGLQ, IORGQR, IPHI, ITAUP1, ITAUP2, ITAUQ1, J, LBBCSD, LORBDB, LORGLQ, LORGLQMIN, LORGLQOPT, LORGQR, LORGQRMIN, LORGQROPT, LWORKMIN, LWORKOPT, R;
      bool               LQUERY, WANTU1, WANTU2, WANTV1T;
      // ..
      // .. Local Arrays ..
      REAL               DUM1(1), DUM2(1,1)
      // ..
      // .. External Subroutines ..
      // EXTERNAL SBBCSD, SCOPY, SLACPY, SLAPMR, SLAPMT, SORBDB1, SORBDB2, SORBDB3, SORBDB4, SORGLQ, SORGQR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTV1T = LSAME( JOBV1T, 'Y' )
      LQUERY = LWORK .EQ. -1

      if ( M .LT. 0 ) {
         INFO = -4
      } else if ( P .LT. 0 .OR. P .GT. M ) {
         INFO = -5
      } else if ( Q .LT. 0 .OR. Q .GT. M ) {
         INFO = -6
      } else if ( LDX11 .LT. MAX( 1, P ) ) {
         INFO = -8
      } else if ( LDX21 .LT. MAX( 1, M-P ) ) {
         INFO = -10
      } else if ( WANTU1 .AND. LDU1 .LT. MAX( 1, P ) ) {
         INFO = -13
      } else if ( WANTU2 .AND. LDU2 .LT. MAX( 1, M - P ) ) {
         INFO = -15
      } else if ( WANTV1T .AND. LDV1T .LT. MAX( 1, Q ) ) {
         INFO = -17
      }

      R = MIN( P, M-P, Q, M-Q )

      // Compute workspace

        // WORK layout:
      // |-------------------------------------------------------|
      // | LWORKOPT (1)                                          |
      // |-------------------------------------------------------|
      // | PHI (MAX(1,R-1))                                      |
      // |-------------------------------------------------------|
      // | TAUP1 (MAX(1,P))                        | B11D (R)    |
      // | TAUP2 (MAX(1,M-P))                      | B11E (R-1)  |
      // | TAUQ1 (MAX(1,Q))                        | B12D (R)    |
      // |-----------------------------------------| B12E (R-1)  |
      // | SORBDB WORK | SORGQR WORK | SORGLQ WORK | B21D (R)    |
      // |             |             |             | B21E (R-1)  |
      // |             |             |             | B22D (R)    |
      // |             |             |             | B22E (R-1)  |
      // |             |             |             | SBBCSD WORK |
      // |-------------------------------------------------------|

      if ( INFO .EQ. 0 ) {
         IPHI = 2
         IB11D = IPHI + MAX( 1, R-1 )
         IB11E = IB11D + MAX( 1, R )
         IB12D = IB11E + MAX( 1, R - 1 )
         IB12E = IB12D + MAX( 1, R )
         IB21D = IB12E + MAX( 1, R - 1 )
         IB21E = IB21D + MAX( 1, R )
         IB22D = IB21E + MAX( 1, R - 1 )
         IB22E = IB22D + MAX( 1, R )
         IBBCSD = IB22E + MAX( 1, R - 1 )
         ITAUP1 = IPHI + MAX( 1, R-1 )
         ITAUP2 = ITAUP1 + MAX( 1, P )
         ITAUQ1 = ITAUP2 + MAX( 1, M-P )
         IORBDB = ITAUQ1 + MAX( 1, Q )
         IORGQR = ITAUQ1 + MAX( 1, Q )
         IORGLQ = ITAUQ1 + MAX( 1, Q )
         LORGQRMIN = 1
         LORGQROPT = 1
         LORGLQMIN = 1
         LORGLQOPT = 1
         if ( R .EQ. Q ) {
            sorbdb1(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM1, DUM1, DUM1, DUM1, WORK, -1, CHILDINFO );
            LORBDB = INT( WORK(1) )
            if ( WANTU1 .AND. P .GT. 0 ) {
               sorgqr(P, P, Q, U1, LDU1, DUM1, WORK(1), -1, CHILDINFO );
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            ENDIF
            if ( WANTU2 .AND. M-P .GT. 0 ) {
               sorgqr(M-P, M-P, Q, U2, LDU2, DUM1, WORK(1), -1, CHILDINFO );
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            }
            if ( WANTV1T .AND. Q .GT. 0 ) {
               sorglq(Q-1, Q-1, Q-1, V1T, LDV1T, DUM1, WORK(1), -1, CHILDINFO );
               LORGLQMIN = MAX( LORGLQMIN, Q-1 )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            }
            sbbcsd(JOBU1, JOBU2, JOBV1T, 'N', 'N', M, P, Q, THETA, DUM1, U1, LDU1, U2, LDU2, V1T, LDV1T, DUM2, 1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO );
            LBBCSD = INT( WORK(1) )
         } else if ( R .EQ. P ) {
            sorbdb2(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM1, DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO );
            LORBDB = INT( WORK(1) )
            if ( WANTU1 .AND. P .GT. 0 ) {
               sorgqr(P-1, P-1, P-1, U1(2,2), LDU1, DUM1, WORK(1), -1, CHILDINFO );
               LORGQRMIN = MAX( LORGQRMIN, P-1 )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            }
            if ( WANTU2 .AND. M-P .GT. 0 ) {
               sorgqr(M-P, M-P, Q, U2, LDU2, DUM1, WORK(1), -1, CHILDINFO );
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            }
            if ( WANTV1T .AND. Q .GT. 0 ) {
               sorglq(Q, Q, R, V1T, LDV1T, DUM1, WORK(1), -1, CHILDINFO );
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            }
            sbbcsd(JOBV1T, 'N', JOBU1, JOBU2, 'T', M, Q, P, THETA, DUM1, V1T, LDV1T, DUM2, 1, U1, LDU1, U2, LDU2, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO );
            LBBCSD = INT( WORK(1) )
         } else if ( R .EQ. M-P ) {
            sorbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM1, DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO );
            LORBDB = INT( WORK(1) )
            if ( WANTU1 .AND. P .GT. 0 ) {
               sorgqr(P, P, Q, U1, LDU1, DUM1, WORK(1), -1, CHILDINFO );
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            }
            if ( WANTU2 .AND. M-P .GT. 0 ) {
               sorgqr(M-P-1, M-P-1, M-P-1, U2(2,2), LDU2, DUM1, WORK(1), -1, CHILDINFO );
               LORGQRMIN = MAX( LORGQRMIN, M-P-1 )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            }
            if ( WANTV1T .AND. Q .GT. 0 ) {
               sorglq(Q, Q, R, V1T, LDV1T, DUM1, WORK(1), -1, CHILDINFO );
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            }
            sbbcsd('N', JOBV1T, JOBU2, JOBU1, 'T', M, M-Q, M-P, THETA, DUM1, DUM2, 1, V1T, LDV1T, U2, LDU2, U1, LDU1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO );
            LBBCSD = INT( WORK(1) )
         } else {
            sorbdb4(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM1, DUM1, DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO );
            LORBDB = M + INT( WORK(1) )
            if ( WANTU1 .AND. P .GT. 0 ) {
               sorgqr(P, P, M-Q, U1, LDU1, DUM1, WORK(1), -1, CHILDINFO );
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            }
            if ( WANTU2 .AND. M-P .GT. 0 ) {
               sorgqr(M-P, M-P, M-Q, U2, LDU2, DUM1, WORK(1), -1, CHILDINFO );
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            }
            if ( WANTV1T .AND. Q .GT. 0 ) {
               sorglq(Q, Q, Q, V1T, LDV1T, DUM1, WORK(1), -1, CHILDINFO );
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            }
            sbbcsd(JOBU2, JOBU1, 'N', JOBV1T, 'N', M, M-P, M-Q, THETA, DUM1, U2, LDU2, U1, LDU1, DUM2, 1, V1T, LDV1T, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, DUM1, WORK(1), -1, CHILDINFO );
            LBBCSD = INT( WORK(1) )
         }
         LWORKMIN = MAX( IORBDB+LORBDB-1, IORGQR+LORGQRMIN-1, IORGLQ+LORGLQMIN-1, IBBCSD+LBBCSD-1 )          LWORKOPT = MAX( IORBDB+LORBDB-1, IORGQR+LORGQROPT-1, IORGLQ+LORGLQOPT-1, IBBCSD+LBBCSD-1 )
         WORK(1) = LWORKOPT
         if ( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) {
            INFO = -19
         }
      }
      if ( INFO .NE. 0 ) {
         xerbla('SORCSD2BY1', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }
      LORGQR = LWORK-IORGQR+1
      LORGLQ = LWORK-IORGLQ+1

      // Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q,
      // in which R = MIN(P,M-P,Q,M-Q)

      if ( R .EQ. Q ) {

         // Case 1: R = Q

         // Simultaneously bidiagonalize X11 and X21

         sorbdb1(M, P, Q, X11, LDX11, X21, LDX21, THETA, WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO );

         // Accumulate Householder reflectors

         if ( WANTU1 .AND. P .GT. 0 ) {
            slacpy('L', P, Q, X11, LDX11, U1, LDU1 );
            sorgqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTU2 .AND. M-P .GT. 0 ) {
            slacpy('L', M-P, Q, X21, LDX21, U2, LDU2 );
            sorgqr(M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTV1T .AND. Q .GT. 0 ) {
            V1T(1,1) = ONE
            DO J = 2, Q
               V1T(1,J) = ZERO
               V1T(J,1) = ZERO
            END DO
            slacpy('U', Q-1, Q-1, X21(1,2), LDX21, V1T(2,2), LDV1T )             CALL SORGLQ( Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO );
         }

         // Simultaneously diagonalize X11 and X21.

         sbbcsd(JOBU1, JOBU2, JOBV1T, 'N', 'N', M, P, Q, THETA, WORK(IPHI), U1, LDU1, U2, LDU2, V1T, LDV1T, DUM2, 1, WORK(IB11D), WORK(IB11E), WORK(IB12D), WORK(IB12E), WORK(IB21D), WORK(IB21E), WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSD, CHILDINFO );

         // Permute rows and columns to place zero submatrices in
         // preferred positions

         if ( Q .GT. 0 .AND. WANTU2 ) {
            DO I = 1, Q
               IWORK(I) = M - P - Q + I
            END DO
            DO I = Q + 1, M - P
               IWORK(I) = I - Q
            END DO
            slapmt(.FALSE., M-P, M-P, U2, LDU2, IWORK );
         }
      } else if ( R .EQ. P ) {

         // Case 2: R = P

         // Simultaneously bidiagonalize X11 and X21

         sorbdb2(M, P, Q, X11, LDX11, X21, LDX21, THETA, WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO );

         // Accumulate Householder reflectors

         if ( WANTU1 .AND. P .GT. 0 ) {
            U1(1,1) = ONE
            DO J = 2, P
               U1(1,J) = ZERO
               U1(J,1) = ZERO
            END DO
            slacpy('L', P-1, P-1, X11(2,1), LDX11, U1(2,2), LDU1 );
            sorgqr(P-1, P-1, P-1, U1(2,2), LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTU2 .AND. M-P .GT. 0 ) {
            slacpy('L', M-P, Q, X21, LDX21, U2, LDU2 );
            sorgqr(M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTV1T .AND. Q .GT. 0 ) {
            slacpy('U', P, Q, X11, LDX11, V1T, LDV1T );
            sorglq(Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO );
         }

         // Simultaneously diagonalize X11 and X21.

         sbbcsd(JOBV1T, 'N', JOBU1, JOBU2, 'T', M, Q, P, THETA, WORK(IPHI), V1T, LDV1T, DUM1, 1, U1, LDU1, U2, LDU2, WORK(IB11D), WORK(IB11E), WORK(IB12D), WORK(IB12E), WORK(IB21D), WORK(IB21E), WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSD, CHILDINFO );

         // Permute rows and columns to place identity submatrices in
         // preferred positions

         if ( Q .GT. 0 .AND. WANTU2 ) {
            DO I = 1, Q
               IWORK(I) = M - P - Q + I
            END DO
            DO I = Q + 1, M - P
               IWORK(I) = I - Q
            END DO
            slapmt(.FALSE., M-P, M-P, U2, LDU2, IWORK );
         }
      } else if ( R .EQ. M-P ) {

         // Case 3: R = M-P

         // Simultaneously bidiagonalize X11 and X21

         sorbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO );

         // Accumulate Householder reflectors

         if ( WANTU1 .AND. P .GT. 0 ) {
            slacpy('L', P, Q, X11, LDX11, U1, LDU1 );
            sorgqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTU2 .AND. M-P .GT. 0 ) {
            U2(1,1) = ONE
            DO J = 2, M-P
               U2(1,J) = ZERO
               U2(J,1) = ZERO
            END DO
            slacpy('L', M-P-1, M-P-1, X21(2,1), LDX21, U2(2,2), LDU2 )             CALL SORGQR( M-P-1, M-P-1, M-P-1, U2(2,2), LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTV1T .AND. Q .GT. 0 ) {
            slacpy('U', M-P, Q, X21, LDX21, V1T, LDV1T );
            sorglq(Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO );
         }

         // Simultaneously diagonalize X11 and X21.

         sbbcsd('N', JOBV1T, JOBU2, JOBU1, 'T', M, M-Q, M-P, THETA, WORK(IPHI), DUM1, 1, V1T, LDV1T, U2, LDU2, U1, LDU1, WORK(IB11D), WORK(IB11E), WORK(IB12D), WORK(IB12E), WORK(IB21D), WORK(IB21E), WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSD, CHILDINFO );

         // Permute rows and columns to place identity submatrices in
         // preferred positions

         if ( Q .GT. R ) {
            DO I = 1, R
               IWORK(I) = Q - R + I
            END DO
            DO I = R + 1, Q
               IWORK(I) = I - R
            END DO
            if ( WANTU1 ) {
               slapmt(.FALSE., P, Q, U1, LDU1, IWORK );
            }
            if ( WANTV1T ) {
               slapmr(.FALSE., Q, Q, V1T, LDV1T, IWORK );
            }
         }
      } else {

         // Case 4: R = M-Q

         // Simultaneously bidiagonalize X11 and X21

         sorbdb4(M, P, Q, X11, LDX11, X21, LDX21, THETA, WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), WORK(IORBDB+M), LORBDB-M, CHILDINFO );

         // Accumulate Householder reflectors

         if ( WANTU2 .AND. M-P .GT. 0 ) {
            scopy(M-P, WORK(IORBDB+P), 1, U2, 1 );
         }
         if ( WANTU1 .AND. P .GT. 0 ) {
            scopy(P, WORK(IORBDB), 1, U1, 1 );
            DO J = 2, P
               U1(1,J) = ZERO
            END DO
            slacpy('L', P-1, M-Q-1, X11(2,1), LDX11, U1(2,2), LDU1 )             CALL SORGQR( P, P, M-Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTU2 .AND. M-P .GT. 0 ) {
            DO J = 2, M-P
               U2(1,J) = ZERO
            END DO
            slacpy('L', M-P-1, M-Q-1, X21(2,1), LDX21, U2(2,2), LDU2 )             CALL SORGQR( M-P, M-P, M-Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTV1T .AND. Q .GT. 0 ) {
            slacpy('U', M-Q, Q, X21, LDX21, V1T, LDV1T );
            slacpy('U', P-(M-Q), Q-(M-Q), X11(M-Q+1,M-Q+1), LDX11, V1T(M-Q+1,M-Q+1), LDV1T )             CALL SLACPY( 'U', -P+Q, Q-P, X21(M-Q+1,P+1), LDX21, V1T(P+1,P+1), LDV1T )             CALL SORGLQ( Q, Q, Q, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO );
         }

         // Simultaneously diagonalize X11 and X21.

         sbbcsd(JOBU2, JOBU1, 'N', JOBV1T, 'N', M, M-P, M-Q, THETA, WORK(IPHI), U2, LDU2, U1, LDU1, DUM1, 1, V1T, LDV1T, WORK(IB11D), WORK(IB11E), WORK(IB12D), WORK(IB12E), WORK(IB21D), WORK(IB21E), WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSD, CHILDINFO );

         // Permute rows and columns to place identity submatrices in
         // preferred positions

         if ( P .GT. R ) {
            DO I = 1, R
               IWORK(I) = P - R + I
            END DO
            DO I = R + 1, P
               IWORK(I) = I - R
            END DO
            if ( WANTU1 ) {
               slapmt(.FALSE., P, P, U1, LDU1, IWORK );
            }
            if ( WANTV1T ) {
               slapmr(.FALSE., P, Q, V1T, LDV1T, IWORK );
            }
         }
      }

      RETURN

      // End of SORCSD2BY1

      }
