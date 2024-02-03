      RECURSIVE SUBROUTINE CUNCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK, LWORK, RWORK, LRWORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS;
      int                INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12, LDX21, LDX22, LRWORK, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               THETA( * )
      REAL               RWORK( * )
      COMPLEX            U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ), X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, * )
      // ..

*  ===================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO
      const              ONE = (1.0E0,0.0E0), ZERO = (0.0E0,0.0E0) ;
      // ..
      // .. Local Scalars ..
      String             TRANST, SIGNST;
      int                CHILDINFO, I, IB11D, IB11E, IB12D, IB12E, IB21D, IB21E, IB22D, IB22E, IBBCSD, IORBDB, IORGLQ, IORGQR, IPHI, ITAUP1, ITAUP2, ITAUQ1, ITAUQ2, J, LBBCSDWORK, LBBCSDWORKMIN, LBBCSDWORKOPT, LORBDBWORK, LORBDBWORKMIN, LORBDBWORKOPT, LORGLQWORK, LORGLQWORKMIN, LORGLQWORKOPT, LORGQRWORK, LORGQRWORKMIN, LORGQRWORKOPT, LWORKMIN, LWORKOPT, P1, Q1;
      bool               COLMAJOR, DEFAULTSIGNS, LQUERY, WANTU1, WANTU2, WANTV1T, WANTV2T;
      int                LRWORKMIN, LRWORKOPT;
      bool               LRQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CBBCSD, CLACPY, CLAPMR, CLAPMT, CUNBDB, CUNGLQ, CUNGQR
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTV1T = LSAME( JOBV1T, 'Y' )
      WANTV2T = LSAME( JOBV2T, 'Y' )
      COLMAJOR = .NOT. LSAME( TRANS, 'T' )
      DEFAULTSIGNS = .NOT. LSAME( SIGNS, 'O' )
      LQUERY = LWORK == -1
      LRQUERY = LRWORK == -1
      if ( M .LT. 0 ) {
         INFO = -7
      } else if ( P .LT. 0 .OR. P .GT. M ) {
         INFO = -8
      } else if ( Q .LT. 0 .OR. Q .GT. M ) {
         INFO = -9
      } else if ( COLMAJOR && LDX11 .LT. MAX( 1, P ) ) {
        INFO = -11
      } else if (.NOT. COLMAJOR && LDX11 .LT. MAX( 1, Q ) ) {
        INFO = -11
      } else if (COLMAJOR && LDX12 .LT. MAX( 1, P ) ) {
        INFO = -13
      } else if (.NOT. COLMAJOR && LDX12 .LT. MAX( 1, M-Q ) ) {
        INFO = -13
      } else if (COLMAJOR && LDX21 .LT. MAX( 1, M-P ) ) {
        INFO = -15
      } else if (.NOT. COLMAJOR && LDX21 .LT. MAX( 1, Q ) ) {
        INFO = -15
      } else if (COLMAJOR && LDX22 .LT. MAX( 1, M-P ) ) {
        INFO = -17
      } else if (.NOT. COLMAJOR && LDX22 .LT. MAX( 1, M-Q ) ) {
        INFO = -17
      } else if ( WANTU1 && LDU1 .LT. P ) {
         INFO = -20
      } else if ( WANTU2 && LDU2 .LT. M-P ) {
         INFO = -22
      } else if ( WANTV1T && LDV1T .LT. Q ) {
         INFO = -24
      } else if ( WANTV2T && LDV2T .LT. M-Q ) {
         INFO = -26
      }

      // Work with transpose if convenient

      if ( INFO == 0 && MIN( P, M-P ) .LT. MIN( Q, M-Q ) ) {
         if ( COLMAJOR ) {
            TRANST = 'T'
         } else {
            TRANST = 'N'
         }
         if ( DEFAULTSIGNS ) {
            SIGNST = 'O'
         } else {
            SIGNST = 'D'
         }
         cuncsd(JOBV1T, JOBV2T, JOBU1, JOBU2, TRANST, SIGNST, M, Q, P, X11, LDX11, X21, LDX21, X12, LDX12, X22, LDX22, THETA, V1T, LDV1T, V2T, LDV2T, U1, LDU1, U2, LDU2, WORK, LWORK, RWORK, LRWORK, IWORK, INFO );
         RETURN
      }

      // Work with permutation [ 0 I; I 0 ] * X * [ 0 I; I 0 ] if
      // convenient

      if ( INFO == 0 && M-Q .LT. Q ) {
         if ( DEFAULTSIGNS ) {
            SIGNST = 'O'
         } else {
            SIGNST = 'D'
         }
         cuncsd(JOBU2, JOBU1, JOBV2T, JOBV1T, TRANS, SIGNST, M, M-P, M-Q, X22, LDX22, X21, LDX21, X12, LDX12, X11, LDX11, THETA, U2, LDU2, U1, LDU1, V2T, LDV2T, V1T, LDV1T, WORK, LWORK, RWORK, LRWORK, IWORK, INFO );
         RETURN
      }

      // Compute workspace

      if ( INFO == 0 ) {

         // Real workspace

         IPHI = 2
         IB11D = IPHI + MAX( 1, Q - 1 )
         IB11E = IB11D + MAX( 1, Q )
         IB12D = IB11E + MAX( 1, Q - 1 )
         IB12E = IB12D + MAX( 1, Q )
         IB21D = IB12E + MAX( 1, Q - 1 )
         IB21E = IB21D + MAX( 1, Q )
         IB22D = IB21E + MAX( 1, Q - 1 )
         IB22E = IB22D + MAX( 1, Q )
         IBBCSD = IB22E + MAX( 1, Q - 1 )
         cbbcsd(JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, THETA, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, THETA, THETA, THETA, THETA, THETA, THETA, THETA, THETA, RWORK, -1, CHILDINFO );
         LBBCSDWORKOPT = INT( RWORK(1) )
         LBBCSDWORKMIN = LBBCSDWORKOPT
         LRWORKOPT = IBBCSD + LBBCSDWORKOPT - 1
         LRWORKMIN = IBBCSD + LBBCSDWORKMIN - 1
         RWORK(1) = LRWORKOPT

         // Complex workspace

         ITAUP1 = 2
         ITAUP2 = ITAUP1 + MAX( 1, P )
         ITAUQ1 = ITAUP2 + MAX( 1, M - P )
         ITAUQ2 = ITAUQ1 + MAX( 1, Q )
         IORGQR = ITAUQ2 + MAX( 1, M - Q )
         cungqr(M-Q, M-Q, M-Q, U1, MAX(1,M-Q), U1, WORK, -1, CHILDINFO );
         LORGQRWORKOPT = INT( WORK(1) )
         LORGQRWORKMIN = MAX( 1, M - Q )
         IORGLQ = ITAUQ2 + MAX( 1, M - Q )
         cunglq(M-Q, M-Q, M-Q, U1, MAX(1,M-Q), U1, WORK, -1, CHILDINFO );
         LORGLQWORKOPT = INT( WORK(1) )
         LORGLQWORKMIN = MAX( 1, M - Q )
         IORBDB = ITAUQ2 + MAX( 1, M - Q )
         cunbdb(TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, THETA, U1, U2, V1T, V2T, WORK, -1, CHILDINFO );
         LORBDBWORKOPT = INT( WORK(1) )
         LORBDBWORKMIN = LORBDBWORKOPT
         LWORKOPT = MAX( IORGQR + LORGQRWORKOPT, IORGLQ + LORGLQWORKOPT, IORBDB + LORBDBWORKOPT ) - 1          LWORKMIN = MAX( IORGQR + LORGQRWORKMIN, IORGLQ + LORGLQWORKMIN, IORBDB + LORBDBWORKMIN ) - 1
         LWORKOPT = MAX(LWORKOPT,LWORKMIN)
         WORK(1) = SROUNDUP_LWORK(LWORKOPT)

         if ( LWORK .LT. LWORKMIN && .NOT. ( LQUERY .OR. LRQUERY ) ) {
            INFO = -22
         } else if ( LRWORK .LT. LRWORKMIN && .NOT. ( LQUERY .OR. LRQUERY ) ) {
            INFO = -24
         } else {
            LORGQRWORK = LWORK - IORGQR + 1
            LORGLQWORK = LWORK - IORGLQ + 1
            LORBDBWORK = LWORK - IORBDB + 1
            LBBCSDWORK = LRWORK - IBBCSD + 1
         }
      }

      // Abort if any illegal arguments

      if ( INFO != 0 ) {
         xerbla('CUNCSD', -INFO );
         RETURN
      } else if ( LQUERY .OR. LRQUERY ) {
         RETURN
      }

      // Transform to bidiagonal block form

      cunbdb(TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(ITAUQ2), WORK(IORBDB), LORBDBWORK, CHILDINFO );

      // Accumulate Householder reflectors

      if ( COLMAJOR ) {
         if ( WANTU1 && P .GT. 0 ) {
            clacpy('L', P, Q, X11, LDX11, U1, LDU1 );
            cungqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQRWORK, INFO);
         }
         if ( WANTU2 && M-P .GT. 0 ) {
            clacpy('L', M-P, Q, X21, LDX21, U2, LDU2 );
            cungqr(M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQRWORK, INFO );
         }
         if ( WANTV1T && Q .GT. 0 ) {
            clacpy('U', Q-1, Q-1, X11(1,2), LDX11, V1T(2,2), LDV1T );
            V1T(1, 1) = ONE
            for (J = 2; J <= Q; J++) {
               V1T(1,J) = ZERO
               V1T(J,1) = ZERO
            }
            cunglq(Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQWORK, INFO );
         }
         if ( WANTV2T && M-Q .GT. 0 ) {
            clacpy('U', P, M-Q, X12, LDX12, V2T, LDV2T );
            if ( M-P .GT. Q ) {
               clacpy('U', M-P-Q, M-P-Q, X22(Q+1,P+1), LDX22, V2T(P+1,P+1), LDV2T );
            }
            if ( M .GT. Q ) {
               cunglq(M-Q, M-Q, M-Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGLQ), LORGLQWORK, INFO );
            }
         }
      } else {
         if ( WANTU1 && P .GT. 0 ) {
            clacpy('U', Q, P, X11, LDX11, U1, LDU1 );
            cunglq(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGLQ), LORGLQWORK, INFO);
         }
         if ( WANTU2 && M-P .GT. 0 ) {
            clacpy('U', Q, M-P, X21, LDX21, U2, LDU2 );
            cunglq(M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGLQ), LORGLQWORK, INFO );
         }
         if ( WANTV1T && Q .GT. 0 ) {
            clacpy('L', Q-1, Q-1, X11(2,1), LDX11, V1T(2,2), LDV1T );
            V1T(1, 1) = ONE
            for (J = 2; J <= Q; J++) {
               V1T(1,J) = ZERO
               V1T(J,1) = ZERO
            }
            cungqr(Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGQR), LORGQRWORK, INFO );
         }
         if ( WANTV2T && M-Q .GT. 0 ) {
            P1 = MIN( P+1, M )
            Q1 = MIN( Q+1, M )
            clacpy('L', M-Q, P, X12, LDX12, V2T, LDV2T );
            if ( M .GT. P+Q ) {
               clacpy('L', M-P-Q, M-P-Q, X22(P1,Q1), LDX22, V2T(P+1,P+1), LDV2T );
            }
            cungqr(M-Q, M-Q, M-Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGQR), LORGQRWORK, INFO );
         }
      }

      // Compute the CSD of the matrix in bidiagonal-block form

      cbbcsd(JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, THETA, RWORK(IPHI), U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LBBCSDWORK, INFO );

      // Permute rows and columns to place identity submatrices in top-
      // left corner of (1,1)-block and/or bottom-right corner of (1,2)-
      // block and/or bottom-right corner of (2,1)-block and/or top-left
      // corner of (2,2)-block

      if ( Q .GT. 0 && WANTU2 ) {
         for (I = 1; I <= Q; I++) {
            IWORK(I) = M - P - Q + I
         }
         for (I = Q + 1; I <= M - P; I++) {
            IWORK(I) = I - Q
         }
         if ( COLMAJOR ) {
            clapmt( false , M-P, M-P, U2, LDU2, IWORK );
         } else {
            clapmr( false , M-P, M-P, U2, LDU2, IWORK );
         }
      }
      if ( M .GT. 0 && WANTV2T ) {
         for (I = 1; I <= P; I++) {
            IWORK(I) = M - P - Q + I
         }
         for (I = P + 1; I <= M - Q; I++) {
            IWORK(I) = I - P
         }
         if ( .NOT. COLMAJOR ) {
            clapmt( false , M-Q, M-Q, V2T, LDV2T, IWORK );
         } else {
            clapmr( false , M-Q, M-Q, V2T, LDV2T, IWORK );
         }
      }

      RETURN

      // End CUNCSD

      }
