      RECURSIVE SUBROUTINE SORCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK, LWORK, IWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS;
      int                INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               THETA( * );
      REAL               U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ), X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, * );
      // ..

// ===================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Arrays ..
      REAL               DUMMY(1);
      // ..
      // .. Local Scalars ..
      String             TRANST, SIGNST;
      int                CHILDINFO, I, IB11D, IB11E, IB12D, IB12E, IB21D, IB21E, IB22D, IB22E, IBBCSD, IORBDB, IORGLQ, IORGQR, IPHI, ITAUP1, ITAUP2, ITAUQ1, ITAUQ2, J, LBBCSDWORK, LBBCSDWORKMIN, LBBCSDWORKOPT, LORBDBWORK, LORBDBWORKMIN, LORBDBWORKOPT, LORGLQWORK, LORGLQWORKMIN, LORGLQWORKOPT, LORGQRWORK, LORGQRWORKMIN, LORGQRWORKOPT, LWORKMIN, LWORKOPT;
      bool               COLMAJOR, DEFAULTSIGNS, LQUERY, WANTU1, WANTU2, WANTV1T, WANTV2T;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SBBCSD, SLACPY, SLAPMR, SLAPMT, SORBDB, SORGLQ, SORGQR, XERBLA
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0;
      WANTU1 = LSAME( JOBU1, 'Y' );
      WANTU2 = LSAME( JOBU2, 'Y' );
      WANTV1T = LSAME( JOBV1T, 'Y' );
      WANTV2T = LSAME( JOBV2T, 'Y' );
      COLMAJOR = !LSAME( TRANS, 'T' );
      DEFAULTSIGNS = !LSAME( SIGNS, 'O' );
      LQUERY = LWORK == -1;
      if ( M < 0 ) {
         INFO = -7;
      } else if ( P < 0 || P > M ) {
         INFO = -8;
      } else if ( Q < 0 || Q > M ) {
         INFO = -9;
      } else if ( COLMAJOR && LDX11 < max( 1, P ) ) {
        INFO = -11;
      } else if ( !COLMAJOR && LDX11 < max( 1, Q ) ) {
        INFO = -11;
      } else if (COLMAJOR && LDX12 < max( 1, P ) ) {
        INFO = -13;
      } else if ( !COLMAJOR && LDX12 < max( 1, M-Q ) ) {
        INFO = -13;
      } else if (COLMAJOR && LDX21 < max( 1, M-P ) ) {
        INFO = -15;
      } else if ( !COLMAJOR && LDX21 < max( 1, Q ) ) {
        INFO = -15;
      } else if (COLMAJOR && LDX22 < max( 1, M-P ) ) {
        INFO = -17;
      } else if ( !COLMAJOR && LDX22 < max( 1, M-Q ) ) {
        INFO = -17;
      } else if ( WANTU1 && LDU1 < P ) {
         INFO = -20;
      } else if ( WANTU2 && LDU2 < M-P ) {
         INFO = -22;
      } else if ( WANTV1T && LDV1T < Q ) {
         INFO = -24;
      } else if ( WANTV2T && LDV2T < M-Q ) {
         INFO = -26;
      }

      // Work with transpose if convenient

      if ( INFO == 0 && min( P, M-P ) < min( Q, M-Q ) ) {
         if ( COLMAJOR ) {
            TRANST = 'T';
         } else {
            TRANST = 'N';
         }
         if ( DEFAULTSIGNS ) {
            SIGNST = 'O';
         } else {
            SIGNST = 'D';
         }
         sorcsd(JOBV1T, JOBV2T, JOBU1, JOBU2, TRANST, SIGNST, M, Q, P, X11, LDX11, X21, LDX21, X12, LDX12, X22, LDX22, THETA, V1T, LDV1T, V2T, LDV2T, U1, LDU1, U2, LDU2, WORK, LWORK, IWORK, INFO );
         return;
      }

      // Work with permutation [ 0 I; I 0 ] * X * [ 0 I; I 0 ] if
      // convenient

      if ( INFO == 0 && M-Q < Q ) {
         if ( DEFAULTSIGNS ) {
            SIGNST = 'O';
         } else {
            SIGNST = 'D';
         }
         sorcsd(JOBU2, JOBU1, JOBV2T, JOBV1T, TRANS, SIGNST, M, M-P, M-Q, X22, LDX22, X21, LDX21, X12, LDX12, X11, LDX11, THETA, U2, LDU2, U1, LDU1, V2T, LDV2T, V1T, LDV1T, WORK, LWORK, IWORK, INFO );
         return;
      }

      // Compute workspace

      if ( INFO == 0 ) {

         IPHI = 2;
         ITAUP1 = IPHI + max( 1, Q - 1 );
         ITAUP2 = ITAUP1 + max( 1, P );
         ITAUQ1 = ITAUP2 + max( 1, M - P );
         ITAUQ2 = ITAUQ1 + max( 1, Q );
         IORGQR = ITAUQ2 + max( 1, M - Q );
         sorgqr(M-Q, M-Q, M-Q, DUMMY, max(1,M-Q), DUMMY, WORK, -1, CHILDINFO );
         LORGQRWORKOPT = INT( WORK(1) );
         LORGQRWORKMIN = max( 1, M - Q );
         IORGLQ = ITAUQ2 + max( 1, M - Q );
         sorglq(M-Q, M-Q, M-Q, DUMMY, max(1,M-Q), DUMMY, WORK, -1, CHILDINFO );
         LORGLQWORKOPT = INT( WORK(1) );
         LORGLQWORKMIN = max( 1, M - Q );
         IORBDB = ITAUQ2 + max( 1, M - Q );
         sorbdb(TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY,WORK,-1,CHILDINFO );
         LORBDBWORKOPT = INT( WORK(1) );
         LORBDBWORKMIN = LORBDBWORKOPT;
         IB11D = ITAUQ2 + max( 1, M - Q );
         IB11E = IB11D + max( 1, Q );
         IB12D = IB11E + max( 1, Q - 1 );
         IB12E = IB12D + max( 1, Q );
         IB21D = IB12E + max( 1, Q - 1 );
         IB21E = IB21D + max( 1, Q );
         IB22D = IB21E + max( 1, Q - 1 );
         IB22E = IB22D + max( 1, Q );
         IBBCSD = IB22E + max( 1, Q - 1 );
         sbbcsd(JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, DUMMY, DUMMY, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, WORK, -1, CHILDINFO );
         LBBCSDWORKOPT = INT( WORK(1) );
         LBBCSDWORKMIN = LBBCSDWORKOPT;
         LWORKOPT = max( IORGQR + LORGQRWORKOPT, IORGLQ + LORGLQWORKOPT, IORBDB + LORBDBWORKOPT, IBBCSD + LBBCSDWORKOPT ) - 1          LWORKMIN = max( IORGQR + LORGQRWORKMIN, IORGLQ + LORGLQWORKMIN, IORBDB + LORBDBWORKOPT, IBBCSD + LBBCSDWORKMIN ) - 1;
         WORK(1) = max(LWORKOPT,LWORKMIN);

         if ( LWORK < LWORKMIN && !LQUERY ) {
            INFO = -22;
         } else {
            LORGQRWORK = LWORK - IORGQR + 1;
            LORGLQWORK = LWORK - IORGLQ + 1;
            LORBDBWORK = LWORK - IORBDB + 1;
            LBBCSDWORK = LWORK - IBBCSD + 1;
         }
      }

      // Abort if any illegal arguments

      if ( INFO != 0 ) {
         xerbla('SORCSD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Transform to bidiagonal block form

      sorbdb(TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(ITAUQ2), WORK(IORBDB), LORBDBWORK, CHILDINFO );

      // Accumulate Householder reflectors

      if ( COLMAJOR ) {
         if ( WANTU1 && P > 0 ) {
            slacpy('L', P, Q, X11, LDX11, U1, LDU1 );
            sorgqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQRWORK, INFO);
         }
         if ( WANTU2 && M-P > 0 ) {
            slacpy('L', M-P, Q, X21, LDX21, U2, LDU2 );
            sorgqr(M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQRWORK, INFO );
         }
         if ( WANTV1T && Q > 0 ) {
            slacpy('U', Q-1, Q-1, X11(1,2), LDX11, V1T(2,2), LDV1T );
            V1T(1, 1) = ONE;
            for (J = 2; J <= Q; J++) {
               V1T(1,J) = ZERO;
               V1T(J,1) = ZERO;
            }
            sorglq(Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQWORK, INFO );
         }
         if ( WANTV2T && M-Q > 0 ) {
            slacpy('U', P, M-Q, X12, LDX12, V2T, LDV2T );
            slacpy('U', M-P-Q, M-P-Q, X22(Q+1,P+1), LDX22, V2T(P+1,P+1), LDV2T );
            sorglq(M-Q, M-Q, M-Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGLQ), LORGLQWORK, INFO );
         }
      } else {
         if ( WANTU1 && P > 0 ) {
            slacpy('U', Q, P, X11, LDX11, U1, LDU1 );
            sorglq(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGLQ), LORGLQWORK, INFO);
         }
         if ( WANTU2 && M-P > 0 ) {
            slacpy('U', Q, M-P, X21, LDX21, U2, LDU2 );
            sorglq(M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGLQ), LORGLQWORK, INFO );
         }
         if ( WANTV1T && Q > 0 ) {
            slacpy('L', Q-1, Q-1, X11(2,1), LDX11, V1T(2,2), LDV1T );
            V1T(1, 1) = ONE;
            for (J = 2; J <= Q; J++) {
               V1T(1,J) = ZERO;
               V1T(J,1) = ZERO;
            }
            sorgqr(Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGQR), LORGQRWORK, INFO );
         }
         if ( WANTV2T && M-Q > 0 ) {
            slacpy('L', M-Q, P, X12, LDX12, V2T, LDV2T );
            slacpy('L', M-P-Q, M-P-Q, X22(P+1,Q+1), LDX22, V2T(P+1,P+1), LDV2T );
            sorgqr(M-Q, M-Q, M-Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGQR), LORGQRWORK, INFO );
         }
      }

      // Compute the CSD of the matrix in bidiagonal-block form

      sbbcsd(JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, THETA, WORK(IPHI), U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK(IB11D), WORK(IB11E), WORK(IB12D), WORK(IB12E), WORK(IB21D), WORK(IB21E), WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSDWORK, INFO );

      // Permute rows and columns to place identity submatrices in top-
      // left corner of (1,1)-block and/or bottom-right corner of (1,2)-
      // block and/or bottom-right corner of (2,1)-block and/or top-left
      // corner of (2,2)-block

      if ( Q > 0 && WANTU2 ) {
         for (I = 1; I <= Q; I++) {
            IWORK(I) = M - P - Q + I;
         }
         for (I = Q + 1; I <= M - P; I++) {
            IWORK(I) = I - Q;
         }
         if ( COLMAJOR ) {
            slapmt( false , M-P, M-P, U2, LDU2, IWORK );
         } else {
            slapmr( false , M-P, M-P, U2, LDU2, IWORK );
         }
      }
      if ( M > 0 && WANTV2T ) {
         for (I = 1; I <= P; I++) {
            IWORK(I) = M - P - Q + I;
         }
         for (I = P + 1; I <= M - Q; I++) {
            IWORK(I) = I - P;
         }
         if ( !COLMAJOR ) {
            slapmt( false , M-Q, M-Q, V2T, LDV2T, IWORK );
         } else {
            slapmr( false , M-Q, M-Q, V2T, LDV2T, IWORK );
         }
      }

      return;

      // End SORCSD

      }
