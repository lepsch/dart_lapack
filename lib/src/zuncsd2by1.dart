      void zuncsd2by1(final int JOBU1, final int JOBU2, final int JOBV1T, final int M, final int P, final int Q, final Matrix<double> X11, final int LDX11, final Matrix<double> X21, final int LDX21, final int THETA, final Matrix<double> U1, final int LDU1, final Matrix<double> U2, final int LDU2, final Matrix<double> V1T, final int LDV1T, final Array<double> WORK, final int LWORK, final Array<int> RWORK, final int LRWORK, final Array<int> IWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBU1, JOBU2, JOBV1T;
      int                INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21, M, P, Q;
      int                LRWORK, LRWORKMIN, LRWORKOPT;
      double             RWORK(*);
      double             THETA(*);
      Complex         U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*), X11(LDX11,*), X21(LDX21,*);
      int                IWORK(*);
      // ..

      Complex         ONE, ZERO;
      const              ONE = (1.0,0.0), ZERO = (0.0,0.0) ;
      int                CHILDINFO, I, IB11D, IB11E, IB12D, IB12E, IB21D, IB21E, IB22D, IB22E, IBBCSD, IORBDB, IORGLQ, IORGQR, IPHI, ITAUP1, ITAUP2, ITAUQ1, J, LBBCSD, LORBDB, LORGLQ, LORGLQMIN, LORGLQOPT, LORGQR, LORGQRMIN, LORGQROPT, LWORKMIN, LWORKOPT, R;
      bool               LQUERY, WANTU1, WANTU2, WANTV1T;
      double             DUM( 1 );
      Complex         CDUM( 1, 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZBBCSD, ZCOPY, ZLACPY, ZLAPMR, ZLAPMT, ZUNBDB1, ZUNBDB2, ZUNBDB3, ZUNBDB4, ZUNGLQ, ZUNGQR, XERBLA
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC INT, MAX, MIN

      // Test input arguments

      INFO = 0;
      WANTU1 = lsame( JOBU1, 'Y' );
      WANTU2 = lsame( JOBU2, 'Y' );
      WANTV1T = lsame( JOBV1T, 'Y' );
      LQUERY = ( LWORK == -1 ) || ( LRWORK == -1 );

      if ( M < 0 ) {
         INFO = -4;
      } else if ( P < 0 || P > M ) {
         INFO = -5;
      } else if ( Q < 0 || Q > M ) {
         INFO = -6;
      } else if ( LDX11 < max( 1, P ) ) {
         INFO = -8;
      } else if ( LDX21 < max( 1, M-P ) ) {
         INFO = -10;
      } else if ( WANTU1 && LDU1 < max( 1, P ) ) {
         INFO = -13;
      } else if ( WANTU2 && LDU2 < max( 1, M - P ) ) {
         INFO = -15;
      } else if ( WANTV1T && LDV1T < max( 1, Q ) ) {
         INFO = -17;
      }

      R = min( P, M-P, Q, M-Q );

      // Compute workspace

        // WORK layout:
      // |-----------------------------------------|
      // | LWORKOPT (1)                            |
      // |-----------------------------------------|
      // | TAUP1 (max(1,P))                        |
      // | TAUP2 (max(1,M-P))                      |
      // | TAUQ1 (max(1,Q))                        |
      // |-----------------------------------------|
      // | ZUNBDB WORK | ZUNGQR WORK | ZUNGLQ WORK |
      // |             |             |             |
      // |             |             |             |
      // |             |             |             |
      // |             |             |             |
      // |-----------------------------------------|
      //   RWORK layout:
      // |------------------|
      // | LRWORKOPT (1)    |
      // |------------------|
      // | PHI (max(1,R-1)) |
      // |------------------|
      // | B11D (R)         |
      // | B11E (R-1)       |
      // | B12D (R)         |
      // | B12E (R-1)       |
      // | B21D (R)         |
      // | B21E (R-1)       |
      // | B22D (R)         |
      // | B22E (R-1)       |
      // | ZBBCSD RWORK     |
      // |------------------|

      if ( INFO == 0 ) {
         IPHI = 2;
         IB11D = IPHI + max( 1, R-1 );
         IB11E = IB11D + max( 1, R );
         IB12D = IB11E + max( 1, R - 1 );
         IB12E = IB12D + max( 1, R );
         IB21D = IB12E + max( 1, R - 1 );
         IB21E = IB21D + max( 1, R );
         IB22D = IB21E + max( 1, R - 1 );
         IB22E = IB22D + max( 1, R );
         IBBCSD = IB22E + max( 1, R - 1 );
         ITAUP1 = 2;
         ITAUP2 = ITAUP1 + max( 1, P );
         ITAUQ1 = ITAUP2 + max( 1, M-P );
         IORBDB = ITAUQ1 + max( 1, Q );
         IORGQR = ITAUQ1 + max( 1, Q );
         IORGLQ = ITAUQ1 + max( 1, Q );
         LORGQRMIN = 1;
         LORGQROPT = 1;
         LORGLQMIN = 1;
         LORGLQOPT = 1;
         if ( R == Q ) {
            zunbdb1(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM, CDUM, CDUM, WORK, -1, CHILDINFO );
            LORBDB = INT( WORK(1) );
            if ( WANTU1 && P > 0 ) {
               zungqr(P, P, Q, U1, LDU1, CDUM, WORK(1), -1, CHILDINFO );
               LORGQRMIN = max( LORGQRMIN, P );
               LORGQROPT = max( LORGQROPT, INT( WORK(1) ) );
            }
            if ( WANTU2 && M-P > 0 ) {
               zungqr(M-P, M-P, Q, U2, LDU2, CDUM, WORK(1), -1, CHILDINFO );
               LORGQRMIN = max( LORGQRMIN, M-P );
               LORGQROPT = max( LORGQROPT, INT( WORK(1) ) );
            }
            if ( WANTV1T && Q > 0 ) {
               zunglq(Q-1, Q-1, Q-1, V1T, LDV1T, CDUM, WORK(1), -1, CHILDINFO );
               LORGLQMIN = max( LORGLQMIN, Q-1 );
               LORGLQOPT = max( LORGLQOPT, INT( WORK(1) ) );
            }
            zbbcsd(JOBU1, JOBU2, JOBV1T, 'N', 'N', M, P, Q, THETA, DUM, U1, LDU1, U2, LDU2, V1T, LDV1T, CDUM, 1, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, RWORK(1), -1, CHILDINFO );
            LBBCSD = INT( RWORK(1) );
         } else if ( R == P ) {
            zunbdb2(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM, CDUM, CDUM, WORK(1), -1, CHILDINFO );
            LORBDB = INT( WORK(1) );
            if ( WANTU1 && P > 0 ) {
               zungqr(P-1, P-1, P-1, U1(2,2), LDU1, CDUM, WORK(1), -1, CHILDINFO );
               LORGQRMIN = max( LORGQRMIN, P-1 );
               LORGQROPT = max( LORGQROPT, INT( WORK(1) ) );
            }
            if ( WANTU2 && M-P > 0 ) {
               zungqr(M-P, M-P, Q, U2, LDU2, CDUM, WORK(1), -1, CHILDINFO );
               LORGQRMIN = max( LORGQRMIN, M-P );
               LORGQROPT = max( LORGQROPT, INT( WORK(1) ) );
            }
            if ( WANTV1T && Q > 0 ) {
               zunglq(Q, Q, R, V1T, LDV1T, CDUM, WORK(1), -1, CHILDINFO );
               LORGLQMIN = max( LORGLQMIN, Q );
               LORGLQOPT = max( LORGLQOPT, INT( WORK(1) ) );
            }
            zbbcsd(JOBV1T, 'N', JOBU1, JOBU2, 'T', M, Q, P, THETA, DUM, V1T, LDV1T, CDUM, 1, U1, LDU1, U2, LDU2, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, RWORK(1), -1, CHILDINFO );
            LBBCSD = INT( RWORK(1) );
         } else if ( R == M-P ) {
            zunbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM, CDUM, CDUM, WORK(1), -1, CHILDINFO );
            LORBDB = INT( WORK(1) );
            if ( WANTU1 && P > 0 ) {
               zungqr(P, P, Q, U1, LDU1, CDUM, WORK(1), -1, CHILDINFO );
               LORGQRMIN = max( LORGQRMIN, P );
               LORGQROPT = max( LORGQROPT, INT( WORK(1) ) );
            }
            if ( WANTU2 && M-P > 0 ) {
               zungqr(M-P-1, M-P-1, M-P-1, U2(2,2), LDU2, CDUM, WORK(1), -1, CHILDINFO );
               LORGQRMIN = max( LORGQRMIN, M-P-1 );
               LORGQROPT = max( LORGQROPT, INT( WORK(1) ) );
            }
            if ( WANTV1T && Q > 0 ) {
               zunglq(Q, Q, R, V1T, LDV1T, CDUM, WORK(1), -1, CHILDINFO );
               LORGLQMIN = max( LORGLQMIN, Q );
               LORGLQOPT = max( LORGLQOPT, INT( WORK(1) ) );
            }
            zbbcsd('N', JOBV1T, JOBU2, JOBU1, 'T', M, M-Q, M-P, THETA, DUM, CDUM, 1, V1T, LDV1T, U2, LDU2, U1, LDU1, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, RWORK(1), -1, CHILDINFO );
            LBBCSD = INT( RWORK(1) );
         } else {
            zunbdb4(M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM, CDUM, CDUM, CDUM, WORK(1), -1, CHILDINFO );
            LORBDB = M + INT( WORK(1) );
            if ( WANTU1 && P > 0 ) {
               zungqr(P, P, M-Q, U1, LDU1, CDUM, WORK(1), -1, CHILDINFO );
               LORGQRMIN = max( LORGQRMIN, P );
               LORGQROPT = max( LORGQROPT, INT( WORK(1) ) );
            }
            if ( WANTU2 && M-P > 0 ) {
               zungqr(M-P, M-P, M-Q, U2, LDU2, CDUM, WORK(1), -1, CHILDINFO );
               LORGQRMIN = max( LORGQRMIN, M-P );
               LORGQROPT = max( LORGQROPT, INT( WORK(1) ) );
            }
            if ( WANTV1T && Q > 0 ) {
               zunglq(Q, Q, Q, V1T, LDV1T, CDUM, WORK(1), -1, CHILDINFO );
               LORGLQMIN = max( LORGLQMIN, Q );
               LORGLQOPT = max( LORGLQOPT, INT( WORK(1) ) );
            }
            zbbcsd(JOBU2, JOBU1, 'N', JOBV1T, 'N', M, M-P, M-Q, THETA, DUM, U2, LDU2, U1, LDU1, CDUM, 1, V1T, LDV1T, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, RWORK(1), -1, CHILDINFO );
            LBBCSD = INT( RWORK(1) );
         }
         LRWORKMIN = IBBCSD+LBBCSD-1;
         LRWORKOPT = LRWORKMIN;
         RWORK[1] = LRWORKOPT;
         LWORKMIN = max( IORBDB+LORBDB-1, IORGQR+LORGQRMIN-1, IORGLQ+LORGLQMIN-1 )          LWORKOPT = max( IORBDB+LORBDB-1, IORGQR+LORGQROPT-1, IORGLQ+LORGLQOPT-1 );
         WORK[1] = LWORKOPT;
         if ( LWORK < LWORKMIN && !LQUERY ) {
            INFO = -19;
         }
         if ( LRWORK < LRWORKMIN && !LQUERY ) {
            INFO = -21;
         }
      }
      if ( INFO != 0 ) {
         xerbla('ZUNCSD2BY1', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }
      LORGQR = LWORK-IORGQR+1;
      LORGLQ = LWORK-IORGLQ+1;

      // Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q,
      // in which R = min(P,M-P,Q,M-Q)

      if ( R == Q ) {

         // Case 1: R = Q

         // Simultaneously bidiagonalize X11 and X21

         zunbdb1(M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO );

         // Accumulate Householder reflectors

         if ( WANTU1 && P > 0 ) {
            zlacpy('L', P, Q, X11, LDX11, U1, LDU1 );
            zungqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTU2 && M-P > 0 ) {
            zlacpy('L', M-P, Q, X21, LDX21, U2, LDU2 );
            zungqr(M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTV1T && Q > 0 ) {
            V1T[1][1] = ONE;
            for (J = 2; J <= Q; J++) {
               V1T[1][J] = ZERO;
               V1T[J][1] = ZERO;
            }
            zlacpy('U', Q-1, Q-1, X21(1,2), LDX21, V1T(2,2), LDV1T );
            zunglq(Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO );
         }

         // Simultaneously diagonalize X11 and X21.

         zbbcsd(JOBU1, JOBU2, JOBV1T, 'N', 'N', M, P, Q, THETA, RWORK(IPHI), U1, LDU1, U2, LDU2, V1T, LDV1T, CDUM, 1, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LRWORK-IBBCSD+1, CHILDINFO );

         // Permute rows and columns to place zero submatrices in
         // preferred positions

         if ( Q > 0 && WANTU2 ) {
            for (I = 1; I <= Q; I++) {
               IWORK[I] = M - P - Q + I;
            }
            for (I = Q + 1; I <= M - P; I++) {
               IWORK[I] = I - Q;
            }
            zlapmt( false , M-P, M-P, U2, LDU2, IWORK );
         }
      } else if ( R == P ) {

         // Case 2: R = P

         // Simultaneously bidiagonalize X11 and X21

         zunbdb2(M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO );

         // Accumulate Householder reflectors

         if ( WANTU1 && P > 0 ) {
            U1[1][1] = ONE;
            for (J = 2; J <= P; J++) {
               U1[1][J] = ZERO;
               U1[J][1] = ZERO;
            }
            zlacpy('L', P-1, P-1, X11(2,1), LDX11, U1(2,2), LDU1 );
            zungqr(P-1, P-1, P-1, U1(2,2), LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTU2 && M-P > 0 ) {
            zlacpy('L', M-P, Q, X21, LDX21, U2, LDU2 );
            zungqr(M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTV1T && Q > 0 ) {
            zlacpy('U', P, Q, X11, LDX11, V1T, LDV1T );
            zunglq(Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO );
         }

         // Simultaneously diagonalize X11 and X21.

         zbbcsd(JOBV1T, 'N', JOBU1, JOBU2, 'T', M, Q, P, THETA, RWORK(IPHI), V1T, LDV1T, CDUM, 1, U1, LDU1, U2, LDU2, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LBBCSD, CHILDINFO );

         // Permute rows and columns to place identity submatrices in
         // preferred positions

         if ( Q > 0 && WANTU2 ) {
            for (I = 1; I <= Q; I++) {
               IWORK[I] = M - P - Q + I;
            }
            for (I = Q + 1; I <= M - P; I++) {
               IWORK[I] = I - Q;
            }
            zlapmt( false , M-P, M-P, U2, LDU2, IWORK );
         }
      } else if ( R == M-P ) {

         // Case 3: R = M-P

         // Simultaneously bidiagonalize X11 and X21

         zunbdb3(M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO );

         // Accumulate Householder reflectors

         if ( WANTU1 && P > 0 ) {
            zlacpy('L', P, Q, X11, LDX11, U1, LDU1 );
            zungqr(P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTU2 && M-P > 0 ) {
            U2[1][1] = ONE;
            for (J = 2; J <= M-P; J++) {
               U2[1][J] = ZERO;
               U2[J][1] = ZERO;
            }
            zlacpy('L', M-P-1, M-P-1, X21(2,1), LDX21, U2(2,2), LDU2 );
            zungqr(M-P-1, M-P-1, M-P-1, U2(2,2), LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTV1T && Q > 0 ) {
            zlacpy('U', M-P, Q, X21, LDX21, V1T, LDV1T );
            zunglq(Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO );
         }

         // Simultaneously diagonalize X11 and X21.

         zbbcsd('N', JOBV1T, JOBU2, JOBU1, 'T', M, M-Q, M-P, THETA, RWORK(IPHI), CDUM, 1, V1T, LDV1T, U2, LDU2, U1, LDU1, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LBBCSD, CHILDINFO );

         // Permute rows and columns to place identity submatrices in
         // preferred positions

         if ( Q > R ) {
            for (I = 1; I <= R; I++) {
               IWORK[I] = Q - R + I;
            }
            for (I = R + 1; I <= Q; I++) {
               IWORK[I] = I - R;
            }
            if ( WANTU1 ) {
               zlapmt( false , P, Q, U1, LDU1, IWORK );
            }
            if ( WANTV1T ) {
               zlapmr( false , Q, Q, V1T, LDV1T, IWORK );
            }
         }
      } else {

         // Case 4: R = M-Q

         // Simultaneously bidiagonalize X11 and X21

         zunbdb4(M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), WORK(IORBDB+M), LORBDB-M, CHILDINFO );

         // Accumulate Householder reflectors

         if ( WANTU2 && M-P > 0 ) {
            zcopy(M-P, WORK(IORBDB+P), 1, U2, 1 );
         }
         if ( WANTU1 && P > 0 ) {
            zcopy(P, WORK(IORBDB), 1, U1, 1 );
            for (J = 2; J <= P; J++) {
               U1[1][J] = ZERO;
            }
            zlacpy('L', P-1, M-Q-1, X11(2,1), LDX11, U1(2,2), LDU1 );
            zungqr(P, P, M-Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTU2 && M-P > 0 ) {
            for (J = 2; J <= M-P; J++) {
               U2[1][J] = ZERO;
            }
            zlacpy('L', M-P-1, M-Q-1, X21(2,1), LDX21, U2(2,2), LDU2 );
            zungqr(M-P, M-P, M-Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO );
         }
         if ( WANTV1T && Q > 0 ) {
            zlacpy('U', M-Q, Q, X21, LDX21, V1T, LDV1T );
            zlacpy('U', P-(M-Q), Q-(M-Q), X11(M-Q+1,M-Q+1), LDX11, V1T(M-Q+1,M-Q+1), LDV1T );
            zlacpy('U', -P+Q, Q-P, X21(M-Q+1,P+1), LDX21, V1T(P+1,P+1), LDV1T );
            zunglq(Q, Q, Q, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO );
         }

         // Simultaneously diagonalize X11 and X21.

         zbbcsd(JOBU2, JOBU1, 'N', JOBV1T, 'N', M, M-P, M-Q, THETA, RWORK(IPHI), U2, LDU2, U1, LDU1, CDUM, 1, V1T, LDV1T, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LBBCSD, CHILDINFO );

         // Permute rows and columns to place identity submatrices in
         // preferred positions

         if ( P > R ) {
            for (I = 1; I <= R; I++) {
               IWORK[I] = P - R + I;
            }
            for (I = R + 1; I <= P; I++) {
               IWORK[I] = I - R;
            }
            if ( WANTU1 ) {
               zlapmt( false , P, P, U1, LDU1, IWORK );
            }
            if ( WANTV1T ) {
               zlapmr( false , P, Q, V1T, LDV1T, IWORK );
            }
         }
      }

      }
