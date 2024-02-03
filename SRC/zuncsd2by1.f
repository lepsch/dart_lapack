      SUBROUTINE ZUNCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11, X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, WORK, LWORK, RWORK, LRWORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBU1, JOBU2, JOBV1T;
      int                INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21, M, P, Q;
      int                LRWORK, LRWORKMIN, LRWORKOPT;
*     ..
*     .. Array Arguments ..
      double             RWORK(*);
      double             THETA(*);
      COMPLEX*16         U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*), X11(LDX11,*), X21(LDX21,*)
      int                IWORK(*);
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = (1.0D0,0.0D0), ZERO = (0.0D0,0.0D0) )
*     ..
*     .. Local Scalars ..
      int                CHILDINFO, I, IB11D, IB11E, IB12D, IB12E, IB21D, IB21E, IB22D, IB22E, IBBCSD, IORBDB, IORGLQ, IORGQR, IPHI, ITAUP1, ITAUP2, ITAUQ1, J, LBBCSD, LORBDB, LORGLQ, LORGLQMIN, LORGLQOPT, LORGQR, LORGQRMIN, LORGQROPT, LWORKMIN, LWORKOPT, R;
      bool               LQUERY, WANTU1, WANTU2, WANTV1T;
*     ..
*     .. Local Arrays ..
      double             DUM( 1 );
      COMPLEX*16         CDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZBBCSD, ZCOPY, ZLACPY, ZLAPMR, ZLAPMT, ZUNBDB1, ZUNBDB2, ZUNBDB3, ZUNBDB4, ZUNGLQ, ZUNGQR, XERBLA
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Function ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test input arguments
*
      INFO = 0
      WANTU1 = LSAME( JOBU1, 'Y' )
      WANTU2 = LSAME( JOBU2, 'Y' )
      WANTV1T = LSAME( JOBV1T, 'Y' )
      LQUERY = ( LWORK.EQ.-1 ) .OR. ( LRWORK.EQ.-1 )
*
      IF( M .LT. 0 ) THEN
         INFO = -4
      ELSE IF( P .LT. 0 .OR. P .GT. M ) THEN
         INFO = -5
      ELSE IF( Q .LT. 0 .OR. Q .GT. M ) THEN
         INFO = -6
      ELSE IF( LDX11 .LT. MAX( 1, P ) ) THEN
         INFO = -8
      ELSE IF( LDX21 .LT. MAX( 1, M-P ) ) THEN
         INFO = -10
      ELSE IF( WANTU1 .AND. LDU1 .LT. MAX( 1, P ) ) THEN
         INFO = -13
      ELSE IF( WANTU2 .AND. LDU2 .LT. MAX( 1, M - P ) ) THEN
         INFO = -15
      ELSE IF( WANTV1T .AND. LDV1T .LT. MAX( 1, Q ) ) THEN
         INFO = -17
      END IF
*
      R = MIN( P, M-P, Q, M-Q )
*
*     Compute workspace
*
*       WORK layout:
*     |-----------------------------------------|
*     | LWORKOPT (1)                            |
*     |-----------------------------------------|
*     | TAUP1 (MAX(1,P))                        |
*     | TAUP2 (MAX(1,M-P))                      |
*     | TAUQ1 (MAX(1,Q))                        |
*     |-----------------------------------------|
*     | ZUNBDB WORK | ZUNGQR WORK | ZUNGLQ WORK |
*     |             |             |             |
*     |             |             |             |
*     |             |             |             |
*     |             |             |             |
*     |-----------------------------------------|
*       RWORK layout:
*     |------------------|
*     | LRWORKOPT (1)    |
*     |------------------|
*     | PHI (MAX(1,R-1)) |
*     |------------------|
*     | B11D (R)         |
*     | B11E (R-1)       |
*     | B12D (R)         |
*     | B12E (R-1)       |
*     | B21D (R)         |
*     | B21E (R-1)       |
*     | B22D (R)         |
*     | B22E (R-1)       |
*     | ZBBCSD RWORK     |
*     |------------------|
*
      IF( INFO .EQ. 0 ) THEN
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
         ITAUP1 = 2
         ITAUP2 = ITAUP1 + MAX( 1, P )
         ITAUQ1 = ITAUP2 + MAX( 1, M-P )
         IORBDB = ITAUQ1 + MAX( 1, Q )
         IORGQR = ITAUQ1 + MAX( 1, Q )
         IORGLQ = ITAUQ1 + MAX( 1, Q )
         LORGQRMIN = 1
         LORGQROPT = 1
         LORGLQMIN = 1
         LORGLQOPT = 1
         IF( R .EQ. Q ) THEN
            CALL ZUNBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM, CDUM, CDUM, WORK, -1, CHILDINFO )
            LORBDB = INT( WORK(1) )
            IF( WANTU1 .AND. P .GT. 0 ) THEN
               CALL ZUNGQR( P, P, Q, U1, LDU1, CDUM, WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            ENDIF
            IF( WANTU2 .AND. M-P .GT. 0 ) THEN
               CALL ZUNGQR( M-P, M-P, Q, U2, LDU2, CDUM, WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTV1T .AND. Q .GT. 0 ) THEN
               CALL ZUNGLQ( Q-1, Q-1, Q-1, V1T, LDV1T, CDUM, WORK(1), -1, CHILDINFO )
               LORGLQMIN = MAX( LORGLQMIN, Q-1 )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            END IF
            CALL ZBBCSD( JOBU1, JOBU2, JOBV1T, 'N', 'N', M, P, Q, THETA, DUM, U1, LDU1, U2, LDU2, V1T, LDV1T, CDUM, 1, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, RWORK(1), -1, CHILDINFO )
            LBBCSD = INT( RWORK(1) )
         ELSE IF( R .EQ. P ) THEN
            CALL ZUNBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM, CDUM, CDUM, WORK(1), -1, CHILDINFO )
            LORBDB = INT( WORK(1) )
            IF( WANTU1 .AND. P .GT. 0 ) THEN
               CALL ZUNGQR( P-1, P-1, P-1, U1(2,2), LDU1, CDUM, WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, P-1 )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTU2 .AND. M-P .GT. 0 ) THEN
               CALL ZUNGQR( M-P, M-P, Q, U2, LDU2, CDUM, WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTV1T .AND. Q .GT. 0 ) THEN
               CALL ZUNGLQ( Q, Q, R, V1T, LDV1T, CDUM, WORK(1), -1, CHILDINFO )
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            END IF
            CALL ZBBCSD( JOBV1T, 'N', JOBU1, JOBU2, 'T', M, Q, P, THETA, DUM, V1T, LDV1T, CDUM, 1, U1, LDU1, U2, LDU2, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, RWORK(1), -1, CHILDINFO )
            LBBCSD = INT( RWORK(1) )
         ELSE IF( R .EQ. M-P ) THEN
            CALL ZUNBDB3( M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM, CDUM, CDUM, WORK(1), -1, CHILDINFO )
            LORBDB = INT( WORK(1) )
            IF( WANTU1 .AND. P .GT. 0 ) THEN
               CALL ZUNGQR( P, P, Q, U1, LDU1, CDUM, WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTU2 .AND. M-P .GT. 0 ) THEN
               CALL ZUNGQR( M-P-1, M-P-1, M-P-1, U2(2,2), LDU2, CDUM, WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, M-P-1 )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTV1T .AND. Q .GT. 0 ) THEN
               CALL ZUNGLQ( Q, Q, R, V1T, LDV1T, CDUM, WORK(1), -1, CHILDINFO )
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            END IF
            CALL ZBBCSD( 'N', JOBV1T, JOBU2, JOBU1, 'T', M, M-Q, M-P, THETA, DUM, CDUM, 1, V1T, LDV1T, U2, LDU2, U1, LDU1, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, RWORK(1), -1, CHILDINFO )
            LBBCSD = INT( RWORK(1) )
         ELSE
            CALL ZUNBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, DUM, CDUM, CDUM, CDUM, CDUM, WORK(1), -1, CHILDINFO )
            LORBDB = M + INT( WORK(1) )
            IF( WANTU1 .AND. P .GT. 0 ) THEN
               CALL ZUNGQR( P, P, M-Q, U1, LDU1, CDUM, WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTU2 .AND. M-P .GT. 0 ) THEN
               CALL ZUNGQR( M-P, M-P, M-Q, U2, LDU2, CDUM, WORK(1), -1, CHILDINFO )
               LORGQRMIN = MAX( LORGQRMIN, M-P )
               LORGQROPT = MAX( LORGQROPT, INT( WORK(1) ) )
            END IF
            IF( WANTV1T .AND. Q .GT. 0 ) THEN
               CALL ZUNGLQ( Q, Q, Q, V1T, LDV1T, CDUM, WORK(1), -1, CHILDINFO )
               LORGLQMIN = MAX( LORGLQMIN, Q )
               LORGLQOPT = MAX( LORGLQOPT, INT( WORK(1) ) )
            END IF
            CALL ZBBCSD( JOBU2, JOBU1, 'N', JOBV1T, 'N', M, M-P, M-Q, THETA, DUM, U2, LDU2, U1, LDU1, CDUM, 1, V1T, LDV1T, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, RWORK(1), -1, CHILDINFO )
            LBBCSD = INT( RWORK(1) )
         END IF
         LRWORKMIN = IBBCSD+LBBCSD-1
         LRWORKOPT = LRWORKMIN
         RWORK(1) = LRWORKOPT
         LWORKMIN = MAX( IORBDB+LORBDB-1, IORGQR+LORGQRMIN-1, IORGLQ+LORGLQMIN-1 )          LWORKOPT = MAX( IORBDB+LORBDB-1, IORGQR+LORGQROPT-1, IORGLQ+LORGLQOPT-1 )
         WORK(1) = LWORKOPT
         IF( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -19
         END IF
         IF( LRWORK .LT. LRWORKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -21
         END IF
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'ZUNCSD2BY1', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
      LORGQR = LWORK-IORGQR+1
      LORGLQ = LWORK-IORGLQ+1
*
*     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q,
*     in which R = MIN(P,M-P,Q,M-Q)
*
      IF( R .EQ. Q ) THEN
*
*        Case 1: R = Q
*
*        Simultaneously bidiagonalize X11 and X21
*
         CALL ZUNBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO )
*
*        Accumulate Householder reflectors
*
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            CALL ZLACPY( 'L', P, Q, X11, LDX11, U1, LDU1 )
            CALL ZUNGQR( P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            CALL ZLACPY( 'L', M-P, Q, X21, LDX21, U2, LDU2 )
            CALL ZUNGQR( M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            V1T(1,1) = ONE
            DO J = 2, Q
               V1T(1,J) = ZERO
               V1T(J,1) = ZERO
            END DO
            CALL ZLACPY( 'U', Q-1, Q-1, X21(1,2), LDX21, V1T(2,2), LDV1T )             CALL ZUNGLQ( Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO )
         END IF
*
*        Simultaneously diagonalize X11 and X21.
*
         CALL ZBBCSD( JOBU1, JOBU2, JOBV1T, 'N', 'N', M, P, Q, THETA, RWORK(IPHI), U1, LDU1, U2, LDU2, V1T, LDV1T, CDUM, 1, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LRWORK-IBBCSD+1, CHILDINFO )
*
*        Permute rows and columns to place zero submatrices in
*        preferred positions
*
         IF( Q .GT. 0 .AND. WANTU2 ) THEN
            DO I = 1, Q
               IWORK(I) = M - P - Q + I
            END DO
            DO I = Q + 1, M - P
               IWORK(I) = I - Q
            END DO
            CALL ZLAPMT( .FALSE., M-P, M-P, U2, LDU2, IWORK )
         END IF
      ELSE IF( R .EQ. P ) THEN
*
*        Case 2: R = P
*
*        Simultaneously bidiagonalize X11 and X21
*
         CALL ZUNBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO )
*
*        Accumulate Householder reflectors
*
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            U1(1,1) = ONE
            DO J = 2, P
               U1(1,J) = ZERO
               U1(J,1) = ZERO
            END DO
            CALL ZLACPY( 'L', P-1, P-1, X11(2,1), LDX11, U1(2,2), LDU1 )
            CALL ZUNGQR( P-1, P-1, P-1, U1(2,2), LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            CALL ZLACPY( 'L', M-P, Q, X21, LDX21, U2, LDU2 )
            CALL ZUNGQR( M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            CALL ZLACPY( 'U', P, Q, X11, LDX11, V1T, LDV1T )
            CALL ZUNGLQ( Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO )
         END IF
*
*        Simultaneously diagonalize X11 and X21.
*
         CALL ZBBCSD( JOBV1T, 'N', JOBU1, JOBU2, 'T', M, Q, P, THETA, RWORK(IPHI), V1T, LDV1T, CDUM, 1, U1, LDU1, U2, LDU2, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LBBCSD, CHILDINFO )
*
*        Permute rows and columns to place identity submatrices in
*        preferred positions
*
         IF( Q .GT. 0 .AND. WANTU2 ) THEN
            DO I = 1, Q
               IWORK(I) = M - P - Q + I
            END DO
            DO I = Q + 1, M - P
               IWORK(I) = I - Q
            END DO
            CALL ZLAPMT( .FALSE., M-P, M-P, U2, LDU2, IWORK )
         END IF
      ELSE IF( R .EQ. M-P ) THEN
*
*        Case 3: R = M-P
*
*        Simultaneously bidiagonalize X11 and X21
*
         CALL ZUNBDB3( M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), LORBDB, CHILDINFO )
*
*        Accumulate Householder reflectors
*
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            CALL ZLACPY( 'L', P, Q, X11, LDX11, U1, LDU1 )
            CALL ZUNGQR( P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            U2(1,1) = ONE
            DO J = 2, M-P
               U2(1,J) = ZERO
               U2(J,1) = ZERO
            END DO
            CALL ZLACPY( 'L', M-P-1, M-P-1, X21(2,1), LDX21, U2(2,2), LDU2 )             CALL ZUNGQR( M-P-1, M-P-1, M-P-1, U2(2,2), LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            CALL ZLACPY( 'U', M-P, Q, X21, LDX21, V1T, LDV1T )
            CALL ZUNGLQ( Q, Q, R, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO )
         END IF
*
*        Simultaneously diagonalize X11 and X21.
*
         CALL ZBBCSD( 'N', JOBV1T, JOBU2, JOBU1, 'T', M, M-Q, M-P, THETA, RWORK(IPHI), CDUM, 1, V1T, LDV1T, U2, LDU2, U1, LDU1, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LBBCSD, CHILDINFO )
*
*        Permute rows and columns to place identity submatrices in
*        preferred positions
*
         IF( Q .GT. R ) THEN
            DO I = 1, R
               IWORK(I) = Q - R + I
            END DO
            DO I = R + 1, Q
               IWORK(I) = I - R
            END DO
            IF( WANTU1 ) THEN
               CALL ZLAPMT( .FALSE., P, Q, U1, LDU1, IWORK )
            END IF
            IF( WANTV1T ) THEN
               CALL ZLAPMR( .FALSE., Q, Q, V1T, LDV1T, IWORK )
            END IF
         END IF
      ELSE
*
*        Case 4: R = M-Q
*
*        Simultaneously bidiagonalize X11 and X21
*
         CALL ZUNBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, RWORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(IORBDB), WORK(IORBDB+M), LORBDB-M, CHILDINFO )
*
*        Accumulate Householder reflectors
*
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            CALL ZCOPY( M-P, WORK(IORBDB+P), 1, U2, 1 )
         END IF
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            CALL ZCOPY( P, WORK(IORBDB), 1, U1, 1 )
            DO J = 2, P
               U1(1,J) = ZERO
            END DO
            CALL ZLACPY( 'L', P-1, M-Q-1, X11(2,1), LDX11, U1(2,2), LDU1 )             CALL ZUNGQR( P, P, M-Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            DO J = 2, M-P
               U2(1,J) = ZERO
            END DO
            CALL ZLACPY( 'L', M-P-1, M-Q-1, X21(2,1), LDX21, U2(2,2), LDU2 )             CALL ZUNGQR( M-P, M-P, M-Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQR, CHILDINFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            CALL ZLACPY( 'U', M-Q, Q, X21, LDX21, V1T, LDV1T )
            CALL ZLACPY( 'U', P-(M-Q), Q-(M-Q), X11(M-Q+1,M-Q+1), LDX11, V1T(M-Q+1,M-Q+1), LDV1T )             CALL ZLACPY( 'U', -P+Q, Q-P, X21(M-Q+1,P+1), LDX21, V1T(P+1,P+1), LDV1T )             CALL ZUNGLQ( Q, Q, Q, V1T, LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQ, CHILDINFO )
         END IF
*
*        Simultaneously diagonalize X11 and X21.
*
         CALL ZBBCSD( JOBU2, JOBU1, 'N', JOBV1T, 'N', M, M-P, M-Q, THETA, RWORK(IPHI), U2, LDU2, U1, LDU1, CDUM, 1, V1T, LDV1T, RWORK(IB11D), RWORK(IB11E), RWORK(IB12D), RWORK(IB12E), RWORK(IB21D), RWORK(IB21E), RWORK(IB22D), RWORK(IB22E), RWORK(IBBCSD), LBBCSD, CHILDINFO )
*
*        Permute rows and columns to place identity submatrices in
*        preferred positions
*
         IF( P .GT. R ) THEN
            DO I = 1, R
               IWORK(I) = P - R + I
            END DO
            DO I = R + 1, P
               IWORK(I) = I - R
            END DO
            IF( WANTU1 ) THEN
               CALL ZLAPMT( .FALSE., P, P, U1, LDU1, IWORK )
            END IF
            IF( WANTV1T ) THEN
               CALL ZLAPMR( .FALSE., P, Q, V1T, LDV1T, IWORK )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ZUNCSD2BY1
*
      END
