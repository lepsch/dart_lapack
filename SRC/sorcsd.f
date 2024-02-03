      RECURSIVE SUBROUTINE SORCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK, LWORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS;
      int                INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, Q
*     ..
*     .. Array Arguments ..
      int                IWORK( * )
      REAL               THETA( * )
      REAL               U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ), X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, * )
*     ..
*
*  ===================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Arrays ..
      REAL               DUMMY(1)
*     ..
*     .. Local Scalars ..
      String             TRANST, SIGNST;
      int                CHILDINFO, I, IB11D, IB11E, IB12D, IB12E, IB21D, IB21E, IB22D, IB22E, IBBCSD, IORBDB, IORGLQ, IORGQR, IPHI, ITAUP1, ITAUP2, ITAUQ1, ITAUQ2, J, LBBCSDWORK, LBBCSDWORKMIN, LBBCSDWORKOPT, LORBDBWORK, LORBDBWORKMIN, LORBDBWORKOPT, LORGLQWORK, LORGLQWORKMIN, LORGLQWORKOPT, LORGQRWORK, LORGQRWORKMIN, LORGQRWORKOPT, LWORKMIN, LWORKOPT
      LOGICAL            COLMAJOR, DEFAULTSIGNS, LQUERY, WANTU1, WANTU2, WANTV1T, WANTV2T
*     ..
*     .. External Subroutines ..
      EXTERNAL           SBBCSD, SLACPY, SLAPMR, SLAPMT, SORBDB, SORGLQ, SORGQR, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions
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
      WANTV2T = LSAME( JOBV2T, 'Y' )
      COLMAJOR = .NOT. LSAME( TRANS, 'T' )
      DEFAULTSIGNS = .NOT. LSAME( SIGNS, 'O' )
      LQUERY = LWORK .EQ. -1
      IF( M .LT. 0 ) THEN
         INFO = -7
      ELSE IF( P .LT. 0 .OR. P .GT. M ) THEN
         INFO = -8
      ELSE IF( Q .LT. 0 .OR. Q .GT. M ) THEN
         INFO = -9
      ELSE IF ( COLMAJOR .AND.  LDX11 .LT. MAX( 1, P ) ) THEN
        INFO = -11
      ELSE IF (.NOT. COLMAJOR .AND. LDX11 .LT. MAX( 1, Q ) ) THEN
        INFO = -11
      ELSE IF (COLMAJOR .AND. LDX12 .LT. MAX( 1, P ) ) THEN
        INFO = -13
      ELSE IF (.NOT. COLMAJOR .AND. LDX12 .LT. MAX( 1, M-Q ) ) THEN
        INFO = -13
      ELSE IF (COLMAJOR .AND. LDX21 .LT. MAX( 1, M-P ) ) THEN
        INFO = -15
      ELSE IF (.NOT. COLMAJOR .AND. LDX21 .LT. MAX( 1, Q ) ) THEN
        INFO = -15
      ELSE IF (COLMAJOR .AND. LDX22 .LT. MAX( 1, M-P ) ) THEN
        INFO = -17
      ELSE IF (.NOT. COLMAJOR .AND. LDX22 .LT. MAX( 1, M-Q ) ) THEN
        INFO = -17
      ELSE IF( WANTU1 .AND. LDU1 .LT. P ) THEN
         INFO = -20
      ELSE IF( WANTU2 .AND. LDU2 .LT. M-P ) THEN
         INFO = -22
      ELSE IF( WANTV1T .AND. LDV1T .LT. Q ) THEN
         INFO = -24
      ELSE IF( WANTV2T .AND. LDV2T .LT. M-Q ) THEN
         INFO = -26
      END IF
*
*     Work with transpose if convenient
*
      IF( INFO .EQ. 0 .AND. MIN( P, M-P ) .LT. MIN( Q, M-Q ) ) THEN
         IF( COLMAJOR ) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF
         IF( DEFAULTSIGNS ) THEN
            SIGNST = 'O'
         ELSE
            SIGNST = 'D'
         END IF
         CALL SORCSD( JOBV1T, JOBV2T, JOBU1, JOBU2, TRANST, SIGNST, M, Q, P, X11, LDX11, X21, LDX21, X12, LDX12, X22, LDX22, THETA, V1T, LDV1T, V2T, LDV2T, U1, LDU1, U2, LDU2, WORK, LWORK, IWORK, INFO )
         RETURN
      END IF
*
*     Work with permutation [ 0 I; I 0 ] * X * [ 0 I; I 0 ] if
*     convenient
*
      IF( INFO .EQ. 0 .AND. M-Q .LT. Q ) THEN
         IF( DEFAULTSIGNS ) THEN
            SIGNST = 'O'
         ELSE
            SIGNST = 'D'
         END IF
         CALL SORCSD( JOBU2, JOBU1, JOBV2T, JOBV1T, TRANS, SIGNST, M, M-P, M-Q, X22, LDX22, X21, LDX21, X12, LDX12, X11, LDX11, THETA, U2, LDU2, U1, LDU1, V2T, LDV2T, V1T, LDV1T, WORK, LWORK, IWORK, INFO )
         RETURN
      END IF
*
*     Compute workspace
*
      IF( INFO .EQ. 0 ) THEN
*
         IPHI = 2
         ITAUP1 = IPHI + MAX( 1, Q - 1 )
         ITAUP2 = ITAUP1 + MAX( 1, P )
         ITAUQ1 = ITAUP2 + MAX( 1, M - P )
         ITAUQ2 = ITAUQ1 + MAX( 1, Q )
         IORGQR = ITAUQ2 + MAX( 1, M - Q )
         CALL SORGQR( M-Q, M-Q, M-Q, DUMMY, MAX(1,M-Q), DUMMY, WORK, -1, CHILDINFO )
         LORGQRWORKOPT = INT( WORK(1) )
         LORGQRWORKMIN = MAX( 1, M - Q )
         IORGLQ = ITAUQ2 + MAX( 1, M - Q )
         CALL SORGLQ( M-Q, M-Q, M-Q, DUMMY, MAX(1,M-Q), DUMMY, WORK, -1, CHILDINFO )
         LORGLQWORKOPT = INT( WORK(1) )
         LORGLQWORKMIN = MAX( 1, M - Q )
         IORBDB = ITAUQ2 + MAX( 1, M - Q )
         CALL SORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY,WORK,-1,CHILDINFO )
         LORBDBWORKOPT = INT( WORK(1) )
         LORBDBWORKMIN = LORBDBWORKOPT
         IB11D = ITAUQ2 + MAX( 1, M - Q )
         IB11E = IB11D + MAX( 1, Q )
         IB12D = IB11E + MAX( 1, Q - 1 )
         IB12E = IB12D + MAX( 1, Q )
         IB21D = IB12E + MAX( 1, Q - 1 )
         IB21E = IB21D + MAX( 1, Q )
         IB22D = IB21E + MAX( 1, Q - 1 )
         IB22E = IB22D + MAX( 1, Q )
         IBBCSD = IB22E + MAX( 1, Q - 1 )
         CALL SBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, DUMMY, DUMMY, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, WORK, -1, CHILDINFO )
         LBBCSDWORKOPT = INT( WORK(1) )
         LBBCSDWORKMIN = LBBCSDWORKOPT
         LWORKOPT = MAX( IORGQR + LORGQRWORKOPT, IORGLQ + LORGLQWORKOPT, IORBDB + LORBDBWORKOPT, IBBCSD + LBBCSDWORKOPT ) - 1          LWORKMIN = MAX( IORGQR + LORGQRWORKMIN, IORGLQ + LORGLQWORKMIN, IORBDB + LORBDBWORKOPT, IBBCSD + LBBCSDWORKMIN ) - 1
         WORK(1) = MAX(LWORKOPT,LWORKMIN)
*
         IF( LWORK .LT. LWORKMIN .AND. .NOT. LQUERY ) THEN
            INFO = -22
         ELSE
            LORGQRWORK = LWORK - IORGQR + 1
            LORGLQWORK = LWORK - IORGLQ + 1
            LORBDBWORK = LWORK - IORBDB + 1
            LBBCSDWORK = LWORK - IBBCSD + 1
         END IF
      END IF
*
*     Abort if any illegal arguments
*
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'SORCSD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Transform to bidiagonal block form
*
      CALL SORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, WORK(IPHI), WORK(ITAUP1), WORK(ITAUP2), WORK(ITAUQ1), WORK(ITAUQ2), WORK(IORBDB), LORBDBWORK, CHILDINFO )
*
*     Accumulate Householder reflectors
*
      IF( COLMAJOR ) THEN
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            CALL SLACPY( 'L', P, Q, X11, LDX11, U1, LDU1 )
            CALL SORGQR( P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGQR), LORGQRWORK, INFO)
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            CALL SLACPY( 'L', M-P, Q, X21, LDX21, U2, LDU2 )
            CALL SORGQR( M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGQR), LORGQRWORK, INFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            CALL SLACPY( 'U', Q-1, Q-1, X11(1,2), LDX11, V1T(2,2), LDV1T )
            V1T(1, 1) = ONE
            DO J = 2, Q
               V1T(1,J) = ZERO
               V1T(J,1) = ZERO
            END DO
            CALL SORGLQ( Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGLQ), LORGLQWORK, INFO )
         END IF
         IF( WANTV2T .AND. M-Q .GT. 0 ) THEN
            CALL SLACPY( 'U', P, M-Q, X12, LDX12, V2T, LDV2T )
            CALL SLACPY( 'U', M-P-Q, M-P-Q, X22(Q+1,P+1), LDX22, V2T(P+1,P+1), LDV2T )             CALL SORGLQ( M-Q, M-Q, M-Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGLQ), LORGLQWORK, INFO )
         END IF
      ELSE
         IF( WANTU1 .AND. P .GT. 0 ) THEN
            CALL SLACPY( 'U', Q, P, X11, LDX11, U1, LDU1 )
            CALL SORGLQ( P, P, Q, U1, LDU1, WORK(ITAUP1), WORK(IORGLQ), LORGLQWORK, INFO)
         END IF
         IF( WANTU2 .AND. M-P .GT. 0 ) THEN
            CALL SLACPY( 'U', Q, M-P, X21, LDX21, U2, LDU2 )
            CALL SORGLQ( M-P, M-P, Q, U2, LDU2, WORK(ITAUP2), WORK(IORGLQ), LORGLQWORK, INFO )
         END IF
         IF( WANTV1T .AND. Q .GT. 0 ) THEN
            CALL SLACPY( 'L', Q-1, Q-1, X11(2,1), LDX11, V1T(2,2), LDV1T )
            V1T(1, 1) = ONE
            DO J = 2, Q
               V1T(1,J) = ZERO
               V1T(J,1) = ZERO
            END DO
            CALL SORGQR( Q-1, Q-1, Q-1, V1T(2,2), LDV1T, WORK(ITAUQ1), WORK(IORGQR), LORGQRWORK, INFO )
         END IF
         IF( WANTV2T .AND. M-Q .GT. 0 ) THEN
            CALL SLACPY( 'L', M-Q, P, X12, LDX12, V2T, LDV2T )
            CALL SLACPY( 'L', M-P-Q, M-P-Q, X22(P+1,Q+1), LDX22, V2T(P+1,P+1), LDV2T )             CALL SORGQR( M-Q, M-Q, M-Q, V2T, LDV2T, WORK(ITAUQ2), WORK(IORGQR), LORGQRWORK, INFO )
         END IF
      END IF
*
*     Compute the CSD of the matrix in bidiagonal-block form
*
      CALL SBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q, THETA, WORK(IPHI), U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK(IB11D), WORK(IB11E), WORK(IB12D), WORK(IB12E), WORK(IB21D), WORK(IB21E), WORK(IB22D), WORK(IB22E), WORK(IBBCSD), LBBCSDWORK, INFO )
*
*     Permute rows and columns to place identity submatrices in top-
*     left corner of (1,1)-block and/or bottom-right corner of (1,2)-
*     block and/or bottom-right corner of (2,1)-block and/or top-left
*     corner of (2,2)-block
*
      IF( Q .GT. 0 .AND. WANTU2 ) THEN
         DO I = 1, Q
            IWORK(I) = M - P - Q + I
         END DO
         DO I = Q + 1, M - P
            IWORK(I) = I - Q
         END DO
         IF( COLMAJOR ) THEN
            CALL SLAPMT( .FALSE., M-P, M-P, U2, LDU2, IWORK )
         ELSE
            CALL SLAPMR( .FALSE., M-P, M-P, U2, LDU2, IWORK )
         END IF
      END IF
      IF( M .GT. 0 .AND. WANTV2T ) THEN
         DO I = 1, P
            IWORK(I) = M - P - Q + I
         END DO
         DO I = P + 1, M - Q
            IWORK(I) = I - P
         END DO
         IF( .NOT. COLMAJOR ) THEN
            CALL SLAPMT( .FALSE., M-Q, M-Q, V2T, LDV2T, IWORK )
         ELSE
            CALL SLAPMR( .FALSE., M-Q, M-Q, V2T, LDV2T, IWORK )
         END IF
      END IF
*
      RETURN
*
*     End SORCSD
*
      END
