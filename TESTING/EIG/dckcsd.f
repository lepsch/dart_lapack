      SUBROUTINE DCKCSD( NM, MVAL, PVAL, QVAL, NMATS, ISEED, THRESH, MMAX, X, XF, U1, U2, V1T, V2T, THETA, IWORK, WORK, RWORK, NIN, NOUT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, NIN, NM, NMATS, MMAX, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 ), IWORK( * ), MVAL( * ), PVAL( * ), QVAL( * )
      DOUBLE PRECISION   RWORK( * ), THETA( * )
      DOUBLE PRECISION   U1( * ), U2( * ), V1T( * ), V2T( * ), WORK( * ), X( * ), XF( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NTESTS
      PARAMETER          ( NTESTS = 15 )
      int                NTYPES
      PARAMETER          ( NTYPES = 4 )
      DOUBLE PRECISION   GAPDIGIT, ONE, ORTH, TEN, ZERO
      PARAMETER          ( GAPDIGIT = 18.0D0, ONE = 1.0D0, ORTH = 1.0D-12, TEN = 10.0D0, ZERO = 0.0D0 )
      DOUBLE PRECISION   PIOVER2
      PARAMETER ( PIOVER2 = 1.57079632679489661923132169163975144210D0 )
*     ..
*     .. Local Scalars ..
      bool               FIRSTT;
      String             PATH;
      int                I, IINFO, IM, IMAT, J, LDU1, LDU2, LDV1T, LDV2T, LDX, LWORK, M, NFAIL, NRUN, NT, P, Q, R
*     ..
*     .. Local Arrays ..
      bool               DOTYPE( NTYPES );
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALAREQ, ALASUM, DCSDTS, DLACSG, DLAROR, DLASET, DROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLARAN, DLARND
      EXTERNAL           DLARAN, DLARND
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 3 ) = 'CSD'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
      LDX = MMAX
      LDU1 = MMAX
      LDU2 = MMAX
      LDV1T = MMAX
      LDV2T = MMAX
      LWORK = MMAX*MMAX
*
*     Do for each value of M in MVAL.
*
      DO 30 IM = 1, NM
         M = MVAL( IM )
         P = PVAL( IM )
         Q = QVAL( IM )
*
         DO 20 IMAT = 1, NTYPES
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) ) GO TO 20
*
*           Generate X
*
            IF( IMAT.EQ.1 ) THEN
               CALL DLAROR( 'L', 'I', M, M, X, LDX, ISEED, WORK, IINFO )
               IF( M .NE. 0 .AND. IINFO .NE. 0 ) THEN
                  WRITE( NOUT, FMT = 9999 ) M, IINFO
                  INFO = ABS( IINFO )
                  GO TO 20
               END IF
            ELSE IF( IMAT.EQ.2 ) THEN
               R = MIN( P, M-P, Q, M-Q )
               DO I = 1, R
                  THETA(I) = PIOVER2 * DLARND( 1, ISEED )
               END DO
               CALL DLACSG( M, P, Q, THETA, ISEED, X, LDX, WORK )
               DO I = 1, M
                  DO J = 1, M
                     X(I+(J-1)*LDX) = X(I+(J-1)*LDX) + ORTH*DLARND(2,ISEED)
                  END DO
               END DO
            ELSE IF( IMAT.EQ.3 ) THEN
               R = MIN( P, M-P, Q, M-Q )
               DO I = 1, R+1
                  THETA(I) = TEN**(-DLARND(1,ISEED)*GAPDIGIT)
               END DO
               DO I = 2, R+1
                  THETA(I) = THETA(I-1) + THETA(I)
               END DO
               DO I = 1, R
                  THETA(I) = PIOVER2 * THETA(I) / THETA(R+1)
               END DO
               CALL DLACSG( M, P, Q, THETA, ISEED, X, LDX, WORK )
            ELSE
               CALL DLASET( 'F', M, M, ZERO, ONE, X, LDX )
               DO I = 1, M
                  J = INT( DLARAN( ISEED ) * M ) + 1
                  IF( J .NE. I ) THEN
                     CALL DROT( M, X(1+(I-1)*LDX), 1, X(1+(J-1)*LDX), 1, ZERO, ONE )
                  END IF
               END DO
            END IF
*
            NT = 15
*
            CALL DCSDTS( M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK, RWORK, RESULT )
*
*           Print information about the tests that did not
*           pass the threshold.
*
            DO 10 I = 1, NT
               IF( RESULT( I ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                     FIRSTT = .FALSE.
                     CALL ALAHDG( NOUT, PATH )
                  END IF
                  WRITE( NOUT, FMT = 9998 )M, P, Q, IMAT, I, RESULT( I )
                  NFAIL = NFAIL + 1
               END IF
   10       CONTINUE
            NRUN = NRUN + NT
   20    CONTINUE
   30 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
*
 9999 FORMAT( ' DLAROR in DCKCSD: M = ', I5, ', INFO = ', I15 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', Q=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
      RETURN
*
*     End of DCKCSD
*
      END
*
*
*
      SUBROUTINE DLACSG( M, P, Q, THETA, ISEED, X, LDX, WORK )
      IMPLICIT NONE
*
      int                LDX, M, P, Q
      int                ISEED( 4 )
      DOUBLE PRECISION   THETA( * )
      DOUBLE PRECISION   WORK( * ), X( LDX, * )
*
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
*
      int                I, INFO, R
*
      R = MIN( P, M-P, Q, M-Q )
*
      CALL DLASET( 'Full', M, M, ZERO, ZERO, X, LDX )
*
      DO I = 1, MIN(P,Q)-R
         X(I,I) = ONE
      END DO
      DO I = 1, R
         X(MIN(P,Q)-R+I,MIN(P,Q)-R+I) = COS(THETA(I))
      END DO
      DO I = 1, MIN(P,M-Q)-R
         X(P-I+1,M-I+1) = -ONE
      END DO
      DO I = 1, R
         X(P-(MIN(P,M-Q)-R)+1-I,M-(MIN(P,M-Q)-R)+1-I) = -SIN(THETA(R-I+1))
      END DO
      DO I = 1, MIN(M-P,Q)-R
         X(M-I+1,Q-I+1) = ONE
      END DO
      DO I = 1, R
         X(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) = SIN(THETA(R-I+1))
      END DO
      DO I = 1, MIN(M-P,M-Q)-R
         X(P+I,Q+I) = ONE
      END DO
      DO I = 1, R
         X(P+(MIN(M-P,M-Q)-R)+I,Q+(MIN(M-P,M-Q)-R)+I) = COS(THETA(I))
      END DO
      CALL DLAROR( 'Left', 'No init', P, M, X, LDX, ISEED, WORK, INFO )
      CALL DLAROR( 'Left', 'No init', M-P, M, X(P+1,1), LDX, ISEED, WORK, INFO )       CALL DLAROR( 'Right', 'No init', M, Q, X, LDX, ISEED, WORK, INFO )       CALL DLAROR( 'Right', 'No init', M, M-Q, X(1,Q+1), LDX, ISEED, WORK, INFO )
*
      END
