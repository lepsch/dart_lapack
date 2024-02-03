      SUBROUTINE DCKGSV( NM, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH, NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R, IWORK, WORK, RWORK, NIN, NOUT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, NIN, NM, NMATS, NMAX, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ), PVAL( * )       DOUBLE PRECISION   A( * ), AF( * ), ALPHA( * ), B( * ), BETA( * ), BF( * ), Q( * ), R( * ), RWORK( * ), U( * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NTESTS
      PARAMETER          ( NTESTS = 12 )
      int                NTYPES
      PARAMETER          ( NTYPES = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRSTT
      CHARACTER          DISTA, DISTB, TYPE
      String             PATH;
      int                I, IINFO, IM, IMAT, KLA, KLB, KUA, KUB, LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, NT, P
      DOUBLE PRECISION   ANORM, BNORM, CNDNMA, CNDNMB
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( NTYPES )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALAREQ, ALASUM, DGSVTS3, DLATB9, DLATMS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 3 ) = 'GSV'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
      LDA = NMAX
      LDB = NMAX
      LDU = NMAX
      LDV = NMAX
      LDQ = NMAX
      LDR = NMAX
      LWORK = NMAX*NMAX
*
*     Do for each value of M in MVAL.
*
      DO 30 IM = 1, NM
         M = MVAL( IM )
         P = PVAL( IM )
         N = NVAL( IM )
*
         DO 20 IMAT = 1, NTYPES
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) ) GO TO 20
*
*           Set up parameters with DLATB9 and generate test
*           matrices A and B with DLATMS.
*
            CALL DLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB )
*
*           Generate M by N matrix A
*
            CALL DLATMS( M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 20
            END IF
*
            CALL DLATMS( P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 20
            END IF
*
            NT = 6
*
            CALL DGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, LWORK, RWORK, RESULT )
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
                  WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I, RESULT( I )
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
 9999 FORMAT( ' DLATMS in DCKGSV   INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
      RETURN
*
*     End of DCKGSV
*
      END
