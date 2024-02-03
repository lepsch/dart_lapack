      SUBROUTINE CCKGQR( NM, MVAL, NP, PVAL, NN, NVAL, NMATS, ISEED, THRESH, NMAX, A, AF, AQ, AR, TAUA, B, BF, BZ, BT, BWK, TAUB, WORK, RWORK, NIN, NOUT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, NIN, NM, NMATS, NMAX, NN, NOUT, NP
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * )
      REAL               RWORK( * )
      COMPLEX            A( * ), AF( * ), AQ( * ), AR( * ), B( * ), BF( * ), BT( * ), BWK( * ), BZ( * ), TAUA( * ), TAUB( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NTESTS
      PARAMETER          ( NTESTS = 7 )
      int                NTYPES
      PARAMETER          ( NTYPES = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRSTT
      CHARACTER          DISTA, DISTB, TYPE
      CHARACTER*3        PATH
      int                I, IINFO, IM, IMAT, IN, IP, KLA, KLB, KUA, KUB, LDA, LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, NT, P
      REAL               ANORM, BNORM, CNDNMA, CNDNMB
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( NTYPES )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALAREQ, ALASUM, CGQRTS, CGRQTS, CLATMS, SLATB9
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Initialize constants.
*
      PATH( 1: 3 ) = 'GQR'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
      LDA = NMAX
      LDB = NMAX
      LWORK = NMAX*NMAX
*
*     Do for each value of M in MVAL.
*
      DO 60 IM = 1, NM
         M = MVAL( IM )
*
*        Do for each value of P in PVAL.
*
         DO 50 IP = 1, NP
            P = PVAL( IP )
*
*           Do for each value of N in NVAL.
*
            DO 40 IN = 1, NN
               N = NVAL( IN )
*
               DO 30 IMAT = 1, NTYPES
*
*                 Do the tests only if DOTYPE( IMAT ) is true.
*
                  IF( .NOT.DOTYPE( IMAT ) ) GO TO 30
*
*                 Test CGGRQF
*
*                 Set up parameters with SLATB9 and generate test
*                 matrices A and B with CLATMS.
*
                  CALL SLATB9( 'GRQ', IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB )
*
                  CALL CLATMS( M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9999 )IINFO
                     INFO = ABS( IINFO )
                     GO TO 30
                  END IF
*
                  CALL CLATMS( P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9999 )IINFO
                     INFO = ABS( IINFO )
                     GO TO 30
                  END IF
*
                  NT = 4
*
                  CALL CGRQTS( M, P, N, A, AF, AQ, AR, LDA, TAUA, B, BF, BZ, BT, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
*
*                 Print information about the tests that did not
*                 pass the threshold.
*
                  DO 10 I = 1, NT
                     IF( RESULT( I ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                           FIRSTT = .FALSE.
                           CALL ALAHDG( NOUT, 'GRQ' )
                        END IF
                        WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I, RESULT( I )
                        NFAIL = NFAIL + 1
                     END IF
   10             CONTINUE
                  NRUN = NRUN + NT
*
*                 Test CGGQRF
*
*                 Set up parameters with SLATB9 and generate test
*                 matrices A and B with CLATMS.
*
                  CALL SLATB9( 'GQR', IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB )
*
                  CALL CLATMS( N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9999 )IINFO
                     INFO = ABS( IINFO )
                     GO TO 30
                  END IF
*
                  CALL CLATMS( N, P, DISTB, ISEED, TYPE, RWORK, MODEA, CNDNMA, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9999 )IINFO
                     INFO = ABS( IINFO )
                     GO TO 30
                  END IF
*
                  NT = 4
*
                  CALL CGQRTS( N, M, P, A, AF, AQ, AR, LDA, TAUA, B, BF, BZ, BT, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
*
*                 Print information about the tests that did not
*                 pass the threshold.
*
                  DO 20 I = 1, NT
                     IF( RESULT( I ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                           FIRSTT = .FALSE.
                           CALL ALAHDG( NOUT, PATH )
                        END IF
                        WRITE( NOUT, FMT = 9997 )N, M, P, IMAT, I, RESULT( I )
                        NFAIL = NFAIL + 1
                     END IF
   20             CONTINUE
                  NRUN = NRUN + NT
*
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
*
 9999 FORMAT( ' CLATMS in CCKGQR:    INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' N=', I4, ' M=', I4, ', P=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
      RETURN
*
*     End of CCKGQR
*
      END
