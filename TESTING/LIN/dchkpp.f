      SUBROUTINE DCHKPP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      int                NMAX, NN, NNS, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      int                IWORK( * ), NSVAL( * ), NVAL( * )
      DOUBLE PRECISION   A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      int                NTYPES
      PARAMETER          ( NTYPES = 9 )
      int                NTESTS
      PARAMETER          ( NTESTS = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      CHARACTER          DIST, PACKIT, TYPE, UPLO, XTYPE
      String             PATH;
      int                I, IMAT, IN, INFO, IOFF, IRHS, IUPLO, IZERO, K, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NPP, NRHS, NRUN
      DOUBLE PRECISION   ANORM, CNDNUM, RCOND, RCONDC
*     ..
*     .. Local Arrays ..
      CHARACTER          PACKS( 2 ), UPLOS( 2 )
      int                ISEED( 4 ), ISEEDY( 4 )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DGET06, DLANSP
      EXTERNAL           DGET06, DLANSP
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, DCOPY, DERRPO, DGET04, DLACPY, DLARHS, DLATB4, DLATMS, DPPCON, DPPRFS, DPPT01, DPPT02, DPPT03, DPPT05, DPPTRF, DPPTRI, DPPTRS
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      String             SRNAMT;
      int                INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , PACKS / 'C', 'R' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'PP'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR ) CALL DERRPO( PATH, NOUT )
      INFOT = 0
*
*     Do for each value of N in NVAL
*
      DO 110 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1
*
         DO 100 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) ) GO TO 100
*
*           Skip types 3, 4, or 5 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 100
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 90 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
               PACKIT = PACKS( IUPLO )
*
*              Set up parameters with DLATB4 and generate a test matrix
*              with DLATMS.
*
               CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
               SRNAMT = 'DLATMS'
               CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO )
*
*              Check error code from DLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 90
               END IF
*
*              For types 3-5, zero one row and column of the matrix to
*              test that INFO is returned correctly.
*
               IF( ZEROT ) THEN
                  IF( IMAT.EQ.3 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT.EQ.4 ) THEN
                     IZERO = N
                  ELSE
                     IZERO = N / 2 + 1
                  END IF
*
*                 Set row and column IZERO of A to 0.
*
                  IF( IUPLO.EQ.1 ) THEN
                     IOFF = ( IZERO-1 )*IZERO / 2
                     DO 20 I = 1, IZERO - 1
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                     IOFF = IOFF + IZERO
                     DO 30 I = IZERO, N
                        A( IOFF ) = ZERO
                        IOFF = IOFF + I
   30                CONTINUE
                  ELSE
                     IOFF = IZERO
                     DO 40 I = 1, IZERO - 1
                        A( IOFF ) = ZERO
                        IOFF = IOFF + N - I
   40                CONTINUE
                     IOFF = IOFF - IZERO
                     DO 50 I = IZERO, N
                        A( IOFF+I ) = ZERO
   50                CONTINUE
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
*              Compute the L*L' or U'*U factorization of the matrix.
*
               NPP = N*( N+1 ) / 2
               CALL DCOPY( NPP, A, 1, AFAC, 1 )
               SRNAMT = 'DPPTRF'
               CALL DPPTRF( UPLO, N, AFAC, INFO )
*
*              Check error code from DPPTRF.
*
               IF( INFO.NE.IZERO ) THEN
                  CALL ALAERH( PATH, 'DPPTRF', INFO, IZERO, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 90
               END IF
*
*              Skip the tests if INFO is not 0.
*
               IF( INFO.NE.0 ) GO TO 90
*
*+    TEST 1
*              Reconstruct matrix from factors and compute residual.
*
               CALL DCOPY( NPP, AFAC, 1, AINV, 1 )
               CALL DPPT01( UPLO, N, A, AINV, RWORK, RESULT( 1 ) )
*
*+    TEST 2
*              Form the inverse and compute the residual.
*
               CALL DCOPY( NPP, AFAC, 1, AINV, 1 )
               SRNAMT = 'DPPTRI'
               CALL DPPTRI( UPLO, N, AINV, INFO )
*
*              Check error code from DPPTRI.
*
               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DPPTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
               CALL DPPT03( UPLO, N, A, AINV, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) )
*
*              Print information about the tests that did not pass
*              the threshold.
*
               DO 60 K = 1, 2
                  IF( RESULT( K ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, K, RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
   60          CONTINUE
               NRUN = NRUN + 2
*
               DO 80 IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
*
*+    TEST 3
*              Solve and compute residual for  A * X = B.
*
                  SRNAMT = 'DLARHS'
                  CALL DLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                  CALL DLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                  SRNAMT = 'DPPTRS'
                  CALL DPPTRS( UPLO, N, NRHS, AFAC, X, LDA, INFO )
*
*              Check error code from DPPTRS.
*
                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DPPTRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                  CALL DLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                  CALL DPPT02( UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) )
*
*+    TEST 4
*              Check solution from generated exact solution.
*
                  CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )
*
*+    TESTS 5, 6, and 7
*              Use iterative refinement to improve the solution.
*
                  SRNAMT = 'DPPRFS'
                  CALL DPPRFS( UPLO, N, NRHS, A, AFAC, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO )
*
*              Check error code from DPPRFS.
*
                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DPPRFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                  CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) )                   CALL DPPT05( UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 6 ) )
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 70 K = 3, 7
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
   70             CONTINUE
                  NRUN = NRUN + 5
   80          CONTINUE
*
*+    TEST 8
*              Get an estimate of RCOND = 1/CNDNUM.
*
               ANORM = DLANSP( '1', UPLO, N, A, RWORK )
               SRNAMT = 'DPPCON'
               CALL DPPCON( UPLO, N, AFAC, ANORM, RCOND, WORK, IWORK, INFO )
*
*              Check error code from DPPCON.
*
               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DPPCON', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
               RESULT( 8 ) = DGET06( RCOND, RCONDC )
*
*              Print the test ratio if greater than or equal to THRESH.
*
               IF( RESULT( 8 ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                   WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, 8, RESULT( 8 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', type ', I2, ', test ',
     $      I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
      RETURN
*
*     End of DCHKPP
*
      END
