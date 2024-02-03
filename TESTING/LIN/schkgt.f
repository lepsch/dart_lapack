      SUBROUTINE SCHKGT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A, AF, B, X, XACT, WORK, RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      int                NN, NNS, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      int                IWORK( * ), NSVAL( * ), NVAL( * )
      REAL               A( * ), AF( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      int                NTYPES
      PARAMETER          ( NTYPES = 12 )
      int                NTESTS
      PARAMETER          ( NTESTS = 7 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TRFCON, ZEROT
      CHARACTER          DIST, NORM, TRANS, TYPE
      String             PATH;
      int                I, IMAT, IN, INFO, IRHS, ITRAN, IX, IZERO, J, K, KL, KOFF, KU, LDA, M, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN
      REAL               AINVNM, ANORM, COND, RCOND, RCONDC, RCONDI, RCONDO
*     ..
*     .. Local Arrays ..
      CHARACTER          TRANSS( 3 )
      int                ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS ), Z( 3 )
*     ..
*     .. External Functions ..
      REAL               SASUM, SGET06, SLANGT
      EXTERNAL           SASUM, SGET06, SLANGT
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, SCOPY, SERRGE, SGET04, SGTCON, SGTRFS, SGTT01, SGTT02, SGTT05, SGTTRF, SGTTRS, SLACPY, SLAGTM, SLARNV, SLATB4, SLATMS, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
*     .. Data statements ..
      DATA               ISEEDY / 0, 0, 0, 1 / , TRANSS / 'N', 'T', 'C' /
*     ..
*     .. Executable Statements ..
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'GT'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR ) CALL SERRGE( PATH, NOUT )
      INFOT = 0
*
      DO 110 IN = 1, NN
*
*        Do for each value of N in NVAL.
*
         N = NVAL( IN )
         M = MAX( N-1, 0 )
         LDA = MAX( 1, N )
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1
*
         DO 100 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) ) GO TO 100
*
*           Set up parameters with SLATB4.
*
            CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST )
*
            ZEROT = IMAT.GE.8 .AND. IMAT.LE.10
            IF( IMAT.LE.6 ) THEN
*
*              Types 1-6:  generate matrices of known condition number.
*
               KOFF = MAX( 2-KU, 3-MAX( 1, N ) )
               SRNAMT = 'SLATMS'
               CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'Z', AF( KOFF ), 3, WORK, INFO )
*
*              Check the error code from SLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'SLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 100
               END IF
               IZERO = 0
*
               IF( N.GT.1 ) THEN
                  CALL SCOPY( N-1, AF( 4 ), 3, A, 1 )
                  CALL SCOPY( N-1, AF( 3 ), 3, A( N+M+1 ), 1 )
               END IF
               CALL SCOPY( N, AF( 2 ), 3, A( M+1 ), 1 )
            ELSE
*
*              Types 7-12:  generate tridiagonal matrices with
*              unknown condition numbers.
*
               IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 7 ) ) THEN
*
*                 Generate a matrix with elements from [-1,1].
*
                  CALL SLARNV( 2, ISEED, N+2*M, A )
                  IF( ANORM.NE.ONE ) CALL SSCAL( N+2*M, ANORM, A, 1 )
               ELSE IF( IZERO.GT.0 ) THEN
*
*                 Reuse the last matrix by copying back the zeroed out
*                 elements.
*
                  IF( IZERO.EQ.1 ) THEN
                     A( N ) = Z( 2 )
                     IF( N.GT.1 ) A( 1 ) = Z( 3 )
                  ELSE IF( IZERO.EQ.N ) THEN
                     A( 3*N-2 ) = Z( 1 )
                     A( 2*N-1 ) = Z( 2 )
                  ELSE
                     A( 2*N-2+IZERO ) = Z( 1 )
                     A( N-1+IZERO ) = Z( 2 )
                     A( IZERO ) = Z( 3 )
                  END IF
               END IF
*
*              If IMAT > 7, set one column of the matrix to 0.
*
               IF( .NOT.ZEROT ) THEN
                  IZERO = 0
               ELSE IF( IMAT.EQ.8 ) THEN
                  IZERO = 1
                  Z( 2 ) = A( N )
                  A( N ) = ZERO
                  IF( N.GT.1 ) THEN
                     Z( 3 ) = A( 1 )
                     A( 1 ) = ZERO
                  END IF
               ELSE IF( IMAT.EQ.9 ) THEN
                  IZERO = N
                  Z( 1 ) = A( 3*N-2 )
                  Z( 2 ) = A( 2*N-1 )
                  A( 3*N-2 ) = ZERO
                  A( 2*N-1 ) = ZERO
               ELSE
                  IZERO = ( N+1 ) / 2
                  DO 20 I = IZERO, N - 1
                     A( 2*N-2+I ) = ZERO
                     A( N-1+I ) = ZERO
                     A( I ) = ZERO
   20             CONTINUE
                  A( 3*N-2 ) = ZERO
                  A( 2*N-1 ) = ZERO
               END IF
            END IF
*
*+    TEST 1
*           Factor A as L*U and compute the ratio
*              norm(L*U - A) / (n * norm(A) * EPS )
*
            CALL SCOPY( N+2*M, A, 1, AF, 1 )
            SRNAMT = 'SGTTRF'
            CALL SGTTRF( N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, INFO )
*
*           Check error code from SGTTRF.
*
            IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'SGTTRF', INFO, IZERO, ' ', N, N, 1, 1, -1, IMAT, NFAIL, NERRS, NOUT )
            TRFCON = INFO.NE.0
*
            CALL SGTT01( N, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, WORK, LDA, RWORK, RESULT( 1 ) )
*
*           Print the test ratio if it is .GE. THRESH.
*
            IF( RESULT( 1 ).GE.THRESH ) THEN
               IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 )
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
*
            DO 50 ITRAN = 1, 2
               TRANS = TRANSS( ITRAN )
               IF( ITRAN.EQ.1 ) THEN
                  NORM = 'O'
               ELSE
                  NORM = 'I'
               END IF
               ANORM = SLANGT( NORM, N, A, A( M+1 ), A( N+M+1 ) )
*
               IF( .NOT.TRFCON ) THEN
*
*                 Use SGTTRS to solve for one column at a time of inv(A)
*                 or inv(A^T), computing the maximum column sum as we
*                 go.
*
                  AINVNM = ZERO
                  DO 40 I = 1, N
                     DO 30 J = 1, N
                        X( J ) = ZERO
   30                CONTINUE
                     X( I ) = ONE
                     CALL SGTTRS( TRANS, N, 1, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO )
                     AINVNM = MAX( AINVNM, SASUM( N, X, 1 ) )
   40             CONTINUE
*
*                 Compute RCONDC = 1 / (norm(A) * norm(inv(A))
*
                  IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                     RCONDC = ONE
                  ELSE
                     RCONDC = ( ONE / ANORM ) / AINVNM
                  END IF
                  IF( ITRAN.EQ.1 ) THEN
                     RCONDO = RCONDC
                  ELSE
                     RCONDI = RCONDC
                  END IF
               ELSE
                  RCONDC = ZERO
               END IF
*
*+    TEST 7
*              Estimate the reciprocal of the condition number of the
*              matrix.
*
               SRNAMT = 'SGTCON'
               CALL SGTCON( NORM, N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, ANORM, RCOND, WORK, IWORK( N+1 ), INFO )
*
*              Check error code from SGTCON.
*
               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'SGTCON', INFO, 0, NORM, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
               RESULT( 7 ) = SGET06( RCOND, RCONDC )
*
*              Print the test ratio if it is .GE. THRESH.
*
               IF( RESULT( 7 ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                   WRITE( NOUT, FMT = 9997 )NORM, N, IMAT, 7, RESULT( 7 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1
   50       CONTINUE
*
*           Skip the remaining tests if the matrix is singular.
*
            IF( TRFCON ) GO TO 100
*
            DO 90 IRHS = 1, NNS
               NRHS = NSVAL( IRHS )
*
*              Generate NRHS random solution vectors.
*
               IX = 1
               DO 60 J = 1, NRHS
                  CALL SLARNV( 2, ISEED, N, XACT( IX ) )
                  IX = IX + LDA
   60          CONTINUE
*
               DO 80 ITRAN = 1, 3
                  TRANS = TRANSS( ITRAN )
                  IF( ITRAN.EQ.1 ) THEN
                     RCONDC = RCONDO
                  ELSE
                     RCONDC = RCONDI
                  END IF
*
*                 Set the right hand side.
*
                  CALL SLAGTM( TRANS, N, NRHS, ONE, A, A( M+1 ), A( N+M+1 ), XACT, LDA, ZERO, B, LDA )
*
*+    TEST 2
*                 Solve op(A) * X = B and compute the residual.
*
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
                  SRNAMT = 'SGTTRS'
                  CALL SGTTRS( TRANS, N, NRHS, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO )
*
*                 Check error code from SGTTRS.
*
                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'SGTTRS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                  CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                  CALL SGTT02( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), X, LDA, WORK, LDA, RESULT( 2 ) )
*
*+    TEST 3
*                 Check solution from generated exact solution.
*
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
*
*+    TESTS 4, 5, and 6
*                 Use iterative refinement to improve the solution.
*
                  SRNAMT = 'SGTRFS'
                  CALL SGTRFS( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO )
*
*                 Check error code from SGTRFS.
*
                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'SGTRFS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                  CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )                   CALL SGTT05( TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) )
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 70 K = 2, 6
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
   70             CONTINUE
                  NRUN = NRUN + 5
   80          CONTINUE
   90       CONTINUE
*
  100    CONTINUE
  110 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 12X, 'N =', I5, ',', 10X, ' type ', I2, ', test(', I2,
     $      ') = ', G12.5 )
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') = ', G12.5 )
 9997 FORMAT( ' NORM =''', A1, ''', N =', I5, ',', 10X, ' type ', I2,
     $      ', test(', I2, ') = ', G12.5 )
      RETURN
*
*     End of SCHKGT
*
      END
