      SUBROUTINE DCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * )       double             A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );;
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      int                NTYPES;
      PARAMETER          ( NTYPES = 11 )
      int                NTESTS;
      PARAMETER          ( NTESTS = 8 )
      int                NTRAN;
      PARAMETER          ( NTRAN = 3 )
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, NORM, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IM, IMAT, IN, INB, INFO, IOFF, IRHS, ITRAN, IZERO, K, KL, KU, LDA, LWORK, M, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      double             AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, DUMMY, RCOND, RCONDC, RCONDI, RCONDO;
      // ..
      // .. Local Arrays ..
      String             TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, DLANGE;
      // EXTERNAL DGET06, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DERRGE, DGECON, DGERFS, DGET01, DGET02, DGET03, DGET04, DGET07, DGETRF, DGETRI, DGETRS, DLACPY, DLARHS, DLASET, DLATB4, DLATMS, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 / , TRANSS / 'N', 'T', 'C' /
      // ..
      // .. Executable Statements ..
*
      // Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'GE'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
      // Test the error exits
*
      CALL XLAENV( 1, 1 )
      IF( TSTERR ) CALL DERRGE( PATH, NOUT )
      INFOT = 0
      CALL XLAENV( 2, 2 )
*
      // Do for each value of M in MVAL
*
      DO 120 IM = 1, NM
         M = MVAL( IM )
         LDA = MAX( 1, M )
*
         // Do for each value of N in NVAL
*
         DO 110 IN = 1, NN
            N = NVAL( IN )
            XTYPE = 'N'
            NIMAT = NTYPES
            IF( M.LE.0 .OR. N.LE.0 ) NIMAT = 1
*
            DO 100 IMAT = 1, NIMAT
*
               // Do the tests only if DOTYPE( IMAT ) is true.
*
               IF( .NOT.DOTYPE( IMAT ) ) GO TO 100
*
               // Skip types 5, 6, or 7 if the matrix size is too small.
*
               ZEROT = IMAT.GE.5 .AND. IMAT.LE.7
               IF( ZEROT .AND. N.LT.IMAT-4 ) GO TO 100
*
               // Set up parameters with DLATB4 and generate a test matrix
               // with DLATMS.
*
               CALL DLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
               SRNAMT = 'DLATMS'
               CALL DLATMS( M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO )
*
               // Check error code from DLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'DLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 100
               END IF
*
               // For types 5-7, zero one or more columns of the matrix to
              t // est that INFO is returned correctly.
*
               IF( ZEROT ) THEN
                  IF( IMAT.EQ.5 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT.EQ.6 ) THEN
                     IZERO = MIN( M, N )
                  ELSE
                     IZERO = MIN( M, N ) / 2 + 1
                  END IF
                  IOFF = ( IZERO-1 )*LDA
                  IF( IMAT.LT.7 ) THEN
                     DO 20 I = 1, M
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                  ELSE
                     CALL DLASET( 'Full', M, N-IZERO+1, ZERO, ZERO, A( IOFF+1 ), LDA )
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
               // These lines, if used in place of the calls in the DO 60
               // loop, cause the code to bomb on a Sun SPARCstation.
*
                // ANORMO = DLANGE( 'O', M, N, A, LDA, RWORK )
                // ANORMI = DLANGE( 'I', M, N, A, LDA, RWORK )
*
               // Do for each blocksize in NBVAL
*
               DO 90 INB = 1, NNB
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
*
                  // Compute the LU factorization of the matrix.
*
                  CALL DLACPY( 'Full', M, N, A, LDA, AFAC, LDA )
                  SRNAMT = 'DGETRF'
                  CALL DGETRF( M, N, AFAC, LDA, IWORK, INFO )
*
                  // Check error code from DGETRF.
*
                  IF( INFO.NE.IZERO ) CALL ALAERH( PATH, 'DGETRF', INFO, IZERO, ' ', M, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )
                  TRFCON = .FALSE.
*
*+    TEST 1
                  // Reconstruct matrix from factors and compute residual.
*
                  CALL DLACPY( 'Full', M, N, AFAC, LDA, AINV, LDA )
                  CALL DGET01( M, N, A, LDA, AINV, LDA, IWORK, RWORK, RESULT( 1 ) )
                  NT = 1
*
*+    TEST 2
                  // Form the inverse if the factorization was successful
                  // and compute the residual.
*
                  IF( M.EQ.N .AND. INFO.EQ.0 ) THEN
                     CALL DLACPY( 'Full', N, N, AFAC, LDA, AINV, LDA )
                     SRNAMT = 'DGETRI'
                     NRHS = NSVAL( 1 )
                     LWORK = NMAX*MAX( 3, NRHS )
                     CALL DGETRI( N, AINV, LDA, IWORK, WORK, LWORK, INFO )
*
                     // Check error code from DGETRI.
*
                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGETRI', INFO, 0, ' ', N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )
*
                     // Compute the residual for the matrix times its
                     // inverse.  Also compute the 1-norm condition number
                     // of A.
*
                     CALL DGET03( N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDO, RESULT( 2 ) )
                     ANORMO = DLANGE( 'O', M, N, A, LDA, RWORK )
*
                     // Compute the infinity-norm condition number of A.
*
                     ANORMI = DLANGE( 'I', M, N, A, LDA, RWORK )
                     AINVNM = DLANGE( 'I', N, N, AINV, LDA, RWORK )
                     IF( ANORMI.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                        RCONDI = ONE
                     ELSE
                        RCONDI = ( ONE / ANORMI ) / AINVNM
                     END IF
                     NT = 2
                  ELSE
*
                     // Do only the condition estimate if INFO > 0.
*
                     TRFCON = .TRUE.
                     ANORMO = DLANGE( 'O', M, N, A, LDA, RWORK )
                     ANORMI = DLANGE( 'I', M, N, A, LDA, RWORK )
                     RCONDO = ZERO
                     RCONDI = ZERO
                  END IF
*
                  // Print information about the tests so far that did not
                  // pass the threshold.
*
                  DO 30 K = 1, NT
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )M, N, NB, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
   30             CONTINUE
                  NRUN = NRUN + NT
*
                  // Skip the remaining tests if this is not the first
                  // block size or if M .ne. N.  Skip the solve tests if
                 t // he matrix is singular.
*
                  IF( INB.GT.1 .OR. M.NE.N ) GO TO 90                   IF( TRFCON ) GO TO 70
*
                  DO 60 IRHS = 1, NNS
                     NRHS = NSVAL( IRHS )
                     XTYPE = 'N'
*
                     DO 50 ITRAN = 1, NTRAN
                        TRANS = TRANSS( ITRAN )
                        IF( ITRAN.EQ.1 ) THEN
                           RCONDC = RCONDO
                        ELSE
                           RCONDC = RCONDI
                        END IF
*
*+    TEST 3
                        // Solve and compute residual for A * X = B.
*
                        SRNAMT = 'DLARHS'
                        CALL DLARHS( PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                        XTYPE = 'C'
*
                        CALL DLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
                        SRNAMT = 'DGETRS'
                        CALL DGETRS( TRANS, N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO )
*
                        // Check error code from DGETRS.
*
                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGETRS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                        CALL DLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )                         CALL DGET02( TRANS, N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) )
*
*+    TEST 4
                        // Check solution from generated exact solution.
*
                        CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )
*
*+    TESTS 5, 6, and 7
                        // Use iterative refinement to improve the
                        // solution.
*
                        SRNAMT = 'DGERFS'
                        CALL DGERFS( TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO )
*
                        // Check error code from DGERFS.
*
                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGERFS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
                        CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) )                         CALL DGET07( TRANS, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, .TRUE., RWORK( NRHS+1 ), RESULT( 6 ) )
*
                        // Print information about the tests that did not
                        // pass the threshold.
*
                        DO 40 K = 3, 7
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   40                   CONTINUE
                        NRUN = NRUN + 5
   50                CONTINUE
   60             CONTINUE
*
*+    TEST 8
                     // Get an estimate of RCOND = 1/CNDNUM.
*
   70             CONTINUE
                  DO 80 ITRAN = 1, 2
                     IF( ITRAN.EQ.1 ) THEN
                        ANORM = ANORMO
                        RCONDC = RCONDO
                        NORM = 'O'
                     ELSE
                        ANORM = ANORMI
                        RCONDC = RCONDI
                        NORM = 'I'
                     END IF
                     SRNAMT = 'DGECON'
                     CALL DGECON( NORM, N, AFAC, LDA, ANORM, RCOND, WORK, IWORK( N+1 ), INFO )
*
                        // Check error code from DGECON.
*
                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DGECON', INFO, 0, NORM, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
                        // This line is needed on a Sun SPARCstation.
*
                     DUMMY = RCOND
*
                     RESULT( 8 ) = DGET06( RCOND, RCONDC )
*
                     // Print information about the tests that did not pass
                    t // he threshold.
*
                     IF( RESULT( 8 ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )NORM, N, IMAT, 8, RESULT( 8 )
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
*
      // Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' M = ', I5, ', N =', I5, ', NB =', I4, ', type ', I2,
     $      ', test(', I2, ') =', G12.5 )
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' NORM =''', A1, ''', N =', I5, ',', 10X, ' type ', I2,
     $      ', test(', I2, ') =', G12.5 )
      RETURN
*
      // End of DCHKGE
*
      END
