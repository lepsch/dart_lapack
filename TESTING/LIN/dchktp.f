      SUBROUTINE DCHKTP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, AP, AINVP, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NSVAL( * ), NVAL( * );
      double             AINVP( * ), AP( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTYPE1, NTYPES;
      const              NTYPE1 = 10, NTYPES = 18 ;
      int                NTESTS;
      const              NTESTS = 9 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      String             DIAG, NORM, TRANS, UPLO, XTYPE;
      String             PATH;
      int                I, IDIAG, IMAT, IN, INFO, IRHS, ITRAN, IUPLO, K, LAP, LDA, N, NERRS, NFAIL, NRHS, NRUN;
      double             AINVNM, ANORM, RCOND, RCONDC, RCONDI, RCONDO, SCALE;
      // ..
      // .. Local Arrays ..
      String             TRANSS( NTRAN ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLANTP;
      // EXTERNAL LSAME, DLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DCOPY, DERRTR, DGET04, DLACPY, DLARHS, DLATPS, DLATTP, DTPCON, DTPRFS, DTPT01, DTPT02, DTPT03, DTPT05, DTPT06, DTPTRI, DTPTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, IOUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , TRANSS / 'N', 'T', 'C' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'TP'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL DERRTR( PATH, NOUT )
      INFOT = 0

      DO 110 IN = 1, NN

         // Do for each value of N in NVAL

         N = NVAL( IN )
         LDA = MAX( 1, N )
         LAP = LDA*( LDA+1 ) / 2
         XTYPE = 'N'

         DO 70 IMAT = 1, NTYPE1

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 70

            DO 60 IUPLO = 1, 2

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO )

               // Call DLATTP to generate a triangular test matrix.

               SRNAMT = 'DLATTP'
               CALL DLATTP( IMAT, UPLO, 'No transpose', DIAG, ISEED, N, AP, X, WORK, INFO )

               // Set IDIAG = 1 for non-unit matrices, 2 for unit.

               IF( LSAME( DIAG, 'N' ) ) THEN
                  IDIAG = 1
               ELSE
                  IDIAG = 2
               END IF

*+    TEST 1
               // Form the inverse of A.

               IF( N.GT.0 ) CALL DCOPY( LAP, AP, 1, AINVP, 1 )
               SRNAMT = 'DTPTRI'
               CALL DTPTRI( UPLO, DIAG, N, AINVP, INFO )

               // Check error code from DTPTRI.

               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DTPTRI', INFO, 0, UPLO // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

               // Compute the infinity-norm condition number of A.

               ANORM = DLANTP( 'I', UPLO, DIAG, N, AP, RWORK )
               AINVNM = DLANTP( 'I', UPLO, DIAG, N, AINVP, RWORK )
               IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                  RCONDI = ONE
               ELSE
                  RCONDI = ( ONE / ANORM ) / AINVNM
               END IF

               // Compute the residual for the triangular matrix times its
               // inverse.  Also compute the 1-norm condition number of A.

               CALL DTPT01( UPLO, DIAG, N, AP, AINVP, RCONDO, RWORK, RESULT( 1 ) )

               // Print the test ratio if it is .GE. THRESH.

               IF( RESULT( 1 ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                   WRITE( NOUT, FMT = 9999 )UPLO, DIAG, N, IMAT, 1, RESULT( 1 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1

               DO 40 IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
                  XTYPE = 'N'

                  DO 30 ITRAN = 1, NTRAN

                  // Do for op(A) = A, A**T, or A**H.

                     TRANS = TRANSS( ITRAN )
                     IF( ITRAN.EQ.1 ) THEN
                        NORM = 'O'
                        RCONDC = RCONDO
                     ELSE
                        NORM = 'I'
                        RCONDC = RCONDI
                     END IF

*+    TEST 2
                  // Solve and compute residual for op(A)*x = b.

                     SRNAMT = 'DLARHS'
                     CALL DLARHS( PATH, XTYPE, UPLO, TRANS, N, N, 0, IDIAG, NRHS, AP, LAP, XACT, LDA, B, LDA, ISEED, INFO )
                     XTYPE = 'C'
                     CALL DLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                     SRNAMT = 'DTPTRS'
                     CALL DTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, X, LDA, INFO )

                  // Check error code from DTPTRS.

                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DTPTRS', INFO, 0, UPLO // TRANS // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                     CALL DTPT02( UPLO, TRANS, DIAG, N, NRHS, AP, X, LDA, B, LDA, WORK, RESULT( 2 ) )

*+    TEST 3
                  // Check solution from generated exact solution.

                     CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )

*+    TESTS 4, 5, and 6
                  // Use iterative refinement to improve the solution and
                  // compute error bounds.

                     SRNAMT = 'DTPRFS'
                     CALL DTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO )

                  // Check error code from DTPRFS.

                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DTPRFS', INFO, 0, UPLO // TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                     CALL DGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )                      CALL DTPT05( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) )

                     // Print information about the tests that did not pass
                    t // he threshold.

                     DO 20 K = 2, 6
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9998 )UPLO, TRANS, DIAG, N, NRHS, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
   20                CONTINUE
                     NRUN = NRUN + 5
   30             CONTINUE
   40          CONTINUE

*+    TEST 7
                  // Get an estimate of RCOND = 1/CNDNUM.

               DO 50 ITRAN = 1, 2
                  IF( ITRAN.EQ.1 ) THEN
                     NORM = 'O'
                     RCONDC = RCONDO
                  ELSE
                     NORM = 'I'
                     RCONDC = RCONDI
                  END IF

                  SRNAMT = 'DTPCON'
                  CALL DTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK, INFO )

                  // Check error code from DTPCON.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DTPCON', INFO, 0, NORM // UPLO // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                  CALL DTPT06( RCOND, RCONDC, UPLO, DIAG, N, AP, RWORK, RESULT( 7 ) )

                  // Print the test ratio if it is .GE. THRESH.

                  IF( RESULT( 7 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9997 ) 'DTPCON', NORM, UPLO, DIAG, N, IMAT, 7, RESULT( 7 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE

         // Use pathological test matrices to test DLATPS.

         DO 100 IMAT = NTYPE1 + 1, NTYPES

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 100

            DO 90 IUPLO = 1, 2

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO )
               DO 80 ITRAN = 1, NTRAN

                  // Do for op(A) = A, A**T, or A**H.

                  TRANS = TRANSS( ITRAN )

                  // Call DLATTP to generate a triangular test matrix.

                  SRNAMT = 'DLATTP'
                  CALL DLATTP( IMAT, UPLO, TRANS, DIAG, ISEED, N, AP, X, WORK, INFO )

*+    TEST 8
                  // Solve the system op(A)*x = b.

                  SRNAMT = 'DLATPS'
                  CALL DCOPY( N, X, 1, B, 1 )
                  CALL DLATPS( UPLO, TRANS, DIAG, 'N', N, AP, B, SCALE, RWORK, INFO )

                  // Check error code from DLATPS.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DLATPS', INFO, 0, UPLO // TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                  CALL DTPT03( UPLO, TRANS, DIAG, N, 1, AP, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) )

*+    TEST 9
                  // Solve op(A)*x = b again with NORMIN = 'Y'.

                  CALL DCOPY( N, X, 1, B( N+1 ), 1 )
                  CALL DLATPS( UPLO, TRANS, DIAG, 'Y', N, AP, B( N+1 ), SCALE, RWORK, INFO )

                  // Check error code from DLATPS.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DLATPS', INFO, 0, UPLO // TRANS // DIAG // 'Y', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                  CALL DTPT03( UPLO, TRANS, DIAG, N, 1, AP, SCALE, RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RESULT( 9 ) )

                  // Print information about the tests that did not pass
                 t // he threshold.

                  IF( RESULT( 8 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9996 )'DLATPS', UPLO, TRANS, DIAG, 'N', N, IMAT, 8, RESULT( 8 )
                     NFAIL = NFAIL + 1
                  END IF
                  IF( RESULT( 9 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9996 )'DLATPS', UPLO, TRANS, DIAG, 'Y', N, IMAT, 9, RESULT( 9 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 2
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( ' UPLO=''', A1, ''', DIAG=''', A1, ''', N=', I5,
     $      ', type ', I2, ', test(', I2, ')= ', G12.5 )
 9998 FORMAT( ' UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1,
     $      ''', N=', I5, ''', NRHS=', I5, ', type ', I2, ', test(',
     $      I2, ')= ', G12.5 )
 9997 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''',',
     $      I5, ', ... ), type ', I2, ', test(', I2, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''', ''',
     $      A1, ''',', I5, ', ... ), type ', I2, ', test(', I2, ')=',
     $      G12.5 )
      RETURN

      // End of DCHKTP

      }
