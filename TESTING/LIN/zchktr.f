      SUBROUTINE ZCHKTR( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AINV, B, X, XACT, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NBVAL( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      COMPLEX*16         A( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTYPE1, NTYPES;
      const              NTYPE1 = 10, NTYPES = 18 ;
      int                NTESTS;
      const              NTESTS = 10 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      double             ONE, ZERO;
      const              ONE = 1.0D0, ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      String             DIAG, NORM, TRANS, UPLO, XTYPE;
      String             PATH;
      int                I, IDIAG, IMAT, IN, INB, INFO, IRHS, ITRAN, IUPLO, K, LDA, N, NB, NERRS, NFAIL, NRHS, NRUN       double             AINVNM, ANORM, BIGNUM, DUMMY, RCOND, RCONDC, RCONDI, RCONDO, RES, SCALE, DLAMCH;;
      // ..
      // .. Local Arrays ..
      String             TRANSS( NTRAN ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS ), RWORK2( 2*NMAX ), SCALE3( 2 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             ZLANTR;
      // EXTERNAL LSAME, ZLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DLAMCH, XLAENV, ZCOPY, ZDSCAL, ZERRTR, ZGET04, ZLACPY, ZLARHS, ZLATRS, ZLATRS3, ZLATTR, ZTRCON, ZTRRFS, ZTRT01, ZTRT02, ZTRT03, ZTRT05, ZTRT06, ZTRTRI, ZTRTRS
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

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'TR'
      BIGNUM = DLAMCH('Overflow') / DLAMCH('Precision')
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL ZERRTR( PATH, NOUT )
      INFOT = 0

      DO 120 IN = 1, NN

         // Do for each value of N in NVAL

         N = NVAL( IN )
         LDA = MAX( 1, N )
         XTYPE = 'N'

         DO 80 IMAT = 1, NTYPE1

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 80

            DO 70 IUPLO = 1, 2

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO )

               // Call ZLATTR to generate a triangular test matrix.

               SRNAMT = 'ZLATTR'
               CALL ZLATTR( IMAT, UPLO, 'No transpose', DIAG, ISEED, N, A, LDA, X, WORK, RWORK, INFO )

               // Set IDIAG = 1 for non-unit matrices, 2 for unit.

               IF( LSAME( DIAG, 'N' ) ) THEN
                  IDIAG = 1
               ELSE
                  IDIAG = 2
               END IF

               DO 60 INB = 1, NNB

                  // Do for each blocksize in NBVAL

                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )

*+    TEST 1
                  // Form the inverse of A.

                  CALL ZLACPY( UPLO, N, N, A, LDA, AINV, LDA )
                  SRNAMT = 'ZTRTRI'
                  CALL ZTRTRI( UPLO, DIAG, N, AINV, LDA, INFO )

                  // Check error code from ZTRTRI.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZTRTRI', INFO, 0, UPLO // DIAG, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT )

                  // Compute the infinity-norm condition number of A.

                  ANORM = ZLANTR( 'I', UPLO, DIAG, N, N, A, LDA, RWORK )
                  AINVNM = ZLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, RWORK )
                  IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
                     RCONDI = ONE
                  ELSE
                     RCONDI = ( ONE / ANORM ) / AINVNM
                  END IF

                  // Compute the residual for the triangular matrix times
                  // its inverse.  Also compute the 1-norm condition number
                  // of A.

                  CALL ZTRT01( UPLO, DIAG, N, A, LDA, AINV, LDA, RCONDO, RWORK, RESULT( 1 ) )
                  // Print the test ratio if it is .GE. THRESH.

                  IF( RESULT( 1 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9999 )UPLO, DIAG, N, NB, IMAT, 1, RESULT( 1 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1

                  // Skip remaining tests if not the first block size.

                  IF( INB.NE.1 ) GO TO 60

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

                        SRNAMT = 'ZLARHS'
                        CALL ZLARHS( PATH, XTYPE, UPLO, TRANS, N, N, 0, IDIAG, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                        XTYPE = 'C'
                        CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                        SRNAMT = 'ZTRTRS'
                        CALL ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDA, INFO )

                        // Check error code from ZTRTRS.

                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZTRTRS', INFO, 0, UPLO // TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                        // This line is needed on a Sun SPARCstation.

                        IF( N.GT.0 ) DUMMY = DBLE( A( 1 ) )

                        CALL ZTRT02( UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDA, B, LDA, WORK, RWORK, RESULT( 2 ) )

*+    TEST 3
                        // Check solution from generated exact solution.

                        CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )

*+    TESTS 4, 5, and 6
                        // Use iterative refinement to improve the solution
                        // and compute error bounds.

                        SRNAMT = 'ZTRRFS'
                        CALL ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )

                        // Check error code from ZTRRFS.

                        IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZTRRFS', INFO, 0, UPLO // TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                        CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )                         CALL ZTRT05( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) )

                        // Print information about the tests that did not
                        // pass the threshold.

                        DO 20 K = 2, 6
                           IF( RESULT( K ).GE.THRESH ) THEN
                              IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )UPLO, TRANS, DIAG, N, NRHS, IMAT, K, RESULT( K )
                              NFAIL = NFAIL + 1
                           END IF
   20                   CONTINUE
                        NRUN = NRUN + 5
   30                CONTINUE
   40             CONTINUE

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
                     SRNAMT = 'ZTRCON'
                     CALL ZTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, RWORK, INFO )

                        // Check error code from ZTRCON.

                     IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZTRCON', INFO, 0, NORM // UPLO // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                     CALL ZTRT06( RCOND, RCONDC, UPLO, DIAG, N, A, LDA, RWORK, RESULT( 7 ) )

                     // Print the test ratio if it is .GE. THRESH.

                     IF( RESULT( 7 ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )NORM, UPLO, N, IMAT, 7, RESULT( 7 )
                        NFAIL = NFAIL + 1
                     END IF
                     NRUN = NRUN + 1
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE

         // Use pathological test matrices to test ZLATRS.

         DO 110 IMAT = NTYPE1 + 1, NTYPES

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 110

            DO 100 IUPLO = 1, 2

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO )
               DO 90 ITRAN = 1, NTRAN

                  // Do for op(A) = A, A**T, and A**H.

                  TRANS = TRANSS( ITRAN )

                  // Call ZLATTR to generate a triangular test matrix.

                  SRNAMT = 'ZLATTR'
                  CALL ZLATTR( IMAT, UPLO, TRANS, DIAG, ISEED, N, A, LDA, X, WORK, RWORK, INFO )

*+    TEST 8
                  // Solve the system op(A)*x = b.

                  SRNAMT = 'ZLATRS'
                  CALL ZCOPY( N, X, 1, B, 1 )
                  CALL ZLATRS( UPLO, TRANS, DIAG, 'N', N, A, LDA, B, SCALE, RWORK, INFO )

                  // Check error code from ZLATRS.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZLATRS', INFO, 0, UPLO // TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                  CALL ZTRT03( UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) )

*+    TEST 9
                  // Solve op(A)*X = b again with NORMIN = 'Y'.

                  CALL ZCOPY( N, X, 1, B( N+1 ), 1 )
                  CALL ZLATRS( UPLO, TRANS, DIAG, 'Y', N, A, LDA, B( N+1 ), SCALE, RWORK, INFO )

                  // Check error code from ZLATRS.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZLATRS', INFO, 0, UPLO // TRANS // DIAG // 'Y', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                  CALL ZTRT03( UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RESULT( 9 ) )

*+    TEST 10
                  // Solve op(A)*X = B

                  SRNAMT = 'ZLATRS3'
                  CALL ZCOPY( N, X, 1, B, 1 )
                  CALL ZCOPY( N, X, 1, B( N+1 ), 1 )
                  CALL ZDSCAL( N, BIGNUM, B( N+1 ), 1 )
                  CALL ZLATRS3( UPLO, TRANS, DIAG, 'N', N, 2, A, LDA, B, MAX(1, N), SCALE3, RWORK, RWORK2, 2*NMAX, INFO )

                  // Check error code from ZLATRS3.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'ZLATRS3', INFO, 0, UPLO // TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  CALL ZTRT03( UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE3( 1 ), RWORK, ONE, B( 1 ), LDA, X, LDA, WORK, RESULT( 10 ) )
                  CALL ZDSCAL( N, BIGNUM, X, 1 )
                  CALL ZTRT03( UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE3( 2 ), RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RES )
                  RESULT( 10 ) = MAX( RESULT( 10 ), RES )

                  // Print information about the tests that did not pass
                 t // he threshold.

                  IF( RESULT( 8 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9996 )'ZLATRS', UPLO, TRANS, DIAG, 'N', N, IMAT, 8, RESULT( 8 )
                     NFAIL = NFAIL + 1
                  END IF
                  IF( RESULT( 9 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9996 )'ZLATRS', UPLO, TRANS, DIAG, 'Y', N, IMAT, 9, RESULT( 9 )
                     NFAIL = NFAIL + 1
                  END IF
                  IF( RESULT( 10 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9996 )'ZLATRS3', UPLO, TRANS, DIAG, 'N', N, IMAT, 10, RESULT( 10 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 3
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( ' UPLO=''', A1, ''', DIAG=''', A1, ''', N=', I5, ', NB=',
     $      I4, ', type ', I2, ', test(', I2, ')= ', G12.5 )
 9998 FORMAT( ' UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1,
     $      ''', N=', I5, ', NB=', I4, ', type ', I2, ', test(',
     $      I2, ')= ', G12.5 )
 9997 FORMAT( ' NORM=''', A1, ''', UPLO =''', A1, ''', N=', I5, ',',
     $      11X, ' type ', I2, ', test(', I2, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''', ''',
     $      A1, ''',', I5, ', ... ), type ', I2, ', test(', I2, ')=',
     $      G12.5 )
      RETURN

      // End of ZCHKTR

      }
