      SUBROUTINE DCHKTR( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

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
      int                IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double             A( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
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
      int                I, IDIAG, IMAT, IN, INB, INFO, IRHS, ITRAN, IUPLO, K, LDA, N, NB, NERRS, NFAIL, NRHS, NRUN;
      double             AINVNM, ANORM, BIGNUM, DLAMCH, DUMMY, RCOND, RCONDC, RCONDI, RCONDO, RES, SCALE;
      // ..
      // .. Local Arrays ..
      String             TRANSS( NTRAN ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS ), SCALE3( 2 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLANTR;
      // EXTERNAL LSAME, DLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DCOPY, DERRTR, DGET04, DLACPY, DLAMCH, DSCAL, DLARHS, DLATRS, DLATRS3, DLATTR, DTRCON, DTRRFS, DTRT01, DTRT02, DTRT03, DTRT05, DTRT06, DTRTRI, DTRTRS, XLAENV
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, IOUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
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
      PATH( 2: 3 ) = 'TR'
      BIGNUM = DLAMCH('Overflow') / DLAMCH('Precision')
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      if (TSTERR) CALL DERRTR( PATH, NOUT );
      INFOT = 0
      xlaenv(2, 2 );

      for (IN = 1; IN <= NN; IN++) { // 120

         // Do for each value of N in NVAL

         N = NVAL( IN )
         LDA = MAX( 1, N )
         XTYPE = 'N'

         for (IMAT = 1; IMAT <= NTYPE1; IMAT++) { // 80

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 80

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 70

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO )

               // Call DLATTR to generate a triangular test matrix.

               SRNAMT = 'DLATTR'
               dlattr(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, A, LDA, X, WORK, INFO );

               // Set IDIAG = 1 for non-unit matrices, 2 for unit.

               if ( LSAME( DIAG, 'N' ) ) {
                  IDIAG = 1
               } else {
                  IDIAG = 2
               }

               for (INB = 1; INB <= NNB; INB++) { // 60

                  // Do for each blocksize in NBVAL

                  NB = NBVAL( INB )
                  xlaenv(1, NB );

*+    TEST 1
                  // Form the inverse of A.

                  dlacpy(UPLO, N, N, A, LDA, AINV, LDA );
                  SRNAMT = 'DTRTRI'
                  dtrtri(UPLO, DIAG, N, AINV, LDA, INFO );

                  // Check error code from DTRTRI.

                  if (INFO != 0) CALL ALAERH( PATH, 'DTRTRI', INFO, 0, UPLO // DIAG, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                  // Compute the infinity-norm condition number of A.

                  ANORM = DLANTR( 'I', UPLO, DIAG, N, N, A, LDA, RWORK )
                  AINVNM = DLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, RWORK )
                  if ( ANORM.LE.ZERO || AINVNM.LE.ZERO ) {
                     RCONDI = ONE
                  } else {
                     RCONDI = ( ONE / ANORM ) / AINVNM
                  }

                  // Compute the residual for the triangular matrix times
                  // its inverse.  Also compute the 1-norm condition number
                  // of A.

                  dtrt01(UPLO, DIAG, N, A, LDA, AINV, LDA, RCONDO, RWORK, RESULT( 1 ) );

                  // Print the test ratio if it is >= THRESH.

                  if ( RESULT( 1 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9999 )UPLO, DIAG, N, NB, IMAT, 1, RESULT( 1 );
                     NFAIL = NFAIL + 1
                  }
                  NRUN = NRUN + 1

                  // Skip remaining tests if not the first block size.

                  if (INB != 1) GO TO 60;

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 40
                     NRHS = NSVAL( IRHS )
                     XTYPE = 'N'

                     for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 30

                     // Do for op(A) = A, A**T, or A**H.

                        TRANS = TRANSS( ITRAN )
                        if ( ITRAN == 1 ) {
                           NORM = 'O'
                           RCONDC = RCONDO
                        } else {
                           NORM = 'I'
                           RCONDC = RCONDI
                        }

*+    TEST 2
                        // Solve and compute residual for op(A)*x = b.

                        SRNAMT = 'DLARHS'
                        dlarhs(PATH, XTYPE, UPLO, TRANS, N, N, 0, IDIAG, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                        XTYPE = 'C'
                        dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                        SRNAMT = 'DTRTRS'
                        dtrtrs(UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDA, INFO );

                        // Check error code from DTRTRS.

                        if (INFO != 0) CALL ALAERH( PATH, 'DTRTRS', INFO, 0, UPLO // TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        // This line is needed on a Sun SPARCstation.

                        if (N > 0) DUMMY = A( 1 );

                        dtrt02(UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDA, B, LDA, WORK, RESULT( 2 ) );

*+    TEST 3
                        // Check solution from generated exact solution.

                        dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

*+    TESTS 4, 5, and 6
                        // Use iterative refinement to improve the solution
                        // and compute error bounds.

                        SRNAMT = 'DTRRFS'
                        dtrrfs(UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO );

                        // Check error code from DTRRFS.

                        if (INFO != 0) CALL ALAERH( PATH, 'DTRRFS', INFO, 0, UPLO // TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                        dtrt05(UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 2; K <= 6; K++) { // 20
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )UPLO, TRANS, DIAG, N, NRHS, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1
                           }
                        } // 20
                        NRUN = NRUN + 5
                     } // 30
                  } // 40

*+    TEST 7
                        // Get an estimate of RCOND = 1/CNDNUM.

                  for (ITRAN = 1; ITRAN <= 2; ITRAN++) { // 50
                     if ( ITRAN == 1 ) {
                        NORM = 'O'
                        RCONDC = RCONDO
                     } else {
                        NORM = 'I'
                        RCONDC = RCONDI
                     }
                     SRNAMT = 'DTRCON'
                     dtrcon(NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, IWORK, INFO );

                        // Check error code from DTRCON.

                     if (INFO != 0) CALL ALAERH( PATH, 'DTRCON', INFO, 0, NORM // UPLO // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     dtrt06(RCOND, RCONDC, UPLO, DIAG, N, A, LDA, RWORK, RESULT( 7 ) );

                     // Print the test ratio if it is >= THRESH.

                     if ( RESULT( 7 ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )NORM, UPLO, N, IMAT, 7, RESULT( 7 );
                        NFAIL = NFAIL + 1
                     }
                     NRUN = NRUN + 1
                  } // 50
               } // 60
            } // 70
         } // 80

         // Use pathological test matrices to test DLATRS.

         for (IMAT = NTYPE1 + 1; IMAT <= NTYPES; IMAT++) { // 110

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 110

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 100

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO )
               for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 90

                  // Do for op(A) = A, A**T, and A**H.

                  TRANS = TRANSS( ITRAN )

                  // Call DLATTR to generate a triangular test matrix.

                  SRNAMT = 'DLATTR'
                  dlattr(IMAT, UPLO, TRANS, DIAG, ISEED, N, A, LDA, X, WORK, INFO );

*+    TEST 8
                  // Solve the system op(A)*x = b.

                  SRNAMT = 'DLATRS'
                  dcopy(N, X, 1, B, 1 );
                  dlatrs(UPLO, TRANS, DIAG, 'N', N, A, LDA, B, SCALE, RWORK, INFO );

                  // Check error code from DLATRS.

                  if (INFO != 0) CALL ALAERH( PATH, 'DLATRS', INFO, 0, UPLO // TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  dtrt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) );

*+    TEST 9
                  // Solve op(A)*X = b again with NORMIN = 'Y'.

                  dcopy(N, X, 1, B( N+1 ), 1 );
                  dlatrs(UPLO, TRANS, DIAG, 'Y', N, A, LDA, B( N+1 ), SCALE, RWORK, INFO );

                  // Check error code from DLATRS.

                  if (INFO != 0) CALL ALAERH( PATH, 'DLATRS', INFO, 0, UPLO // TRANS // DIAG // 'Y', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  dtrt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RESULT( 9 ) );

*+    TEST 10
                  // Solve op(A)*X = B

                  SRNAMT = 'DLATRS3'
                  dcopy(N, X, 1, B, 1 );
                  dcopy(N, X, 1, B( N+1 ), 1 );
                  dscal(N, BIGNUM, B( N+1 ), 1 );
                  dlatrs3(UPLO, TRANS, DIAG, 'N', N, 2, A, LDA, B, MAX(1, N), SCALE3, RWORK, WORK, NMAX, INFO );

                  // Check error code from DLATRS3.

                  if (INFO != 0) CALL ALAERH( PATH, 'DLATRS3', INFO, 0, UPLO // TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  dtrt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE3( 1 ), RWORK, ONE, B( 1 ), LDA, X, LDA, WORK, RESULT( 10 ) );
                  dscal(N, BIGNUM, X, 1 );
                  dtrt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE3( 2 ), RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RES );
                  RESULT( 10 ) = MAX( RESULT( 10 ), RES )

                  // Print information about the tests that did not pass
                  // the threshold.

                  if ( RESULT( 8 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9996 )'DLATRS', UPLO, TRANS, DIAG, 'N', N, IMAT, 8, RESULT( 8 );
                     NFAIL = NFAIL + 1
                  }
                  if ( RESULT( 9 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9996 )'DLATRS', UPLO, TRANS, DIAG, 'Y', N, IMAT, 9, RESULT( 9 );
                     NFAIL = NFAIL + 1
                  }
                  if ( RESULT( 10 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9996 )'DLATRS3', UPLO, TRANS, DIAG, 'N', N, IMAT, 10, RESULT( 10 );
                     NFAIL = NFAIL + 1
                  }
                  NRUN = NRUN + 3
               } // 90
            } // 100
         } // 110
      } // 120

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO=''', A1, ''', DIAG=''', A1, ''', N=', I5, ', NB=', I4, ', type ', I2, ', test(', I2, ')= ', G12.5 )
 9998 FORMAT( ' UPLO=''', A1, ''', TRANS=''', A1, ''', DIAG=''', A1, ''', N=', I5, ', NB=', I4, ', type ', I2, ', test(', I2, ')= ', G12.5 )
 9997 FORMAT( ' NORM=''', A1, ''', UPLO =''', A1, ''', N=', I5, ',', 11X, ' type ', I2, ', test(', I2, ')=', G12.5 )
 9996 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ''', A1, ''', ''', A1, ''',', I5, ', ... ), type ', I2, ', test(', I2, ')=', G12.5 )
      RETURN

      // End of DCHKTR

      }
