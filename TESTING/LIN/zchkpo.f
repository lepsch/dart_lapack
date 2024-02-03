      SUBROUTINE ZCHKPO( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, NOUT )

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
      COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0, 0.0 ) ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, IZERO, K, KL, KU, LDA, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      double             ANORM, CNDNUM, RCOND, RCONDC;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, ZLANHE;
      // EXTERNAL DGET06, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, XLAENV, ZERRPO, ZGET04, ZLACPY, ZLAIPD, ZLARHS, ZLATB4, ZLATMS, ZPOCON, ZPORFS, ZPOT01, ZPOT02, ZPOT03, ZPOT05, ZPOTRF, ZPOTRI, ZPOTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'PO'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      if (TSTERR) CALL ZERRPO( PATH, NOUT );
      INFOT = 0

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 120
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         if (N <= 0) NIMAT = 1;

         IZERO = 0
         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 110

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 110

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 5
            if (ZEROT && N < IMAT-2) GO TO 110;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 100
               UPLO = UPLOS( IUPLO )

               // Set up parameters with ZLATB4 and generate a test matrix
               // with ZLATMS.

               zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'ZLATMS'
               zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from ZLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100
               }

               // For types 3-5, zero one row and column of the matrix to
               // test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT == 3 ) {
                     IZERO = 1
                  } else if ( IMAT == 4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }
                  IOFF = ( IZERO-1 )*LDA

                  // Set row and column IZERO of A to 0.

                  if ( IUPLO == 1 ) {
                     for (I = 1; I <= IZERO - 1; I++) { // 20
                        A( IOFF+I ) = CZERO
                     } // 20
                     IOFF = IOFF + IZERO
                     for (I = IZERO; I <= N; I++) { // 30
                        A( IOFF ) = CZERO
                        IOFF = IOFF + LDA
                     } // 30
                  } else {
                     IOFF = IZERO
                     for (I = 1; I <= IZERO - 1; I++) { // 40
                        A( IOFF ) = CZERO
                        IOFF = IOFF + LDA
                     } // 40
                     IOFF = IOFF - IZERO
                     for (I = IZERO; I <= N; I++) { // 50
                        A( IOFF+I ) = CZERO
                     } // 50
                  }
               } else {
                  IZERO = 0
               }

               // Set the imaginary part of the diagonals.

               zlaipd(N, A, LDA+1, 0 );

               // Do for each value of NB in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 90
                  NB = NBVAL( INB )
                  xlaenv(1, NB );

                  // Compute the L*L' or U'*U factorization of the matrix.

                  zlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                  SRNAMT = 'ZPOTRF'
                  zpotrf(UPLO, N, AFAC, LDA, INFO );

                  // Check error code from ZPOTRF.

                  if ( INFO != IZERO ) {
                     alaerh(PATH, 'ZPOTRF', INFO, IZERO, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 90
                  }

                  // Skip the tests if INFO is not 0.

                  if (INFO != 0) GO TO 90;

*+    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  zlacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                  zpot01(UPLO, N, A, LDA, AINV, LDA, RWORK, RESULT( 1 ) );

*+    TEST 2
                  // Form the inverse and compute the residual.

                  zlacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                  SRNAMT = 'ZPOTRI'
                  zpotri(UPLO, N, AINV, LDA, INFO );

                  // Check error code from ZPOTRI.

                  if (INFO != 0) CALL ALAERH( PATH, 'ZPOTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  zpot03(UPLO, N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= 2; K++) { // 60
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1
                     }
                  } // 60
                  NRUN = NRUN + 2

                  // Skip the rest of the tests unless this is the first
                  // blocksize.

                  if (INB != 1) GO TO 90;

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 80
                     NRHS = NSVAL( IRHS )

*+    TEST 3
                  // Solve and compute residual for A * X = B .

                     SRNAMT = 'ZLARHS'
                     zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                     SRNAMT = 'ZPOTRS'
                     zpotrs(UPLO, N, NRHS, AFAC, LDA, X, LDA, INFO );

                  // Check error code from ZPOTRS.

                     if (INFO != 0) CALL ALAERH( PATH, 'ZPOTRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     zpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

*+    TEST 4
                  // Check solution from generated exact solution.

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );

*+    TESTS 5, 6, and 7
                  // Use iterative refinement to improve the solution.

                     SRNAMT = 'ZPORFS'
                     zporfs(UPLO, N, NRHS, A, LDA, AFAC, LDA, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                  // Check error code from ZPORFS.

                     if (INFO != 0) CALL ALAERH( PATH, 'ZPORFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) );
                     zpot05(UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 6 ) );

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 3; K <= 7; K++) { // 70
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1
                        }
                     } // 70
                     NRUN = NRUN + 5
                  } // 80

*+    TEST 8
                  // Get an estimate of RCOND = 1/CNDNUM.

                  ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
                  SRNAMT = 'ZPOCON'
                  zpocon(UPLO, N, AFAC, LDA, ANORM, RCOND, WORK, RWORK, INFO );

                  // Check error code from ZPOCON.

                  if (INFO != 0) CALL ALAERH( PATH, 'ZPOCON', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  RESULT( 8 ) = DGET06( RCOND, RCONDC )

                  // Print the test ratio if it is >= THRESH.

                  if ( RESULT( 8 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 8, RESULT( 8 );
                     NFAIL = NFAIL + 1
                  }
                  NRUN = NRUN + 1
               } // 90
            } // 100
         } // 110
      } // 120

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ',', 10X, ' type ', I2, ', test(', I2, ') =', G12.5 )
      RETURN

      // End of ZCHKPO

      }
