      SUBROUTINE CCHKPO( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NBVAL( * ), NSVAL( * ), NVAL( * );
      REAL               RWORK( * )
      COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CZERO
      const              CZERO = ( 0.0E+0, 0.0E+0 ) ;
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
      REAL               ANORM, CNDNUM, RCOND, RCONDC
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Functions ..
      REAL               CLANHE, SGET06
      // EXTERNAL CLANHE, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CERRPO, CGET04, CLACPY, CLAIPD, CLARHS, CLATB4, CLATMS, CPOCON, CPORFS, CPOT01, CPOT02, CPOT03, CPOT05, CPOTRF, CPOTRI, CPOTRS, XLAENV
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

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'PO'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      if (TSTERR) CALL CERRPO( PATH, NOUT );
      INFOT = 0

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 120
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         if (N.LE.0) NIMAT = 1;

         IZERO = 0
         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 110

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 110

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            if (ZEROT .AND. N.LT.IMAT-2) GO TO 110;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 100
               UPLO = UPLOS( IUPLO )

               // Set up parameters with CLATB4 and generate a test matrix
               // with CLATMS.

               clatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'CLATMS'
               clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from CLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
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

               claipd(N, A, LDA+1, 0 );

               // Do for each value of NB in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 90
                  NB = NBVAL( INB )
                  xlaenv(1, NB );

                  // Compute the L*L' or U'*U factorization of the matrix.

                  clacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                  SRNAMT = 'CPOTRF'
                  cpotrf(UPLO, N, AFAC, LDA, INFO );

                  // Check error code from CPOTRF.

                  if ( INFO != IZERO ) {
                     alaerh(PATH, 'CPOTRF', INFO, IZERO, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 90
                  }

                  // Skip the tests if INFO is not 0.

                  if (INFO != 0) GO TO 90;

*+    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  clacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                  cpot01(UPLO, N, A, LDA, AINV, LDA, RWORK, RESULT( 1 ) );

*+    TEST 2
                  // Form the inverse and compute the residual.

                  clacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                  SRNAMT = 'CPOTRI'
                  cpotri(UPLO, N, AINV, LDA, INFO );

                  // Check error code from CPOTRI.

                  if (INFO != 0) CALL ALAERH( PATH, 'CPOTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  cpot03(UPLO, N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= 2; K++) { // 60
                     if ( RESULT( K ).GE.THRESH ) {
                        if (NFAIL == 0 .AND. NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K );
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

                     SRNAMT = 'CLARHS'
                     clarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     clacpy('Full', N, NRHS, B, LDA, X, LDA );

                     SRNAMT = 'CPOTRS'
                     cpotrs(UPLO, N, NRHS, AFAC, LDA, X, LDA, INFO );

                  // Check error code from CPOTRS.

                     if (INFO != 0) CALL ALAERH( PATH, 'CPOTRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     clacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     cpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

*+    TEST 4
                  // Check solution from generated exact solution.

                     cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );

*+    TESTS 5, 6, and 7
                  // Use iterative refinement to improve the solution.

                     SRNAMT = 'CPORFS'
                     cporfs(UPLO, N, NRHS, A, LDA, AFAC, LDA, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                  // Check error code from CPORFS.

                     if (INFO != 0) CALL ALAERH( PATH, 'CPORFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) );
                     cpot05(UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 6 ) );

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 3; K <= 7; K++) { // 70
                        if ( RESULT( K ).GE.THRESH ) {
                           if (NFAIL == 0 .AND. NERRS == 0) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1
                        }
                     } // 70
                     NRUN = NRUN + 5
                  } // 80

*+    TEST 8
                  // Get an estimate of RCOND = 1/CNDNUM.

                  ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
                  SRNAMT = 'CPOCON'
                  cpocon(UPLO, N, AFAC, LDA, ANORM, RCOND, WORK, RWORK, INFO );

                  // Check error code from CPOCON.

                  if (INFO != 0) CALL ALAERH( PATH, 'CPOCON', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  RESULT( 8 ) = SGET06( RCOND, RCONDC )

                  // Print the test ratio if it is .GE. THRESH.

                  if ( RESULT( 8 ).GE.THRESH ) {
                     if (NFAIL == 0 .AND. NERRS == 0) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 8, RESULT( 8 );
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

      // End of CCHKPO

      }
