      SUBROUTINE CDRVHE_RK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, E, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      REAL               RWORK( * )
      COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ), E( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 10, NTESTS = 3 ;
      int                NFACT;
      const              NFACT = 2 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, FACT, TYPE, UPLO, XTYPE;
      String             MATPATH, PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT;
      REAL               AINVNM, ANORM, CNDNUM, RCONDC
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS )

      // ..
      // .. External Functions ..
      REAL               CLANHE
      // EXTERNAL CLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, CERRVX, CGET04, CLACPY, CLARHS, CLATB4, CLATMS, CHESV_RK, CHET01_3, CPOT02, CHETRF_RK, CHETRI_3
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
      // INTRINSIC MAX, MIN
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      // Test path

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'HK'

      // Path to generate matrices

      MATPATH( 1: 1 ) = 'Complex precision'
      MATPATH( 2: 3 ) = 'HE'

      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10
      LWORK = MAX( 2*NMAX, NMAX*NRHS )

      // Test the error exits

      if (TSTERR) CALL CERRVX( PATH, NOUT );
      INFOT = 0

      // Set the block size and minimum block size for which the block
      // routine should be used, which will be later returned by ILAENV.

      NB = 1
      NBMIN = 2
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 180
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         if (N.LE.0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 170

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 170

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT.GE.3 && IMAT.LE.6
            if (ZEROT && N < IMAT-2) GO TO 170;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               UPLO = UPLOS( IUPLO )

                  // Begin generate the test matrix A.

                  // Set up parameters with CLATB4 for the matrix generator
                  // based on the type of matrix to be generated.

                  clatb4(MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                  // Generate a matrix with CLATMS.

                  SRNAMT = 'CLATMS'
                  clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

                  // Check error code from CLATMS and handle error.

                  if ( INFO != 0 ) {
                     alaerh(PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 160
                  }

                  // For types 3-6, zero one or more rows and columns of
                  // the matrix to test that INFO is returned correctly.

                  if ( ZEROT ) {
                     if ( IMAT == 3 ) {
                        IZERO = 1
                     } else if ( IMAT == 4 ) {
                        IZERO = N
                     } else {
                        IZERO = N / 2 + 1
                     }

                     if ( IMAT < 6 ) {

                        // Set row and column IZERO to zero.

                        if ( IUPLO == 1 ) {
                           IOFF = ( IZERO-1 )*LDA
                           for (I = 1; I <= IZERO - 1; I++) { // 20
                              A( IOFF+I ) = ZERO
                           } // 20
                           IOFF = IOFF + IZERO
                           for (I = IZERO; I <= N; I++) { // 30
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
                           } // 30
                        } else {
                           IOFF = IZERO
                           for (I = 1; I <= IZERO - 1; I++) { // 40
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
                           } // 40
                           IOFF = IOFF - IZERO
                           for (I = IZERO; I <= N; I++) { // 50
                              A( IOFF+I ) = ZERO
                           } // 50
                        }
                     } else {
                        if ( IUPLO == 1 ) {

                        // Set the first IZERO rows and columns to zero.

                           IOFF = 0
                           for (J = 1; J <= N; J++) { // 70
                              I2 = MIN( J, IZERO )
                              for (I = 1; I <= I2; I++) { // 60
                                 A( IOFF+I ) = ZERO
                              } // 60
                              IOFF = IOFF + LDA
                           } // 70
                        } else {

                        // Set the first IZERO rows and columns to zero.

                           IOFF = 0
                           for (J = 1; J <= N; J++) { // 90
                              I1 = MAX( J, IZERO )
                              for (I = I1; I <= N; I++) { // 80
                                 A( IOFF+I ) = ZERO
                              } // 80
                              IOFF = IOFF + LDA
                           } // 90
                        }
                     }
                  } else {
                     IZERO = 0
                  }

                  // End generate the test matrix A.


               for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 150

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT )

                  // Compute the condition number

                  if ( ZEROT ) {
                     if (IFACT == 1) GO TO 150;
                     RCONDC = ZERO

                  } else if ( IFACT == 1 ) {

                     // Compute the 1-norm of A.

                     ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )

                     // Factor the matrix A.

                     clacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     chetrf_rk(UPLO, N, AFAC, LDA, E, IWORK, WORK, LWORK, INFO );

                     // Compute inv(A) and take its norm.

                     clacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                     LWORK = (N+NB+1)*(NB+3)

                     // We need to compute the inverse to compute
                     // RCONDC that is used later in TEST3.

                     csytri_3(UPLO, N, AINV, LDA, E, IWORK, WORK, LWORK, INFO );
                     AINVNM = CLANHE( '1', UPLO, N, AINV, LDA, RWORK )

                     // Compute the 1-norm condition number of A.

                     if ( ANORM.LE.ZERO || AINVNM.LE.ZERO ) {
                        RCONDC = ONE
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     }
                  }

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'CLARHS'
                  clarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  XTYPE = 'C'

                  // --- Test CHESV_RK  ---

                  if ( IFACT == 2 ) {
                     clacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     clacpy('Full', N, NRHS, B, LDA, X, LDA );

                     // Factor the matrix and solve the system using
                     // CHESV_RK.

                     SRNAMT = 'CHESV_RK'
                     chesv_rk(UPLO, N, NRHS, AFAC, LDA, E, IWORK, X, LDA, WORK, LWORK, INFO );

                     // Adjust the expected value of INFO to account for
                     // pivoting.

                     K = IZERO
                     if ( K > 0 ) {
                        } // 100
                        if ( IWORK( K ) < 0 ) {
                           if ( IWORK( K ) != -K ) {
                              K = -IWORK( K )
                              GO TO 100
                           }
                        } else if ( IWORK( K ) != K ) {
                           K = IWORK( K )
                           GO TO 100
                        }
                     }

                     // Check error code from CHESV_RK and handle error.

                     if ( INFO != K ) {
                        alaerh(PATH, 'CHESV_RK', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 120
                     } else if ( INFO != 0 ) {
                        GO TO 120
                     }

*+    TEST 1      Reconstruct matrix from factors and compute
                  // residual.

                     chet01_3(UPLO, N, A, LDA, AFAC, LDA, E, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );

*+    TEST 2      Compute residual of the computed solution.

                     clacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     cpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

*+    TEST 3
                  // Check solution from generated exact solution.

                     cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                     NT = 3

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 1; K <= NT; K++) { // 110
                        if ( RESULT( K ).GE.THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'CHESV_RK', UPLO, N, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1
                        }
                     } // 110
                     NRUN = NRUN + NT
                     } // 120
                  }

               } // 150

            } // 160
         } // 170
      } // 180

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN

      // End of CDRVHE_RK

      }
