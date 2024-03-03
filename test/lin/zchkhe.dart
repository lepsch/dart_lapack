      void zchkhe(final int DOTYPE, final int NN, final int NVAL, final int NNB, final int NBVAL, final int NNS, final int NSVAL, final int THRESH, final int TSTERR, final int NMAX, final int A, final int AFAC, final int AINV, final int B, final int X, final int XACT, final Array<double> _WORK_, final Array<double> RWORK_, final Array<int> IWORK_, final int NOUT,) {
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      Complex         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      int                NTYPES;
      const              NTYPES = 10 ;
      int                NTESTS;
      const              NTESTS = 9 ;
      bool               TRFCON, ZEROT;
      String             DIST, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      double             ANORM, CNDNUM, RCOND, RCONDC;
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- double             DGET06, ZLANHE;
      // EXTERNAL DGET06, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, XLAENV, ZERRHE, ZGET04, ZHECON, ZHERFS, ZHET01, ZHETRF, ZHETRI2, ZHETRS, ZLACPY, ZLAIPD, ZLARHS, ZLATB4, ZLATMS, ZPOT02, ZPOT03, ZPOT05
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'HE';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) zerrhe( PATH, NOUT );
      infoc.INFOT = 0;

      // Set the minimum block size for which the block routine should
      // be used, which will be later returned by ILAENV

      xlaenv(2, 2 );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 180
         N = NVAL( IN );
         LDA = max( N, 1 );
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         IZERO = 0;
         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 170

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 170;

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 6;
            if (ZEROT && N < IMAT-2) GO TO 170;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               UPLO = UPLOS( IUPLO );

               // Set up parameters with ZLATB4 for the matrix generator
               // based on the type of matrix to be generated.

               zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               // Generate a matrix with ZLATMS.

              srnamc.SRNAMT = 'ZLATMS';
               zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from ZLATMS and handle error.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  // Skip all tests for this generated matrix

                  GO TO 160;
               }

               // For types 3-6, zero one or more rows and columns of
               // the matrix to test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT == 3 ) {
                     IZERO = 1;
                  } else if ( IMAT == 4 ) {
                     IZERO = N;
                  } else {
                     IZERO = N / 2 + 1;
                  }

                  if ( IMAT < 6 ) {

                     // Set row and column IZERO to zero.

                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*LDA;
                        for (I = 1; I <= IZERO - 1; I++) { // 20
                           A[IOFF+I] = CZERO;
                        } // 20
                        IOFF = IOFF + IZERO;
                        for (I = IZERO; I <= N; I++) { // 30
                           A[IOFF] = CZERO;
                           IOFF = IOFF + LDA;
                        } // 30
                     } else {
                        IOFF = IZERO;
                        for (I = 1; I <= IZERO - 1; I++) { // 40
                           A[IOFF] = CZERO;
                           IOFF = IOFF + LDA;
                        } // 40
                        IOFF = IOFF - IZERO;
                        for (I = IZERO; I <= N; I++) { // 50
                           A[IOFF+I] = CZERO;
                        } // 50
                     }
                  } else {
                     if ( IUPLO == 1 ) {

                        // Set the first IZERO rows and columns to zero.

                        IOFF = 0;
                        for (J = 1; J <= N; J++) { // 70
                           I2 = min( J, IZERO );
                           for (I = 1; I <= I2; I++) { // 60
                              A[IOFF+I] = CZERO;
                           } // 60
                           IOFF = IOFF + LDA;
                        } // 70
                     } else {

                        // Set the last IZERO rows and columns to zero.

                        IOFF = 0;
                        for (J = 1; J <= N; J++) { // 90
                           I1 = max( J, IZERO );
                           for (I = I1; I <= N; I++) { // 80
                              A[IOFF+I] = CZERO;
                           } // 80
                           IOFF = IOFF + LDA;
                        } // 90
                     }
                  }
               } else {
                  IZERO = 0;
               }

               // End generate test matrix A.


               // Set the imaginary part of the diagonals.

               zlaipd(N, A, LDA+1, 0 );

               // Do for each value of NB in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 150

                  // Set the optimal blocksize, which will be later
                  // returned by ILAENV.

                  NB = NBVAL( INB );
                  xlaenv(1, NB );

                  // Copy the test matrix A into matrix AFAC which
                  // will be factorized in place. This is needed to
                  // preserve the test matrix A for subsequent tests.

                  zlacpy(UPLO, N, N, A, LDA, AFAC, LDA );

                  // Compute the L*D*L**T or U*D*U**T factorization of the
                  // matrix. IWORK stores details of the interchanges and
                  // the block structure of D. AINV is a work array for
                  // block factorization, LWORK is the length of AINV.

                  LWORK = max( 2, NB )*LDA;
                 srnamc.SRNAMT = 'ZHETRF';
                  zhetrf(UPLO, N, AFAC, LDA, IWORK, AINV, LWORK, INFO );

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                  K = IZERO;
                  if ( K > 0 ) {
                     } // 100
                     if ( IWORK( K ) < 0 ) {
                        if ( IWORK( K ) != -K ) {
                           K = -IWORK( K );
                           GO TO 100;
                        }
                     } else if ( IWORK( K ) != K ) {
                        K = IWORK( K );
                        GO TO 100;
                     }
                  }

                  // Check error code from ZHETRF and handle error.

                  if (INFO != K) alaerh( PATH, 'ZHETRF', INFO, K, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                  // Set the condition estimate flag if the INFO is not 0.

                  if ( INFO != 0 ) {
                     TRFCON = true;
                  } else {
                     TRFCON = false;
                  }

// +    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  zhet01(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );
                  NT = 1;

// +    TEST 2
                  // Form the inverse and compute the residual.

                  if ( INB == 1 && !TRFCON ) {
                     zlacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                    srnamc.SRNAMT = 'ZHETRI2';
                     LWORK = (N+NB+1)*(NB+3);
                     zhetri2(UPLO, N, AINV, LDA, IWORK, WORK, LWORK, INFO );

                     // Check error code from ZHETRI and handle error.

                     if (INFO != 0) alaerh( PATH, 'ZHETRI', INFO, -1, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     // Compute the residual for a symmetric matrix times
                     // its inverse.

                     zpot03(UPLO, N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );
                     NT = 2;
                  }

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= NT; K++) { // 110
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 110
                  NRUN = NRUN + NT;

                  // Skip the other tests if this is not the first block
                  // size.

                  if (INB > 1) GO TO 150;

                  // Do only the condition estimate if INFO is not 0.

                  if ( TRFCON ) {
                     RCONDC = ZERO;
                     GO TO 140;
                  }

                  // Do for each value of NRHS in NSVAL.

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 130
                     NRHS = NSVAL( IRHS );

// +    TEST 3 (Using TRS)
                  // Solve and compute residual for  A * X = B.

                     // Choose a set of NRHS random solution vectors
                     // stored in XACT and set up the right hand side B

                    srnamc.SRNAMT = 'ZLARHS';
                     zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                    srnamc.SRNAMT = 'ZHETRS';
                     zhetrs(UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO );

                     // Check error code from ZHETRS and handle error.

                     if (INFO != 0) alaerh( PATH, 'ZHETRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );

                     // Compute the residual for the solution

                     zpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

// +    TEST 4 (Using TRS2)
                  // Solve and compute residual for  A * X = B.

                     // Choose a set of NRHS random solution vectors
                     // stored in XACT and set up the right hand side B

                    srnamc.SRNAMT = 'ZLARHS';
                     zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                    srnamc.SRNAMT = 'ZHETRS2';
                     zhetrs2(UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, WORK, INFO );

                     // Check error code from ZHETRS2 and handle error.

                     if (INFO != 0) alaerh( PATH, 'ZHETRS2', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );

                     // Compute the residual for the solution

                     zpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 4 ) );

// +    TEST 5
                  // Check solution from generated exact solution.

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) );

// +    TESTS 6, 7, and 8
                  // Use iterative refinement to improve the solution.

                    srnamc.SRNAMT = 'ZHERFS';
                     zherfs(UPLO, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                     // Check error code from ZHERFS.

                     if (INFO != 0) alaerh( PATH, 'ZHERFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 6 ) );
                     zpot05(UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 7 ) );

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 3; K <= 8; K++) { // 120
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                           WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1;
                        }
                     } // 120
                     NRUN = NRUN + 6;

                  // End do for each value of NRHS in NSVAL.

                  } // 130

// +    TEST 9
                  // Get an estimate of RCOND = 1/CNDNUM.

                  } // 140
                  ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK );
                 srnamc.SRNAMT = 'ZHECON';
                  zhecon(UPLO, N, AFAC, LDA, IWORK, ANORM, RCOND, WORK, INFO );

                  // Check error code from ZHECON and handle error.

                  if (INFO != 0) alaerh( PATH, 'ZHECON', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  RESULT[9] = DGET06( RCOND, RCONDC );

                  // Print information about the tests that did not pass
                  // the threshold.

                  if ( RESULT( 9 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 9, RESULT( 9 );
                     NFAIL = NFAIL + 1;
                  }
                  NRUN = NRUN + 1;
               } // 150
            } // 160
         } // 170
      } // 180

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = \'${.a1}\', N =${.i5}, NB =${.i4}, type ${.i2}, test ${.i2}, ratio =${.g12_5};
 9998 FORMAT( ' UPLO = \'${.a1}\', N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) =${.g12_5};
 9997 FORMAT( ' UPLO = \'${.a1}\', N =${.i5},${' ' * 10} type ${.i2}, test(${.i2}) =${.g12_5};
      }
