      void zdrvsy(final Array<bool> DOTYPE_, final int NN, final Array<int> NVAL_, final int NRHS, final double THRESH, final bool TSTERR, final int NMAX, final Array<double> A_, final Array<double> AFAC_, final Array<double> AINV_, final Array<double> B_, final Array<double> X_, final Array<double> XACT_, final Array<double> WORK_, final Array<double> RWORK_, final Array<int> IWORK_, final Nout NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      double             RWORK( * );
      Complex         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 11, NTESTS = 6 ;
      int                NFACT;
      const              NFACT = 2 ;
      bool               ZEROT;
      String             DIST, FACT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, K1, KL, KU, LDA, LWORK, MODE, N, NIMAT, NT;
      double             AINVNM, ANORM, CNDNUM, RCOND, RCONDC;
      String             FACTS( NFACT ), UPLOS( 2 );
      final                ISEED=Array<int>( 4 );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Functions ..
      //- double             DGET06, ZLANSY;
      // EXTERNAL DGET06, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZERRVX, ZGET04, ZLACPY, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZLATSY, ZPOT05, ZSYSV, ZSYSVX, ZSYT01, ZSYT02, ZSYTRF, ZSYTRI2
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
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', FACTS = 'F', 'N';

      // Initialize constants and the random number seed.

      final PATH = '${'Zomplex precision'[0]}SY';
      var NRUN = 0;
      var NFAIL = 0;
      var NERRS = Box(0);
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY[I - 1];
      } // 10
      LWORK = max( 2*NMAX, NMAX*NRHS );

      // Test the error exits

      if (TSTERR) zerrvx( PATH, NOUT );
      infoc.INFOT = 0;

      // Set the block size and minimum block size for testing.

            final NB = 1;
      final NBMIN = 2;
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 180
         final N = NVAL[IN];
         final LDA = max( N, 1 );
         XTYPE = 'N';
            final NIMAT = N <= 0 ? 1 : NTYPES;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 170

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE[IMAT] ) GO TO 170;

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            final ZEROT = IMAT >= 3 && IMAT <= 6;
            if (ZEROT && N < IMAT-2) GO TO 170;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               final UPLO = UPLOS[IUPLO - 1];

               if ( IMAT != NTYPES ) {

                  // Set up parameters with ZLATB4 and generate a test
                  // matrix with ZLATMS.

                  final (:TYPE,:KL,:KU,:ANORM,:MODE,:CNDNUM,:DIST) = zlatb4(PATH, IMAT, N, N);

                 srnamc.SRNAMT = 'ZLATMS';
                  zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

                  // Check error code from ZLATMS.

                  if ( INFO.value != 0 ) {
                     alaerh(PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 160;
                  }

                  // For types 3-6, zero one or more rows and columns of
                  // the matrix to test that INFO is returned correctly.

                  final int IZERO;
                  if ( ZEROT ) {
                     if ( IMAT == 3 ) {
                        IZERO = 1;
                     } else if ( IMAT == 4 ) {
                        IZERO = N;
                     } else {
                        IZERO = N ~/ 2 + 1;
                     }

                     if ( IMAT < 6 ) {

                        // Set row and column IZERO to zero.

                        if ( IUPLO == 1 ) {
                           IOFF = ( IZERO-1 )*LDA;
                           for (I = 1; I <= IZERO - 1; I++) { // 20
                              A[IOFF+I] = ZERO;
                           } // 20
                           IOFF = IOFF + IZERO;
                           for (I = IZERO; I <= N; I++) { // 30
                              A[IOFF] = ZERO;
                              IOFF = IOFF + LDA;
                           } // 30
                        } else {
                           IOFF = IZERO;
                           for (I = 1; I <= IZERO - 1; I++) { // 40
                              A[IOFF] = ZERO;
                              IOFF = IOFF + LDA;
                           } // 40
                           IOFF = IOFF - IZERO;
                           for (I = IZERO; I <= N; I++) { // 50
                              A[IOFF+I] = ZERO;
                           } // 50
                        }
                     } else {
                        if ( IUPLO == 1 ) {

                           // Set the first IZERO rows to zero.

                           IOFF = 0;
                           for (J = 1; J <= N; J++) { // 70
                              I2 = min( J, IZERO );
                              for (I = 1; I <= I2; I++) { // 60
                                 A[IOFF+I] = ZERO;
                              } // 60
                              IOFF = IOFF + LDA;
                           } // 70
                        } else {

                           // Set the last IZERO rows to zero.

                           IOFF = 0;
                           for (J = 1; J <= N; J++) { // 90
                              I1 = max( J, IZERO );
                              for (I = I1; I <= N; I++) { // 80
                                 A[IOFF+I] = ZERO;
                              } // 80
                              IOFF = IOFF + LDA;
                           } // 90
                        }
                     }
                  } else {
                     IZERO = 0;
                  }
               } else {

                  // IMAT = NTYPES:  Use a special block diagonal matrix to
                  // test alternate code for the 2-by-2 blocks.

                  zlatsy(UPLO, N, A, LDA, ISEED );
               }

               for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 150

                  // Do first for FACT = 'F', then for other values.

                  final FACT = FACTS[IFACT - 1];

                  // Compute the condition number for comparison with
                  // the value returned by ZSYSVX.

                  final int IZERO;
                  if ( ZEROT ) {
                     if (IFACT == 1) GO TO 150;
                     RCONDC = ZERO;

                  } else if ( IFACT == 1 ) {

                     // Compute the 1-norm of A.

                     ANORM = ZLANSY( '1', UPLO, N, A, LDA, RWORK );

                     // Factor the matrix A.

                     zlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     zsytrf(UPLO, N, AFAC, LDA, IWORK, WORK, LWORK, INFO );

                     // Compute inv(A) and take its norm.

                     zlacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                     LWORK = (N+NB+1)*(NB+3);
                     zsytri2(UPLO, N, AINV, LDA, IWORK, WORK, LWORK, INFO );
                     AINVNM = ZLANSY( '1', UPLO, N, AINV, LDA, RWORK );

                     // Compute the 1-norm condition number of A.

                     if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                        RCONDC = ONE;
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM;
                     }
                  }

                  // Form an exact solution and set the right hand side.

                 srnamc.SRNAMT = 'ZLARHS';
                  zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  XTYPE = 'C';

                  // --- Test ZSYSV  ---

                  if ( IFACT == 2 ) {
                     zlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                     // Factor the matrix and solve the system using ZSYSV.

                    srnamc.SRNAMT = 'ZSYSV ';
                     zsysv(UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, WORK, LWORK, INFO );

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

                     // Check error code from ZSYSV .

                     if ( INFO.value != K ) {
                        alaerh(PATH, 'ZSYSV ', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 120;
                     } else if ( INFO.value != 0 ) {
                        GO TO 120;
                     }

                     // Reconstruct matrix from factors and compute
                     // residual.

                     zsyt01(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );

                     // Compute residual of the computed solution.

                     zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     zsyt02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                     NT = 3;

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 1; K <= NT; K++) { // 110
                        if ( RESULT[K] >= THRESH ) {
                           if (NFAIL == 0 && NERRS.value == 0) aladhd( NOUT, PATH );
                           NOUT.println( 9999 )'ZSYSV ', UPLO, N, IMAT, K, RESULT[K];
                           NFAIL++;
                        }
                     } // 110
                     NRUN +=  NT;
                     } // 120
                  }

                  // --- Test ZSYSVX ---

                  if (IFACT == 2) zlaset( UPLO, N, N, Complex.zero, Complex.zero, AFAC, LDA );
                  zlaset('Full', N, NRHS, Complex.zero, Complex.zero, X, LDA );

                  // Solve the system and compute the condition number and
                  // error bounds using ZSYSVX.

                 srnamc.SRNAMT = 'ZSYSVX';
                  zsysvx(FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, LWORK, RWORK( 2*NRHS+1 ), INFO );

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                  K = IZERO;
                  if ( K > 0 ) {
                     } // 130
                     if ( IWORK( K ) < 0 ) {
                        if ( IWORK( K ) != -K ) {
                           K = -IWORK( K );
                           GO TO 130;
                        }
                     } else if ( IWORK( K ) != K ) {
                        K = IWORK( K );
                        GO TO 130;
                     }
                  }

                  // Check the error code from ZSYSVX.

                  if ( INFO.value != K ) {
                     alaerh(PATH, 'ZSYSVX', INFO, K, FACT + UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 150;
                  }

                  if ( INFO.value == 0 ) {
                     if ( IFACT >= 2 ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        zsyt01(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1;
                     } else {
                        K1 = 2;
                     }

                     // Compute residual of the computed solution.

                     zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     zsyt02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

                     // Check the error bounds from iterative refinement.

                     zpot05(UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                  } else {
                     K1 = 6;
                  }

                  // Compare RCOND from ZSYSVX with the computed value
                  // in RCONDC.

                  RESULT[6] = DGET06( RCOND, RCONDC );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = K1; K <= 6; K++) { // 140
                     if ( RESULT[K] >= THRESH ) {
                        if (NFAIL == 0 && NERRS.value == 0) aladhd( NOUT, PATH );
                        NOUT.println( 9998 )'ZSYSVX', FACT, UPLO, N, IMAT, K, RESULT[K];
                        NFAIL++;
                     }
                  } // 140
                  NRUN +=  7 - K1;

               } // 150

            } // 160
         } // 170
      } // 180

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT(' ${}, UPLO=\'${.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${.i2}, ratio =${RESULT[].g12_5};
 9998 FORMAT(' ${}, FACT=\'${FACT.a1}\', UPLO=\'${.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${.i2}, ratio =${RESULT[].g12_5};
      }
