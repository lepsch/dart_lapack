      void zchksp(final Array<bool> DOTYPE_, final int NN, final Array<int> NVAL_, final int NNS, final Array<int> NSVAL_, final double THRESH, final bool TSTERR, final int NMAX, final Array<double> A_, final Array<double> AFAC_, final Array<double> AINV_, final Array<double> B_, final Array<double> X_, final Array<double> XACT_, final Array<double> WORK_, final Array<double> RWORK_, final Array<int> IWORK_, final Nout NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      Complex         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      bool               TRFCON, ZEROT;
      String             DIST, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IMAT, IN, INFO, IOFF, IRHS, IUPLO, IZERO, J, K, KL, KU, LDA, MODE, N, NIMAT, NPP, NRHS, NT;
      double             ANORM, CNDNUM, RCOND, RCONDC;
      String             UPLOS( 2 );
      final                ISEED=Array<int>( 4 );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DGET06, ZLANSP;
      // EXTERNAL lsame, DGET06, ZLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, ZCOPY, ZERRSY, ZGET04, ZLACPY, ZLARHS, ZLATB4, ZLATMS, ZLATSP, ZPPT05, ZSPCON, ZSPRFS, ZSPT01, ZSPT02, ZSPT03, ZSPTRF, ZSPTRI, ZSPTRS
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

      final PATH = '${'Zomplex precision'[0]}SP';
      var NRUN = 0;
      var NFAIL = 0;
      var NERRS = Box(0);
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY[I - 1];
      } // 10

      // Test the error exits

      if (TSTERR) zerrsy( PATH, NOUT );
      infoc.INFOT = 0;

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 170
         final N = NVAL[IN];
         final LDA = max( N, 1 );
         XTYPE = 'N';
            final NIMAT = N <= 0 ? 1 : NTYPES;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 160

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE[IMAT] ) GO TO 160;

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            final ZEROT = IMAT >= 3 && IMAT <= 6;
            if (ZEROT && N < IMAT-2) GO TO 160;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 150
               final UPLO = UPLOS[IUPLO - 1];
               if ( lsame( UPLO, 'U' ) ) {
                  PACKIT = 'C';
               } else {
                  PACKIT = 'R';
               }

               if ( IMAT != NTYPES ) {

                  // Set up parameters with ZLATB4 and generate a test
                  // matrix with ZLATMS.

                  final (:TYPE,:KL,:KU,:ANORM,:MODE,:CNDNUM,:DIST) = zlatb4(PATH, IMAT, N, N);

                 srnamc.SRNAMT = 'ZLATMS';
                  zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO );

                  // Check error code from ZLATMS.

                  if ( INFO.value != 0 ) {
                     alaerh(PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 150;
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
                           IOFF = ( IZERO-1 )*IZERO / 2;
                           for (I = 1; I <= IZERO - 1; I++) { // 20
                              A[IOFF+I] = ZERO;
                           } // 20
                           IOFF = IOFF + IZERO;
                           for (I = IZERO; I <= N; I++) { // 30
                              A[IOFF] = ZERO;
                              IOFF = IOFF + I;
                           } // 30
                        } else {
                           IOFF = IZERO;
                           for (I = 1; I <= IZERO - 1; I++) { // 40
                              A[IOFF] = ZERO;
                              IOFF = IOFF + N - I;
                           } // 40
                           IOFF = IOFF - IZERO;
                           for (I = IZERO; I <= N; I++) { // 50
                              A[IOFF+I] = ZERO;
                           } // 50
                        }
                     } else {
                        if ( IUPLO == 1 ) {

                           // Set the first IZERO rows and columns to zero.

                           IOFF = 0;
                           for (J = 1; J <= N; J++) { // 70
                              I2 = min( J, IZERO );
                              for (I = 1; I <= I2; I++) { // 60
                                 A[IOFF+I] = ZERO;
                              } // 60
                              IOFF = IOFF + J;
                           } // 70
                        } else {

                           // Set the last IZERO rows and columns to zero.

                           IOFF = 0;
                           for (J = 1; J <= N; J++) { // 90
                              I1 = max( J, IZERO );
                              for (I = I1; I <= N; I++) { // 80
                                 A[IOFF+I] = ZERO;
                              } // 80
                              IOFF = IOFF + N - J;
                           } // 90
                        }
                     }
                  } else {
                     IZERO = 0;
                  }
               } else {

                  // Use a special block diagonal matrix to test alternate
                  // code for the 2 x 2 blocks.

                  zlatsp(UPLO, N, A, ISEED );
               }

               // Compute the L*D*L' or U*D*U' factorization of the matrix.

               NPP = N*( N+1 ) / 2;
               zcopy(NPP, A, 1, AFAC, 1 );
              srnamc.SRNAMT = 'ZSPTRF';
               zsptrf(UPLO, N, AFAC, IWORK, INFO );

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

               // Check error code from ZSPTRF.

               if (INFO != K) alaerh( PATH, 'ZSPTRF', INFO, K, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               if ( INFO.value != 0 ) {
                  TRFCON = true;
               } else {
                  TRFCON = false;
               }

// +    TEST 1
               // Reconstruct matrix from factors and compute residual.

               zspt01(UPLO, N, A, AFAC, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );
               NT = 1;

// +    TEST 2
               // Form the inverse and compute the residual.

               if ( !TRFCON ) {
                  zcopy(NPP, AFAC, 1, AINV, 1 );
                 srnamc.SRNAMT = 'ZSPTRI';
                  zsptri(UPLO, N, AINV, IWORK, WORK, INFO );

               // Check error code from ZSPTRI.

                  if (INFO != 0) alaerh( PATH, 'ZSPTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  zspt03(UPLO, N, A, AINV, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );
                  NT = 2;
               }

               // Print information about the tests that did not pass
               // the threshold.

               for (K = 1; K <= NT; K++) { // 110
                  if ( RESULT[K] >= THRESH ) {
                     if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                     NOUT.println( 9999 )UPLO, N, IMAT, K, RESULT[K];
                     NFAIL++;
                  }
               } // 110
               NRUN +=  NT;

               // Do only the condition estimate if INFO is not 0.

               if ( TRFCON ) {
                  RCONDC = ZERO;
                  GO TO 140;
               }

               for (IRHS = 1; IRHS <= NNS; IRHS++) { // 130
                  final NRHS = NSVAL[IRHS];

// +    TEST 3
               // Solve and compute residual for  A * X = B.

                 srnamc.SRNAMT = 'ZLARHS';
                  zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                 srnamc.SRNAMT = 'ZSPTRS';
                  zsptrs(UPLO, N, NRHS, AFAC, IWORK, X, LDA, INFO );

               // Check error code from ZSPTRS.

                  if (INFO != 0) alaerh( PATH, 'ZSPTRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                  zspt02(UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

// +    TEST 4
               // Check solution from generated exact solution.

                  zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );

// +    TESTS 5, 6, and 7
               // Use iterative refinement to improve the solution.

                 srnamc.SRNAMT = 'ZSPRFS';
                  zsprfs(UPLO, N, NRHS, A, AFAC, IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

               // Check error code from ZSPRFS.

                  if (INFO != 0) alaerh( PATH, 'ZSPRFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) );
                  zppt05(UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 6 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 3; K <= 7; K++) { // 120
                     if ( RESULT[K] >= THRESH ) {
                        if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                        NOUT.println( 9998 )UPLO, N, NRHS, IMAT, K, RESULT[K];
                        NFAIL++;
                     }
                  } // 120
                  NRUN +=  5;
               } // 130

// +    TEST 8
               // Get an estimate of RCOND = 1/CNDNUM.

               } // 140
               ANORM = ZLANSP( '1', UPLO, N, A, RWORK );
              srnamc.SRNAMT = 'ZSPCON';
               zspcon(UPLO, N, AFAC, IWORK, ANORM, RCOND, WORK, INFO );

               // Check error code from ZSPCON.

               if (INFO != 0) alaerh( PATH, 'ZSPCON', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

               RESULT[8] = DGET06( RCOND, RCONDC );

               // Print the test ratio if it is >= THRESH.

               if ( RESULT[8] >= THRESH ) {
                  if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                  NOUT.println( 9999 )UPLO, N, IMAT, 8, RESULT( 8 );
                  NFAIL++;
               }
               NRUN++;
            } // 150
         } // 160
      } // 170

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = \'${UPLO.a1}\', N =${N.i5}, type ${IMAT.i2}, test ${.i2}, ratio =${RESULT[].g12_5};
 9998 FORMAT( ' UPLO = \'${UPLO.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${.i2}) =${RESULT[].g12_5};
      }
