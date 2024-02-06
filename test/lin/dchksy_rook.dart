import 'common.dart';

      void dchksy_rook(DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double             A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      int                NTYPES;
      const              NTYPES = 10 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      bool               TRFCON, ZEROT;
      String             DIST, TYPE, UPLO, XTYPE;
      String             PATH, MATPATH;
      int                I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      double             ALPHA, ANORM, CNDNUM, CONST, DTEMP, SING_MAX, SING_MIN, RCOND, RCONDC;
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             BLOCK( 2, 2 ), DDUMMY( 1 ), RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- double             DGET06, DLANGE, DLANSY;
      // EXTERNAL DGET06, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DERRSY, DGET04, DLACPY, DLARHS, DLATB4, DLATMS, DPOT02, DPOT03, DGESVD, DSYCON_ROOK, DSYT01_ROOK, DSYTRF_ROOK, DSYTRI_ROOK, DSYTRS_ROOK, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];

      // Initialize constants and the random number seed.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      // Test path

      PATH[1: 1] = 'double          ';
      PATH[2: 3] = 'SR';

      // Path to generate matrices

      MATPATH[1: 1] = 'double          ';
      MATPATH[2: 3] = 'SY';

      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) derrsy( PATH, NOUT );
      infoc.INFOT = 0;

      // Set the minimum block size for which the block routine should
      // be used, which will be later returned by ILAENV

      xlaenv(2, 2 );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 270
         N = NVAL( IN );
         LDA = max( N, 1 );
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         IZERO = 0;

         // Do for each value of matrix type IMAT

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 260

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 260;

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 6;
            if (ZEROT && N < IMAT-2) GO TO 260;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 250
               UPLO = UPLOS( IUPLO );

               // Begin generate the test matrix A.

               // Set up parameters with DLATB4 for the matrix generator
               // based on the type of matrix to be generated.

               dlatb4(MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               // Generate a matrix with DLATMS.

               srnamc.SRNAMT = 'DLATMS';
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from DLATMS and handle error.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  // Skip all tests for this generated matrix

                  GO TO 250;
               }

               // For matrix types 3-6, zero one or more rows and
               // columns of the matrix to test that INFO is returned
               // correctly.

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

                        // Set the first IZERO rows and columns to zero.

                        IOFF = 0;
                        for (J = 1; J <= N; J++) { // 70
                           I2 = min( J, IZERO );
                           for (I = 1; I <= I2; I++) { // 60
                              A[IOFF+I] = ZERO;
                           } // 60
                           IOFF = IOFF + LDA;
                        } // 70
                     } else {

                        // Set the last IZERO rows and columns to zero.

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

               // End generate the test matrix A.


               // Do for each value of NB in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 240

                  // Set the optimal blocksize, which will be later
                  // returned by ILAENV.

                  NB = NBVAL( INB );
                  xlaenv(1, NB );

                  // Copy the test matrix A into matrix AFAC which
                  // will be factorized in place. This is needed to
                  // preserve the test matrix A for subsequent tests.

                  dlacpy(UPLO, N, N, A, LDA, AFAC, LDA );

                  // Compute the L*D*L**T or U*D*U**T factorization of the
                  // matrix. IWORK stores details of the interchanges and
                  // the block structure of D. AINV is a work array for
                  // block factorization, LWORK is the length of AINV.

                  LWORK = max( 2, NB )*LDA;
                  srnamc.SRNAMT = 'DSYTRF_ROOK';
                  dsytrf_rook(UPLO, N, AFAC, LDA, IWORK, AINV, LWORK, INFO );

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

                  // Check error code from DSYTRF_ROOK and handle error.

                  if (INFO != K) alaerh( PATH, 'DSYTRF_ROOK', INFO, K, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                  // Set the condition estimate flag if the INFO is not 0.

                  if ( INFO != 0 ) {
                     TRFCON = true;
                  } else {
                     TRFCON = false;
                  }

// +    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  dsyt01_rook(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );
                  NT = 1;

// +    TEST 2
                  // Form the inverse and compute the residual,
                  // if the factorization was competed without INFO > 0
                  // (i.e. there is no zero rows and columns).
                  // Do it only for the first block size.

                  if ( INB == 1 && !TRFCON ) {
                     dlacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                     srnamc.SRNAMT = 'DSYTRI_ROOK';
                     dsytri_rook(UPLO, N, AINV, LDA, IWORK, WORK, INFO );

                     // Check error code from DSYTRI_ROOK and handle error.

                     if (INFO != 0) alaerh( PATH, 'DSYTRI_ROOK', INFO, -1, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     // Compute the residual for a symmetric matrix times
                     // its inverse.

                     dpot03(UPLO, N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );
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

// +    TEST 3
                  // Compute largest element in U or L

                  RESULT[3] = ZERO;
                  DTEMP = ZERO;

                  CONST = ONE / ( ONE-ALPHA );

                  if ( IUPLO == 1 ) {

                  // Compute largest element in U

                     K = N;
                     } // 120
                     if (K <= 1) GO TO 130;

                     if ( IWORK( K ) > ZERO ) {

                        // Get max absolute value from elements
                        // in column k in in U

                        DTEMP = dlange( 'M', K-1, 1, AFAC( ( K-1 )*LDA+1 ), LDA, RWORK );
                     } else {

                        // Get max absolute value from elements
                        // in columns k and k-1 in U

                        DTEMP = dlange( 'M', K-2, 2, AFAC( ( K-2 )*LDA+1 ), LDA, RWORK );
                        K = K - 1;

                     }

                     // DTEMP should be bounded by CONST

                     DTEMP = DTEMP - CONST + THRESH;
                     if[DTEMP > RESULT( 3 ) ) RESULT( 3] = DTEMP;

                     K = K - 1;

                     GO TO 120;
                     } // 130

                  } else {

                  // Compute largest element in L

                     K = 1;
                     } // 140
                     if (K >= N) GO TO 150;

                     if ( IWORK( K ) > ZERO ) {

                        // Get max absolute value from elements
                        // in column k in in L

                        DTEMP = dlange( 'M', N-K, 1, AFAC( ( K-1 )*LDA+K+1 ), LDA, RWORK );
                     } else {

                        // Get max absolute value from elements
                        // in columns k and k+1 in L

                        DTEMP = dlange( 'M', N-K-1, 2, AFAC( ( K-1 )*LDA+K+2 ), LDA, RWORK );
                        K = K + 1;

                     }

                     // DTEMP should be bounded by CONST

                     DTEMP = DTEMP - CONST + THRESH;
                     if[DTEMP > RESULT( 3 ) ) RESULT( 3] = DTEMP;

                     K = K + 1;

                     GO TO 140;
                     } // 150
                  }


// +    TEST 4
                  // Compute largest 2-Norm (condition number)
                  // of 2-by-2 diag blocks

                  RESULT[4] = ZERO;
                  DTEMP = ZERO;

                  CONST = ( ONE+ALPHA ) / ( ONE-ALPHA );
                  dlacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );

                  if ( IUPLO == 1 ) {

                     // Loop backward for UPLO = 'U'

                     K = N;
                     } // 160
                     if (K <= 1) GO TO 170;

                     if ( IWORK( K ) < ZERO ) {

                        // Get the two singular values
                        // (real and non-negative) of a 2-by-2 block,
                        // store them in RWORK array

                        BLOCK[1][1] = AFAC( ( K-2 )*LDA+K-1 );
                        BLOCK[1][2] = AFAC( (K-1)*LDA+K-1 );
                        BLOCK[2][1] = BLOCK( 1, 2 );
                        BLOCK[2][2] = AFAC( (K-1)*LDA+K );

                        dgesvd('N', 'N', 2, 2, BLOCK, 2, RWORK, DDUMMY, 1, DDUMMY, 1, WORK, 10, INFO );

                        SING_MAX = RWORK( 1 );
                        SING_MIN = RWORK( 2 );

                        DTEMP = SING_MAX / SING_MIN;

                        // DTEMP should be bounded by CONST

                        DTEMP = DTEMP - CONST + THRESH;
                        if[DTEMP > RESULT( 4 ) ) RESULT( 4] = DTEMP;
                        K = K - 1;

                     }

                     K = K - 1;

                     GO TO 160;
                     } // 170

                  } else {

                     // Loop forward for UPLO = 'L'

                     K = 1;
                     } // 180
                     if (K >= N) GO TO 190;

                     if ( IWORK( K ) < ZERO ) {

                        // Get the two singular values
                        // (real and non-negative) of a 2-by-2 block,
                        // store them in RWORK array

                        BLOCK[1][1] = AFAC( ( K-1 )*LDA+K );
                        BLOCK[2][1] = AFAC( ( K-1 )*LDA+K+1 );
                        BLOCK[1][2] = BLOCK( 2, 1 );
                        BLOCK[2][2] = AFAC( K*LDA+K+1 );

                        dgesvd('N', 'N', 2, 2, BLOCK, 2, RWORK, DDUMMY, 1, DDUMMY, 1, WORK, 10, INFO );


                        SING_MAX = RWORK( 1 );
                        SING_MIN = RWORK( 2 );

                        DTEMP = SING_MAX / SING_MIN;

                        // DTEMP should be bounded by CONST

                        DTEMP = DTEMP - CONST + THRESH;
                        if[DTEMP > RESULT( 4 ) ) RESULT( 4] = DTEMP;
                        K = K + 1;

                     }

                     K = K + 1;

                     GO TO 180;
                     } // 190
                  }

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 3; K <= 4; K++) { // 200
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 200
                  NRUN = NRUN + 2;

                  // Skip the other tests if this is not the first block
                  // size.

                  if (INB > 1) GO TO 240;

                  // Do only the condition estimate if INFO is not 0.

                  if ( TRFCON ) {
                     RCONDC = ZERO;
                     GO TO 230;
                  }

                  // Do for each value of NRHS in NSVAL.

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 220
                     NRHS = NSVAL( IRHS );

// +    TEST 5 ( Using TRS_ROOK)
                  // Solve and compute residual for  A * X = B.

                     // Choose a set of NRHS random solution vectors
                     // stored in XACT and set up the right hand side B

                     srnamc.SRNAMT = 'DLARHS';
                     dlarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                     srnamc.SRNAMT = 'DSYTRS_ROOK';
                     dsytrs_rook(UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO );

                     // Check error code from DSYTRS_ROOK and handle error.

                     if (INFO != 0) alaerh( PATH, 'DSYTRS_ROOK', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );

                     // Compute the residual for the solution

                     dpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 5 ) );

// +    TEST 6
                  // Check solution from generated exact solution.

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 6 ) );

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 5; K <= 6; K++) { // 210
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                           WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1;
                        }
                     } // 210
                     NRUN = NRUN + 2;

                  // End do for each value of NRHS in NSVAL.

                  } // 220

// +    TEST 7
                  // Get an estimate of RCOND = 1/CNDNUM.

                  } // 230
                  ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK );
                  srnamc.SRNAMT = 'DSYCON_ROOK';
                  dsycon_rook(UPLO, N, AFAC, LDA, IWORK, ANORM, RCOND, WORK, IWORK( N+1 ), INFO );

                  // Check error code from DSYCON_ROOK and handle error.

                  if (INFO != 0) alaerh( PATH, 'DSYCON_ROOK', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  // Compute the test ratio to compare to values of RCOND

                  RESULT[7] = DGET06( RCOND, RCONDC );

                  // Print information about the tests that did not pass
                  // the threshold.

                  if ( RESULT( 7 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 7, RESULT( 7 );
                     NFAIL = NFAIL + 1;
                  }
                  NRUN = NRUN + 1;
               } // 240

            } // 250
         } // 260
      } // 270

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = ''${.a1}'', N =${.i5}, NB =${.i4}, type ${.i2}, test ${.i2}, ratio =${.g12_5};
 9998 FORMAT( ' UPLO = ''${.a1}'', N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) =${.g12_5};
 9997 FORMAT( ' UPLO = ''${.a1}'', N =${.i5},${' ' * 10} type ${.i2}, test(${.i2}) =${.g12_5};
      return;
      }
