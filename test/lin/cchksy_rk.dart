      void cchksy_rk(DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, E, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double               RWORK( * );
      Complex            A( * ), AFAC( * ), AINV( * ), B( * ), E( * ), WORK( * ), X( * ), XACT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ONEHALF;
      const              ONEHALF = 0.5 ;
      double               EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      Complex            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, TYPE, UPLO, XTYPE;
      String             PATH, MATPATH;
      int                I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, ITEMP, ITEMP2, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      double               ALPHA, ANORM, CNDNUM, CONST, SING_MAX, SING_MIN, RCOND, RCONDC, STEMP;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS );
      Complex            BLOCK( 2, 2 ), CDUMMY( 1 );
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, CLANSY, SGET06;
      // EXTERNAL CLANGE, CLANSY, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CERRSY, CGESVD, CGET04, CLACPY, CLARHS, CLATB4, CLATMS, CLATSY, CSYT02, CSYT03, CSYCON_3, CSYT01_3, CSYTRF_RK, CSYTRI_3, CSYTRS_3, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
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
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      // Test path

      PATH[1: 1] = 'Complex precision';
      PATH[2: 3] = 'SK';

      // Path to generate matrices

      MATPATH[1: 1] = 'Complex precision';
      MATPATH[2: 3] = 'SY';

      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) cerrsy( PATH, NOUT );
      INFOT = 0;

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

               // Begin generate test matrix A.

               if ( IMAT != NTYPES ) {

                  // Set up parameters with CLATB4 for the matrix generator
                  // based on the type of matrix to be generated.

                  clatb4(MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                  // Generate a matrix with CLATMS.

                  SRNAMT = 'CLATMS';
                  clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

                  // Check error code from CLATMS and handle error.

                  if ( INFO != 0 ) {
                     alaerh(PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

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

               } else {

                  // For matrix kind IMAT = 11, generate special block
                  // diagonal matrix to test alternate code
                  // for the 2 x 2 blocks.

                  clatsy(UPLO, N, A, LDA, ISEED );

               }

               // End generate test matrix A.


               // Do for each value of NB in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 240

                  // Set the optimal blocksize, which will be later
                  // returned by ILAENV.

                  NB = NBVAL( INB );
                  xlaenv(1, NB );

                  // Copy the test matrix A into matrix AFAC which
                  // will be factorized in place. This is needed to
                  // preserve the test matrix A for subsequent tests.

                  clacpy(UPLO, N, N, A, LDA, AFAC, LDA );

                  // Compute the L*D*L**T or U*D*U**T factorization of the
                  // matrix. IWORK stores details of the interchanges and
                  // the block structure of D. AINV is a work array for
                  // block factorization, LWORK is the length of AINV.

                  LWORK = max( 2, NB )*LDA;
                  SRNAMT = 'CSYTRF_RK';
                  csytrf_rk(UPLO, N, AFAC, LDA, E, IWORK, AINV, LWORK, INFO );

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

                  // Check error code from CSYTRF_RK and handle error.

                  if (INFO != K) alaerh( PATH, 'CSYTRF_RK', INFO, K, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                  // Set the condition estimate flag if the INFO is not 0.

                  if ( INFO != 0 ) {
                     TRFCON = true;
                  } else {
                     TRFCON = false;
                  }

// +    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  csyt01_3(UPLO, N, A, LDA, AFAC, LDA, E, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );
                  NT = 1;

// +    TEST 2
                  // Form the inverse and compute the residual,
                  // if the factorization was competed without INFO > 0
                  // (i.e. there is no zero rows and columns).
                  // Do it only for the first block size.

                  if ( INB == 1 && !TRFCON ) {
                     clacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                     SRNAMT = 'CSYTRI_3';

                     // Another reason that we need to compute the inverse
                     // is that CSYT03 produces RCONDC which is used later
                     // in TEST6 and TEST7.

                     LWORK = (N+NB+1)*(NB+3);
                     csytri_3(UPLO, N, AINV, LDA, E, IWORK, WORK, LWORK, INFO );

                     // Check error code from CSYTRI_3 and handle error.

                     if (INFO != 0) alaerh( PATH, 'CSYTRI_3', INFO, -1, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     // Compute the residual for a symmetric matrix times
                     // its inverse.

                     csyt03(UPLO, N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );
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
                  STEMP = ZERO;

                  CONST = ( ( ALPHA**2-ONE ) / ( ALPHA**2-ONEHALF ) ) / ( ONE-ALPHA );

                  if ( IUPLO == 1 ) {

                  // Compute largest element in U

                     K = N;
                     } // 120
                     if (K <= 1) GO TO 130;

                     if ( IWORK( K ) > ZERO ) {

                        // Get max absolute value from elements
                        // in column k in in U

                        STEMP = CLANGE( 'M', K-1, 1, AFAC( ( K-1 )*LDA+1 ), LDA, RWORK );
                     } else {

                        // Get max absolute value from elements
                        // in columns k and k-1 in U

                        STEMP = CLANGE( 'M', K-2, 2, AFAC( ( K-2 )*LDA+1 ), LDA, RWORK );
                        K = K - 1;

                     }

                     // STEMP should be bounded by CONST

                     STEMP = STEMP - CONST + THRESH;
                     if[STEMP > RESULT( 3 ) ) RESULT( 3] = STEMP;

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

                        STEMP = CLANGE( 'M', N-K, 1, AFAC( ( K-1 )*LDA+K+1 ), LDA, RWORK );
                     } else {

                        // Get max absolute value from elements
                        // in columns k and k+1 in L

                        STEMP = CLANGE( 'M', N-K-1, 2, AFAC( ( K-1 )*LDA+K+2 ), LDA, RWORK );
                        K = K + 1;

                     }

                     // STEMP should be bounded by CONST

                     STEMP = STEMP - CONST + THRESH;
                     if[STEMP > RESULT( 3 ) ) RESULT( 3] = STEMP;

                     K = K + 1;

                     GO TO 140;
                     } // 150
                  }


// +    TEST 4
                  // Compute largest 2-Norm (condition number)
                  // of 2-by-2 diag blocks

                  RESULT[4] = ZERO;
                  STEMP = ZERO;

                  CONST = ( ( ALPHA**2-ONE ) / ( ALPHA**2-ONEHALF ) )* ( ( ONE + ALPHA ) / ( ONE - ALPHA ) );

                  if ( IUPLO == 1 ) {

                     // Loop backward for UPLO = 'U'

                     K = N;
                     } // 160
                     if (K <= 1) GO TO 170;

                     if ( IWORK( K ) < ZERO ) {

                        // Get the two singular values
                        // (real and non-negative) of a 2-by-2 block,
                        // store them in RWORK array

                        BLOCK[1, 1] = AFAC( ( K-2 )*LDA+K-1 );
                        BLOCK[1, 2] = E( K );
                        BLOCK[2, 1] = BLOCK( 1, 2 );
                        BLOCK[2, 2] = AFAC( (K-1)*LDA+K );

                        cgesvd('N', 'N', 2, 2, BLOCK, 2, RWORK, CDUMMY, 1, CDUMMY, 1, WORK, 6, RWORK( 3 ), INFO );


                        SING_MAX = RWORK( 1 );
                        SING_MIN = RWORK( 2 );

                        STEMP = SING_MAX / SING_MIN;

                        // STEMP should be bounded by CONST

                        STEMP = STEMP - CONST + THRESH;
                        if[STEMP > RESULT( 4 ) ) RESULT( 4] = STEMP;
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

                        BLOCK[1, 1] = AFAC( ( K-1 )*LDA+K );
                        BLOCK[2, 1] = E( K );
                        BLOCK[1, 2] = BLOCK( 2, 1 );
                        BLOCK[2, 2] = AFAC( K*LDA+K+1 );

                        cgesvd('N', 'N', 2, 2, BLOCK, 2, RWORK, CDUMMY, 1, CDUMMY, 1, WORK, 6, RWORK(3), INFO );

                        SING_MAX = RWORK( 1 );
                        SING_MIN = RWORK( 2 );

                        STEMP = SING_MAX / SING_MIN;

                        // STEMP should be bounded by CONST

                        STEMP = STEMP - CONST + THRESH;
                        if[STEMP > RESULT( 4 ) ) RESULT( 4] = STEMP;
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

// +    TEST 5 ( Using TRS_3)
                  // Solve and compute residual for  A * X = B.

                     // Choose a set of NRHS random solution vectors
                     // stored in XACT and set up the right hand side B

                     SRNAMT = 'CLARHS';
                     clarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     clacpy('Full', N, NRHS, B, LDA, X, LDA );

                     SRNAMT = 'CSYTRS_3';
                     csytrs_3(UPLO, N, NRHS, AFAC, LDA, E, IWORK, X, LDA, INFO );

                     // Check error code from CSYTRS_3 and handle error.

                     if (INFO != 0) alaerh( PATH, 'CSYTRS_3', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     clacpy('Full', N, NRHS, B, LDA, WORK, LDA );

                     // Compute the residual for the solution

                     csyt02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 5 ) );

// +    TEST 6
                  // Check solution from generated exact solution.

                     cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 6 ) );

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
                  ANORM = CLANSY( '1', UPLO, N, A, LDA, RWORK );
                  SRNAMT = 'CSYCON_3';
                  csycon_3(UPLO, N, AFAC, LDA, E, IWORK, ANORM, RCOND, WORK, INFO );

                  // Check error code from CSYCON_3 and handle error.

                  if (INFO != 0) alaerh( PATH, 'CSYCON_3', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  // Compute the test ratio to compare values of RCOND

                  RESULT[7] = SGET06( RCOND, RCONDC );

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

 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ', I2, ', test ', I2, ', ratio =', G12.5 );
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') =', G12.5 );
 9997 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ',', 10X, ' type ', I2, ', test(', I2, ') =', G12.5 );
      return;
      }