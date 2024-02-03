      void ddrvsy(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      double             A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 10, NTESTS = 6 ;
      int                NFACT;
      const              NFACT = 2 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, FACT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, K1, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT;
      double             AINVNM, ANORM, CNDNUM, RCOND, RCONDC;
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, DLANSY;
      // EXTERNAL DGET06, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DERRVX, DGET04, DLACPY, DLARHS, DLASET, DLATB4, DLATMS, DPOT02, DPOT05, DSYSV, DSYSVX, DSYT01, DSYTRF, DSYTRI2, XLAENV
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
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', FACTS = 'F', 'N';
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'SY';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10
      LWORK = max( 2*NMAX, NMAX*NRHS );

      // Test the error exits

      if (TSTERR) CALL DERRVX( PATH, NOUT );
      INFOT = 0;

      // Set the block size and minimum block size for testing.

      NB = 1;
      NBMIN = 2;
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 180
         N = NVAL( IN );
         LDA = max( N, 1 );
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 170

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 170;

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 6;
            if (ZEROT && N < IMAT-2) GO TO 170;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               UPLO = UPLOS( IUPLO );

               // Set up parameters with DLATB4 and generate a test matrix
               // with DLATMS.

               dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'DLATMS';
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from DLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 160;
               }

               // For types 3-6, zero one or more rows and columns of the
               // matrix to test that INFO is returned correctly.

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
                           A( IOFF+I ) = ZERO;
                        } // 20
                        IOFF = IOFF + IZERO;
                        for (I = IZERO; I <= N; I++) { // 30
                           A( IOFF ) = ZERO;
                           IOFF = IOFF + LDA;
                        } // 30
                     } else {
                        IOFF = IZERO;
                        for (I = 1; I <= IZERO - 1; I++) { // 40
                           A( IOFF ) = ZERO;
                           IOFF = IOFF + LDA;
                        } // 40
                        IOFF = IOFF - IZERO;
                        for (I = IZERO; I <= N; I++) { // 50
                           A( IOFF+I ) = ZERO;
                        } // 50
                     }
                  } else {
                     IOFF = 0;
                     if ( IUPLO == 1 ) {

                        // Set the first IZERO rows and columns to zero.

                        for (J = 1; J <= N; J++) { // 70
                           I2 = min( J, IZERO );
                           for (I = 1; I <= I2; I++) { // 60
                              A( IOFF+I ) = ZERO;
                           } // 60
                           IOFF = IOFF + LDA;
                        } // 70
                     } else {

                        // Set the last IZERO rows and columns to zero.

                        for (J = 1; J <= N; J++) { // 90
                           I1 = max( J, IZERO );
                           for (I = I1; I <= N; I++) { // 80
                              A( IOFF+I ) = ZERO;
                           } // 80
                           IOFF = IOFF + LDA;
                        } // 90
                     }
                  }
               } else {
                  IZERO = 0;
               }

               for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 150

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT );

                  // Compute the condition number for comparison with
                  // the value returned by DSYSVX.

                  if ( ZEROT ) {
                     if (IFACT == 1) GO TO 150;
                     RCONDC = ZERO;

                  } else if ( IFACT == 1 ) {

                     // Compute the 1-norm of A.

                     ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK );

                     // Factor the matrix A.

                     dlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     dsytrf(UPLO, N, AFAC, LDA, IWORK, WORK, LWORK, INFO );

                     // Compute inv(A) and take its norm.

                     dlacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                     LWORK = (N+NB+1)*(NB+3);
                     dsytri2(UPLO, N, AINV, LDA, IWORK, WORK, LWORK, INFO );
                     AINVNM = DLANSY( '1', UPLO, N, AINV, LDA, RWORK );

                     // Compute the 1-norm condition number of A.

                     if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                        RCONDC = ONE;
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM;
                     }
                  }

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'DLARHS';
                  dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  XTYPE = 'C';

                  // --- Test DSYSV  ---

                  if ( IFACT == 2 ) {
                     dlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                     // Factor the matrix and solve the system using DSYSV.

                     SRNAMT = 'DSYSV ';
                     dsysv(UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, WORK, LWORK, INFO );

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

                     // Check error code from DSYSV .

                     if ( INFO != K ) {
                        alaerh(PATH, 'DSYSV ', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 120;
                     } else if ( INFO != 0 ) {
                        GO TO 120;
                     }

                     // Reconstruct matrix from factors and compute
                     // residual.

                     dsyt01(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );

                     // Compute residual of the computed solution.

                     dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     dpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                     NT = 3;

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 1; K <= NT; K++) { // 110
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                           WRITE( NOUT, FMT = 9999 )'DSYSV ', UPLO, N, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1;
                        }
                     } // 110
                     NRUN = NRUN + NT;
                     } // 120
                  }

                  // --- Test DSYSVX ---

                  if (IFACT == 2) CALL DLASET( UPLO, N, N, ZERO, ZERO, AFAC, LDA );
                  dlaset('Full', N, NRHS, ZERO, ZERO, X, LDA );

                  // Solve the system and compute the condition number and
                  // error bounds using DSYSVX.

                  SRNAMT = 'DSYSVX';
                  dsysvx(FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, LWORK, IWORK( N+1 ), INFO );

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

                  // Check the error code from DSYSVX.

                  if ( INFO != K ) {
                     alaerh(PATH, 'DSYSVX', INFO, K, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 150;
                  }

                  if ( INFO == 0 ) {
                     if ( IFACT >= 2 ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        dsyt01(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1;
                     } else {
                        K1 = 2;
                     }

                     // Compute residual of the computed solution.

                     dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     dpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

                     // Check the error bounds from iterative refinement.

                     dpot05(UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                  } else {
                     K1 = 6;
                  }

                  // Compare RCOND from DSYSVX with the computed value
                  // in RCONDC.

                  RESULT( 6 ) = DGET06( RCOND, RCONDC );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = K1; K <= 6; K++) { // 140
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                        WRITE( NOUT, FMT = 9998 )'DSYSVX', FACT, UPLO, N, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 140
                  NRUN = NRUN + 7 - K1;

               } // 150

            } // 160
         } // 170
      } // 180

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 );
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 );
      return;

      // End of DDRVSY

      }
