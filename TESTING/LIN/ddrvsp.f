      SUBROUTINE DDRVSP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

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

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 10, NTESTS = 6 ;
      int                NFACT;
      const              NFACT = 2 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, FACT, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, K1, KL, KU, LDA, LWORK, MODE, N, NERRS, NFAIL, NIMAT, NPP, NRUN, NT;
      double             AINVNM, ANORM, CNDNUM, RCOND, RCONDC;
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, DLANSP;
      // EXTERNAL DGET06, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DCOPY, DERRVX, DGET04, DLACPY, DLARHS, DLASET, DLATB4, DLATMS, DPPT02, DPPT05, DSPSV, DSPSVX, DSPT01, DSPTRF, DSPTRI
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
      DATA               FACTS / 'F', 'N' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'SP'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10
      LWORK = MAX( 2*NMAX, NMAX*NRHS )

      // Test the error exits

      if (TSTERR) CALL DERRVX( PATH, NOUT );
      INFOT = 0

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 180
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         NPP = N*( N+1 ) / 2
         XTYPE = 'N'
         NIMAT = NTYPES
         if (N.LE.0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 170

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 170

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT.GE.3 && IMAT.LE.6
            if (ZEROT && N.LT.IMAT-2) GO TO 170;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               if ( IUPLO == 1 ) {
                  UPLO = 'U'
                  PACKIT = 'C'
               } else {
                  UPLO = 'L'
                  PACKIT = 'R'
               }

               // Set up parameters with DLATB4 and generate a test matrix
               // with DLATMS.

               dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'DLATMS'
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO );

               // Check error code from DLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 160
               }

               // For types 3-6, zero one or more rows and columns of the
               // matrix to test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT == 3 ) {
                     IZERO = 1
                  } else if ( IMAT == 4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }

                  if ( IMAT.LT.6 ) {

                     // Set row and column IZERO to zero.

                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*IZERO / 2
                        for (I = 1; I <= IZERO - 1; I++) { // 20
                           A( IOFF+I ) = ZERO
                        } // 20
                        IOFF = IOFF + IZERO
                        for (I = IZERO; I <= N; I++) { // 30
                           A( IOFF ) = ZERO
                           IOFF = IOFF + I
                        } // 30
                     } else {
                        IOFF = IZERO
                        for (I = 1; I <= IZERO - 1; I++) { // 40
                           A( IOFF ) = ZERO
                           IOFF = IOFF + N - I
                        } // 40
                        IOFF = IOFF - IZERO
                        for (I = IZERO; I <= N; I++) { // 50
                           A( IOFF+I ) = ZERO
                        } // 50
                     }
                  } else {
                     IOFF = 0
                     if ( IUPLO == 1 ) {

                        // Set the first IZERO rows and columns to zero.

                        for (J = 1; J <= N; J++) { // 70
                           I2 = MIN( J, IZERO )
                           for (I = 1; I <= I2; I++) { // 60
                              A( IOFF+I ) = ZERO
                           } // 60
                           IOFF = IOFF + J
                        } // 70
                     } else {

                        // Set the last IZERO rows and columns to zero.

                        for (J = 1; J <= N; J++) { // 90
                           I1 = MAX( J, IZERO )
                           for (I = I1; I <= N; I++) { // 80
                              A( IOFF+I ) = ZERO
                           } // 80
                           IOFF = IOFF + N - J
                        } // 90
                     }
                  }
               } else {
                  IZERO = 0
               }

               for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 150

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT )

                  // Compute the condition number for comparison with
                  // the value returned by DSPSVX.

                  if ( ZEROT ) {
                     if (IFACT == 1) GO TO 150;
                     RCONDC = ZERO

                  } else if ( IFACT == 1 ) {

                     // Compute the 1-norm of A.

                     ANORM = DLANSP( '1', UPLO, N, A, RWORK )

                     // Factor the matrix A.

                     dcopy(NPP, A, 1, AFAC, 1 );
                     dsptrf(UPLO, N, AFAC, IWORK, INFO );

                     // Compute inv(A) and take its norm.

                     dcopy(NPP, AFAC, 1, AINV, 1 );
                     dsptri(UPLO, N, AINV, IWORK, WORK, INFO );
                     AINVNM = DLANSP( '1', UPLO, N, AINV, RWORK )

                     // Compute the 1-norm condition number of A.

                     if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDC = ONE
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     }
                  }

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'DLARHS'
                  dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  XTYPE = 'C'

                  // --- Test DSPSV  ---

                  if ( IFACT == 2 ) {
                     dcopy(NPP, A, 1, AFAC, 1 );
                     dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                     // Factor the matrix and solve the system using DSPSV.

                     SRNAMT = 'DSPSV '
                     dspsv(UPLO, N, NRHS, AFAC, IWORK, X, LDA, INFO );

                     // Adjust the expected value of INFO to account for
                     // pivoting.

                     K = IZERO
                     if ( K.GT.0 ) {
                        } // 100
                        if ( IWORK( K ).LT.0 ) {
                           if ( IWORK( K ) != -K ) {
                              K = -IWORK( K )
                              GO TO 100
                           }
                        } else if ( IWORK( K ) != K ) {
                           K = IWORK( K )
                           GO TO 100
                        }
                     }

                     // Check error code from DSPSV .

                     if ( INFO != K ) {
                        alaerh(PATH, 'DSPSV ', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 120
                     } else if ( INFO != 0 ) {
                        GO TO 120
                     }

                     // Reconstruct matrix from factors and compute
                     // residual.

                     dspt01(UPLO, N, A, AFAC, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );

                     // Compute residual of the computed solution.

                     dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     dppt02(UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                     NT = 3

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 1; K <= NT; K++) { // 110
                        if ( RESULT( K ).GE.THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'DSPSV ', UPLO, N, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1
                        }
                     } // 110
                     NRUN = NRUN + NT
                     } // 120
                  }

                  // --- Test DSPSVX ---

                  if (IFACT == 2 && NPP.GT.0) CALL DLASET( 'Full', NPP, 1, ZERO, ZERO, AFAC, NPP );
                  dlaset('Full', N, NRHS, ZERO, ZERO, X, LDA );

                  // Solve the system and compute the condition number and
                  // error bounds using DSPSVX.

                  SRNAMT = 'DSPSVX'
                  dspsvx(FACT, UPLO, N, NRHS, A, AFAC, IWORK, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO );

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                  K = IZERO
                  if ( K.GT.0 ) {
                     } // 130
                     if ( IWORK( K ).LT.0 ) {
                        if ( IWORK( K ) != -K ) {
                           K = -IWORK( K )
                           GO TO 130
                        }
                     } else if ( IWORK( K ) != K ) {
                        K = IWORK( K )
                        GO TO 130
                     }
                  }

                  // Check the error code from DSPSVX.

                  if ( INFO != K ) {
                     alaerh(PATH, 'DSPSVX', INFO, K, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 150
                  }

                  if ( INFO == 0 ) {
                     if ( IFACT.GE.2 ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        dspt01(UPLO, N, A, AFAC, IWORK, AINV, LDA, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                     // Compute residual of the computed solution.

                     dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     dppt02(UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

                     // Check the error bounds from iterative refinement.

                     dppt05(UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                  } else {
                     K1 = 6
                  }

                  // Compare RCOND from DSPSVX with the computed value
                  // in RCONDC.

                  RESULT( 6 ) = DGET06( RCOND, RCONDC )

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = K1; K <= 6; K++) { // 140
                     if ( RESULT( K ).GE.THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )'DSPSVX', FACT, UPLO, N, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1
                     }
                  } // 140
                  NRUN = NRUN + 7 - K1

               } // 150

            } // 160
         } // 170
      } // 180

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN

      // End of DDRVSP

      }
