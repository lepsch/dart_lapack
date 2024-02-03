      SUBROUTINE SDRVSY( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

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
      REAL               A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 10, NTESTS = 6 ;
      int                NFACT;
      const              NFACT = 2 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, EQUED, FACT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, K1, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT, N_ERR_BNDS;
      REAL               AINVNM, ANORM, CNDNUM, RCOND, RCONDC, RPVGRW_SVXX
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS ), BERR( NRHS ), ERRBNDS_N( NRHS, 3 ), ERRBNDS_C( NRHS, 3 )
      // ..
      // .. External Functions ..
      REAL               SGET06, SLANSY
      // EXTERNAL SGET06, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, SERRVX, SGET04, SLACPY, SLARHS, SLASET, SLATB4, SLATMS, SPOT02, SPOT05, SSYSV, SSYSVX, SSYT01, SSYTRF, SSYTRI2, XLAENV, SSYSVXX
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

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'SY'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10
      LWORK = MAX( 2*NMAX, NMAX*NRHS )

      // Test the error exits

      if (TSTERR) CALL SERRVX( PATH, NOUT );
      INFOT = 0

      // Set the block size and minimum block size for testing.

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

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            if (ZEROT .AND. N.LT.IMAT-2) GO TO 170;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               UPLO = UPLOS( IUPLO )

               // Set up parameters with SLATB4 and generate a test matrix
               // with SLATMS.

               slatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'SLATMS'
               slatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from SLATMS.

               if ( INFO.NE.0 ) {
                  alaerh(PATH, 'SLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 160
               }

               // For types 3-6, zero one or more rows and columns of the
               // matrix to test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT.EQ.3 ) {
                     IZERO = 1
                  } else if ( IMAT.EQ.4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }

                  if ( IMAT.LT.6 ) {

                     // Set row and column IZERO to zero.

                     if ( IUPLO.EQ.1 ) {
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
                     IOFF = 0
                     if ( IUPLO.EQ.1 ) {

                        // Set the first IZERO rows and columns to zero.

                        for (J = 1; J <= N; J++) { // 70
                           I2 = MIN( J, IZERO )
                           for (I = 1; I <= I2; I++) { // 60
                              A( IOFF+I ) = ZERO
                           } // 60
                           IOFF = IOFF + LDA
                        } // 70
                     } else {

                        // Set the last IZERO rows and columns to zero.

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

               for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 150

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT )

                  // Compute the condition number for comparison with
                  // the value returned by SSYSVX.

                  if ( ZEROT ) {
                     if (IFACT.EQ.1) GO TO 150;
                     RCONDC = ZERO

                  } else if ( IFACT.EQ.1 ) {

                     // Compute the 1-norm of A.

                     ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )

                     // Factor the matrix A.

                     slacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     ssytrf(UPLO, N, AFAC, LDA, IWORK, WORK, LWORK, INFO );

                     // Compute inv(A) and take its norm.

                     slacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                     LWORK = (N+NB+1)*(NB+3)
                     ssytri2(UPLO, N, AINV, LDA, IWORK, WORK, LWORK, INFO );
                     AINVNM = SLANSY( '1', UPLO, N, AINV, LDA, RWORK )

                     // Compute the 1-norm condition number of A.

                     if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDC = ONE
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     }
                  }

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'SLARHS'
                  slarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  XTYPE = 'C'

                  // --- Test SSYSV  ---

                  if ( IFACT.EQ.2 ) {
                     slacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     slacpy('Full', N, NRHS, B, LDA, X, LDA );

                     // Factor the matrix and solve the system using SSYSV.

                     SRNAMT = 'SSYSV '
                     ssysv(UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, WORK, LWORK, INFO );

                     // Adjust the expected value of INFO to account for
                     // pivoting.

                     K = IZERO
                     if ( K.GT.0 ) {
                        } // 100
                        if ( IWORK( K ).LT.0 ) {
                           if ( IWORK( K ).NE.-K ) {
                              K = -IWORK( K )
                              GO TO 100
                           }
                        } else if ( IWORK( K ).NE.K ) {
                           K = IWORK( K )
                           GO TO 100
                        }
                     }

                     // Check error code from SSYSV .

                     if ( INFO.NE.K ) {
                        alaerh(PATH, 'SSYSV ', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 120
                     } else if ( INFO.NE.0 ) {
                        GO TO 120
                     }

                     // Reconstruct matrix from factors and compute
                     // residual.

                     ssyt01(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );

                     // Compute residual of the computed solution.

                     slacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     spot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                     NT = 3

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 1; K <= NT; K++) { // 110
                        if ( RESULT( K ).GE.THRESH ) {
                           if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'SSYSV ', UPLO, N, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1
                        }
                     } // 110
                     NRUN = NRUN + NT
                     } // 120
                  }

                  // --- Test SSYSVX ---

                  if (IFACT.EQ.2) CALL SLASET( UPLO, N, N, ZERO, ZERO, AFAC, LDA );
                  slaset('Full', N, NRHS, ZERO, ZERO, X, LDA );

                  // Solve the system and compute the condition number and
                  // error bounds using SSYSVX.

                  SRNAMT = 'SSYSVX'
                  ssysvx(FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, LWORK, IWORK( N+1 ), INFO );

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                  K = IZERO
                  if ( K.GT.0 ) {
                     } // 130
                     if ( IWORK( K ).LT.0 ) {
                        if ( IWORK( K ).NE.-K ) {
                           K = -IWORK( K )
                           GO TO 130
                        }
                     } else if ( IWORK( K ).NE.K ) {
                        K = IWORK( K )
                        GO TO 130
                     }
                  }

                  // Check the error code from SSYSVX.

                  if ( INFO.NE.K ) {
                     alaerh(PATH, 'SSYSVX', INFO, K, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 150
                  }

                  if ( INFO.EQ.0 ) {
                     if ( IFACT.GE.2 ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        ssyt01(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                     // Compute residual of the computed solution.

                     slacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     spot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

                     // Check the error bounds from iterative refinement.

                     spot05(UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                  } else {
                     K1 = 6
                  }

                  // Compare RCOND from SSYSVX with the computed value
                  // in RCONDC.

                  RESULT( 6 ) = SGET06( RCOND, RCONDC )

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = K1; K <= 6; K++) { // 140
                     if ( RESULT( K ).GE.THRESH ) {
                        if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )'SSYSVX', FACT, UPLO, N, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1
                     }
                  } // 140
                  NRUN = NRUN + 7 - K1

                  // --- Test SSYSVXX ---

                  // Restore the matrices A and B.

                  if (IFACT.EQ.2) CALL SLASET( UPLO, N, N, ZERO, ZERO, AFAC, LDA );
                  slaset('Full', N, NRHS, ZERO, ZERO, X, LDA );

                  // Solve the system and compute the condition number
                  // and error bounds using SSYSVXX.

                  SRNAMT = 'SSYSVXX'
                  N_ERR_BNDS = 3
                  EQUED = 'N'
                  ssysvxx(FACT, UPLO, N, NRHS, A, LDA, AFAC, LDA, IWORK, EQUED, WORK( N+1 ), B, LDA, X, LDA, RCOND, RPVGRW_SVXX, BERR, N_ERR_BNDS, ERRBNDS_N, ERRBNDS_C, 0, ZERO, WORK, IWORK( N+1 ), INFO );

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                  K = IZERO
                  if ( K.GT.0 ) {
                     } // 135
                     if ( IWORK( K ).LT.0 ) {
                        if ( IWORK( K ).NE.-K ) {
                           K = -IWORK( K )
                           GO TO 135
                        }
                     } else if ( IWORK( K ).NE.K ) {
                        K = IWORK( K )
                        GO TO 135
                     }
                  }

                  // Check the error code from SSYSVXX.

                  if ( INFO.NE.K .AND. INFO.LE.N ) {
                     alaerh(PATH, 'SSYSVXX', INFO, K, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                     GO TO 150
                  }

                  if ( INFO.EQ.0 ) {
                     if ( IFACT.GE.2 ) {

                  // Reconstruct matrix from factors and compute
                  // residual.

                        ssyt01(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK(2*NRHS+1), RESULT( 1 ) );
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                  // Compute residual of the computed solution.

                     slacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     spot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                  // Check solution from generated exact solution.

                     sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

                  // Check the error bounds from iterative refinement.

                     spot05(UPLO, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                  } else {
                     K1 = 6
                  }

                  // Compare RCOND from SSYSVXX with the computed value
                  // in RCONDC.

                  RESULT( 6 ) = SGET06( RCOND, RCONDC )

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = K1; K <= 6; K++) { // 85
                     if ( RESULT( K ).GE.THRESH ) {
                        if (NFAIL.EQ.0 .AND. NERRS.EQ.0) CALL ALADHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )'SSYSVXX', FACT, UPLO, N, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1
                     }
                  } // 85
                  NRUN = NRUN + 7 - K1

               } // 150

            } // 160
         } // 170
      } // 180

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );


      // Test Error Bounds from SSYSVXX

      sebchvxx(THRESH, PATH);

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN

      // End of SDRVSYX

      }
