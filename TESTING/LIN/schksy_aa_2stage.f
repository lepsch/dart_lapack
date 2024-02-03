      SUBROUTINE SCHKSY_AA_2STAGE( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // IMPLICIT NONE

      // .. Scalar Arguments ..
      bool         TSTERR;
      int          NN, NNB, NNS, NMAX, NOUT;
      REAL         THRESH;
      // ..
      // .. Array Arguments ..
      bool         DOTYPE( * );
      int          IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      REAL         A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL         ZERO;
      const        ZERO = 0.0 ;
      int          NTYPES;
      const        NTYPES = 10 ;
      int          NTESTS;
      const        NTESTS = 9 ;
      // ..
      // .. Local Scalars ..
      bool         ZEROT;
      String       DIST, TYPE, UPLO, XTYPE;
      String       PATH, MATPATH;
      int          I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      REAL         ANORM, CNDNUM;
      // ..
      // .. Local Arrays ..
      String       UPLOS( 2 );
      int          ISEED( 4 ), ISEEDY( 4 );
      REAL         RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SERRSY, SLACPY, SLARHS, SLATB4, SLATMS, SPOT02, SSYT01_AA, SSYTRF_AA_2STAGE, SSYTRS_AA_2STAGE, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool         LERR, OK;
      String       SRNAMT;
      int          INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA         ISEEDY / 1988, 1989, 1990, 1991 /;
      DATA         UPLOS / 'U', 'L' /;
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.


      // Test path

      PATH( 1: 1 ) = 'Single precision';
      PATH( 2: 3 ) = 'S2';

      // Path to generate matrices

      MATPATH( 1: 1 ) = 'Single precision';
      MATPATH( 2: 3 ) = 'SY';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) CALL SERRSY( PATH, NOUT );
      INFOT = 0;

      // Set the minimum block size for which the block routine should
      // be used, which will be later returned by ILAENV

      xlaenv(2, 2 );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 180
         N = NVAL( IN );
         if ( N > NMAX ) {
            NFAIL = NFAIL + 1;
            WRITE(NOUT, 9995) 'M ', N, NMAX;
            GO TO 180;
         }
         LDA = MAX( N, 1 );
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         IZERO = 0;

         // Do for each value of matrix type IMAT

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 170

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( !DOTYPE( IMAT ) ) GO TO 170;

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 6;
            if (ZEROT && N < IMAT-2) GO TO 170;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               UPLO = UPLOS( IUPLO );

               // Begin generate the test matrix A.


               // Set up parameters with SLATB4 for the matrix generator
               // based on the type of matrix to be generated.

               slatb4(MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               // Generate a matrix with SLATMS.

               SRNAMT = 'SLATMS';
               slatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from SLATMS and handle error.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'SLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     // Skip all tests for this generated matrix

                  GO TO 160;
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
                     if ( IUPLO == 1 ) {

                        // Set the first IZERO rows and columns to zero.

                        IOFF = 0;
                        for (J = 1; J <= N; J++) { // 70
                           I2 = MIN( J, IZERO );
                           for (I = 1; I <= I2; I++) { // 60
                              A( IOFF+I ) = ZERO;
                           } // 60
                           IOFF = IOFF + LDA;
                        } // 70
                        IZERO = 1;
                     } else {

                        // Set the last IZERO rows and columns to zero.

                        IOFF = 0;
                        for (J = 1; J <= N; J++) { // 90
                           I1 = MAX( J, IZERO );
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

               // End generate the test matrix A.

               // Do for each value of NB in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 150

                  // Set the optimal blocksize, which will be later
                  // returned by ILAENV.

                  NB = NBVAL( INB );
                  xlaenv(1, NB );

                  // Copy the test matrix A into matrix AFAC which
                  // will be factorized in place. This is needed to
                  // preserve the test matrix A for subsequent tests.

                  slacpy(UPLO, N, N, A, LDA, AFAC, LDA );

                  // Compute the L*D*L**T or U*D*U**T factorization of the
                  // matrix. IWORK stores details of the interchanges and
                  // the block structure of D. AINV is a work array for
                  // block factorization, LWORK is the length of AINV.

                  SRNAMT = 'SSYTRF_AA_2STAGE';
                  LWORK = MIN( MAX( 1, N*NB ), 3*NMAX*NMAX );
                  ssytrf_aa_2stage(UPLO, N, AFAC, LDA,  AINV, MAX( 1, (3*NB+1)*N ), IWORK, IWORK( 1+N ), WORK, LWORK, INFO );

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                  if ( IZERO > 0 ) {
                     J = 1;
                     K = IZERO;
                     } // 100
                     if ( J == K ) {
                        K = IWORK( J );
                     } else if ( IWORK( J ) == K ) {
                        K = J;
                     }
                     if ( J < K ) {
                        J = J + 1;
                        GO TO 100;
                     }
                  } else {
                     K = 0;
                  }

                  // Check error code from SSYTRF and handle error.

                  if ( INFO != K ) {
                     alaerh(PATH, 'SSYTRF_AA_2STAGE', INFO, K, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );
                  }

// +    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                   // CALL SSYT01_AA( UPLO, N, A, LDA, AFAC, LDA, IWORK,
      // $                            AINV, LDA, RWORK, RESULT( 1 ) )
                   // NT = 1
                   NT = 0;


                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= NT; K++) { // 110
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 110
                  NRUN = NRUN + NT;

                  // Skip solver test if INFO is not 0.

                  if ( INFO != 0 ) {
                     GO TO 140;
                  }

                  // Do for each value of NRHS in NSVAL.

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 130
                     NRHS = NSVAL( IRHS );

// +    TEST 2 (Using TRS)
                  // Solve and compute residual for  A * X = B.

                     // Choose a set of NRHS random solution vectors
                     // stored in XACT and set up the right hand side B

                     SRNAMT = 'SLARHS';
                     slarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     slacpy('Full', N, NRHS, B, LDA, X, LDA );

                     SRNAMT = 'SSYTRS_AA_2STAGE';
                     ssytrs_aa_2stage(UPLO, N, NRHS, AFAC, LDA, AINV, (3*NB+1)*N, IWORK, IWORK( 1+N ), X, LDA, INFO );

                     // Check error code from SSYTRS and handle error.

                     if ( INFO != 0 ) {
                        if ( IZERO == 0 ) {
                           alaerh(PATH, 'SSYTRS_AA_2STAGE', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        }
                     } else {
                        slacpy('Full', N, NRHS, B, LDA, WORK, LDA );

                        // Compute the residual for the solution

                        spot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );


                     // Print information about the tests that did not pass
                     // the threshold.

                        for (K = 2; K <= 2; K++) { // 120
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 120
                     }
                     NRUN = NRUN + 1;

                  // End do for each value of NRHS in NSVAL.

                  } // 130
                  } // 140
               } // 150
            } // 160
         } // 170
      } // 180

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ', I2, ', test ', I2, ', ratio =', G12.5 );
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') =', G12.5 );
 9995 FORMAT( ' Invalid input value: ', A4, '=', I6, '; must be <=', I6 )
      return;

      // End of SCHKSY_AA_2STAGE

      }
