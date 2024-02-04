      void cchkhp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      REAL               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NSVAL( * ), NVAL( * );
      REAL               RWORK( * );
      COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 10 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IMAT, IN, INFO, IOFF, IRHS, IUPLO, IZERO, J, K, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NPP, NRHS, NRUN, NT;
      REAL               ANORM, CNDNUM, RCOND, RCONDC;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- REAL               CLANHP, SGET06;
      // EXTERNAL LSAME, CLANHP, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CCOPY, CERRSY, CGET04, CHPCON, CHPRFS, CHPT01, CHPTRF, CHPTRI, CHPTRS, CLACPY, CLAIPD, CLARHS, CLATB4, CLATMS, CPPT02, CPPT03, CPPT05
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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

      PATH[1: 1] = 'Complex precision';
      PATH[2: 3] = 'HP';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) cerrsy( PATH, NOUT );
      INFOT = 0;

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 170
         N = NVAL( IN );
         LDA = max( N, 1 );
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         IZERO = 0;
         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 160

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 160;

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 6;
            if (ZEROT && N < IMAT-2) GO TO 160;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 150
               UPLO = UPLOS( IUPLO );
               if ( LSAME( UPLO, 'U' ) ) {
                  PACKIT = 'C';
               } else {
                  PACKIT = 'R';
               }

               // Set up parameters with CLATB4 and generate a test matrix
               // with CLATMS.

               clatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'CLATMS';
               clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO );

               // Check error code from CLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 150;
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
                     IOFF = 0;
                     if ( IUPLO == 1 ) {

                        // Set the first IZERO rows and columns to zero.

                        for (J = 1; J <= N; J++) { // 70
                           I2 = min( J, IZERO );
                           for (I = 1; I <= I2; I++) { // 60
                              A[IOFF+I] = ZERO;
                           } // 60
                           IOFF = IOFF + J;
                        } // 70
                     } else {

                        // Set the last IZERO rows and columns to zero.

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

               // Set the imaginary part of the diagonals.

               if ( IUPLO == 1 ) {
                  claipd(N, A, 2, 1 );
               } else {
                  claipd(N, A, N, -1 );
               }

               // Compute the L*D*L' or U*D*U' factorization of the matrix.

               NPP = N*( N+1 ) / 2;
               ccopy(NPP, A, 1, AFAC, 1 );
               SRNAMT = 'CHPTRF';
               chptrf(UPLO, N, AFAC, IWORK, INFO );

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

               // Check error code from CHPTRF.

               if (INFO != K) alaerh( PATH, 'CHPTRF', INFO, K, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               if ( INFO != 0 ) {
                  TRFCON = true;
               } else {
                  TRFCON = false;
               }

// +    TEST 1
               // Reconstruct matrix from factors and compute residual.

               chpt01(UPLO, N, A, AFAC, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );
               NT = 1;

// +    TEST 2
               // Form the inverse and compute the residual.

               if ( !TRFCON ) {
                  ccopy(NPP, AFAC, 1, AINV, 1 );
                  SRNAMT = 'CHPTRI';
                  chptri(UPLO, N, AINV, IWORK, WORK, INFO );

               // Check error code from CHPTRI.

                  if (INFO != 0) alaerh( PATH, 'CHPTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  cppt03(UPLO, N, A, AINV, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );
                  NT = 2;
               }

               // Print information about the tests that did not pass
               // the threshold.

               for (K = 1; K <= NT; K++) { // 110
                  if ( RESULT( K ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, K, RESULT( K );
                     NFAIL = NFAIL + 1;
                  }
               } // 110
               NRUN = NRUN + NT;

               // Do only the condition estimate if INFO is not 0.

               if ( TRFCON ) {
                  RCONDC = ZERO;
                  GO TO 140;
               }

               for (IRHS = 1; IRHS <= NNS; IRHS++) { // 130
                  NRHS = NSVAL( IRHS );

// +    TEST 3
               // Solve and compute residual for  A * X = B.

                  SRNAMT = 'CLARHS';
                  clarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  XTYPE = 'C';
                  clacpy('Full', N, NRHS, B, LDA, X, LDA );

                  SRNAMT = 'CHPTRS';
                  chptrs(UPLO, N, NRHS, AFAC, IWORK, X, LDA, INFO );

               // Check error code from CHPTRS.

                  if (INFO != 0) alaerh( PATH, 'CHPTRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  clacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                  cppt02(UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

// +    TEST 4
               // Check solution from generated exact solution.

                  cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );

// +    TESTS 5, 6, and 7
               // Use iterative refinement to improve the solution.

                  SRNAMT = 'CHPRFS';
                  chprfs(UPLO, N, NRHS, A, AFAC, IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

               // Check error code from CHPRFS.

                  if (INFO != 0) alaerh( PATH, 'CHPRFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) );
                  cppt05(UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 6 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 3; K <= 7; K++) { // 120
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 120
                  NRUN = NRUN + 5;
               } // 130

// +    TEST 8
               // Get an estimate of RCOND = 1/CNDNUM.

               } // 140
               ANORM = CLANHP( '1', UPLO, N, A, RWORK );
               SRNAMT = 'CHPCON';
               chpcon(UPLO, N, AFAC, IWORK, ANORM, RCOND, WORK, INFO );

               // Check error code from CHPCON.

               if (INFO != 0) alaerh( PATH, 'CHPCON', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

               RESULT[8] = SGET06( RCOND, RCONDC );

               // Print the test ratio if it is >= THRESH.

               if ( RESULT( 8 ) >= THRESH ) {
                  if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                  WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, 8, RESULT( 8 );
                  NFAIL = NFAIL + 1;
               }
               NRUN = NRUN + 1;
            } // 150
         } // 160
      } // 170

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 );
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') =', G12.5 );
      return;
      }
