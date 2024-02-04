      void zchkgt(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A, AF, B, X, XACT, WORK, RWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NN, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      Complex         A( * ), AF( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 12 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, NORM, TRANS, TYPE;
      String             PATH;
      int                I, IMAT, IN, INFO, IRHS, ITRAN, IX, IZERO, J, K, KL, KOFF, KU, LDA, M, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      double             AINVNM, ANORM, COND, RCOND, RCONDC, RCONDI, RCONDO;
      // ..
      // .. Local Arrays ..
      String             TRANSS( 3 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      Complex         Z( 3 );
      // ..
      // .. External Functions ..
      //- double             DGET06, DZASUM, ZLANGT;
      // EXTERNAL DGET06, DZASUM, ZLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, ZCOPY, ZDSCAL, ZERRGE, ZGET04, ZGTCON, ZGTRFS, ZGTT01, ZGTT02, ZGTT05, ZGTTRF, ZGTTRS, ZLACPY, ZLAGTM, ZLARNV, ZLATB4, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
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
      const ISEEDY = 0, 0, 0, 1, TRANSS = 'N', 'T', 'C';
      // ..
      // .. Executable Statements ..

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'GT';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) zerrge( PATH, NOUT );
      INFOT = 0;

      for (IN = 1; IN <= NN; IN++) { // 110

         // Do for each value of N in NVAL.

         N = NVAL( IN );
         M = max( N-1, 0 );
         LDA = max( 1, N );
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 100;

            // Set up parameters with ZLATB4.

            zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST );

            ZEROT = IMAT >= 8 && IMAT <= 10;
            if ( IMAT <= 6 ) {

               // Types 1-6:  generate matrices of known condition number.

               KOFF = max( 2-KU, 3-max( 1, N ) );
               SRNAMT = 'ZLATMS';
               zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'Z', AF( KOFF ), 3, WORK, INFO );

               // Check the error code from ZLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'ZLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100;
               }
               IZERO = 0;

               if ( N > 1 ) {
                  zcopy(N-1, AF( 4 ), 3, A, 1 );
                  zcopy(N-1, AF( 3 ), 3, A( N+M+1 ), 1 );
               }
               zcopy(N, AF( 2 ), 3, A( M+1 ), 1 );
            } else {

               // Types 7-12:  generate tridiagonal matrices with
               // unknown condition numbers.

               if ( !ZEROT || !DOTYPE( 7 ) ) {

                  // Generate a matrix with elements whose real and
                  // imaginary parts are from [-1,1].

                  zlarnv(2, ISEED, N+2*M, A );
                  if (ANORM != ONE) zdscal( N+2*M, ANORM, A, 1 );
               } else if ( IZERO > 0 ) {

                  // Reuse the last matrix by copying back the zeroed out
                  // elements.

                  if ( IZERO == 1 ) {
                     A[N] = Z( 2 );
                     if (N > 1) A( 1 ) = Z( 3 );
                  } else if ( IZERO == N ) {
                     A[3*N-2] = Z( 1 );
                     A[2*N-1] = Z( 2 );
                  } else {
                     A[2*N-2+IZERO] = Z( 1 );
                     A[N-1+IZERO] = Z( 2 );
                     A[IZERO] = Z( 3 );
                  }
               }

               // If IMAT > 7, set one column of the matrix to 0.

               if ( !ZEROT ) {
                  IZERO = 0;
               } else if ( IMAT == 8 ) {
                  IZERO = 1;
                  Z[2] = A( N );
                  A[N] = ZERO;
                  if ( N > 1 ) {
                     Z[3] = A( 1 );
                     A[1] = ZERO;
                  }
               } else if ( IMAT == 9 ) {
                  IZERO = N;
                  Z[1] = A( 3*N-2 );
                  Z[2] = A( 2*N-1 );
                  A[3*N-2] = ZERO;
                  A[2*N-1] = ZERO;
               } else {
                  IZERO = ( N+1 ) / 2;
                  for (I = IZERO; I <= N - 1; I++) { // 20
                     A[2*N-2+I] = ZERO;
                     A[N-1+I] = ZERO;
                     A[I] = ZERO;
                  } // 20
                  A[3*N-2] = ZERO;
                  A[2*N-1] = ZERO;
               }
            }

// +    TEST 1
            // Factor A as L*U and compute the ratio
               // norm(L*U - A) / (n * norm(A) * EPS )

            zcopy(N+2*M, A, 1, AF, 1 );
            SRNAMT = 'ZGTTRF';
            zgttrf(N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, INFO );

            // Check error code from ZGTTRF.

            if (INFO != IZERO) alaerh( PATH, 'ZGTTRF', INFO, IZERO, ' ', N, N, 1, 1, -1, IMAT, NFAIL, NERRS, NOUT );
            TRFCON = INFO != 0;

            zgtt01(N, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, WORK, LDA, RWORK, RESULT( 1 ) );

            // Print the test ratio if it is >= THRESH.

            if ( RESULT( 1 ) >= THRESH ) {
               if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
               WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 );
               NFAIL = NFAIL + 1;
            }
            NRUN = NRUN + 1;

            for (ITRAN = 1; ITRAN <= 2; ITRAN++) { // 50
               TRANS = TRANSS( ITRAN );
               if ( ITRAN == 1 ) {
                  NORM = 'O';
               } else {
                  NORM = 'I';
               }
               ANORM = ZLANGT( NORM, N, A, A( M+1 ), A( N+M+1 ) );

               if ( !TRFCON ) {

                  // Use ZGTTRS to solve for one column at a time of
                  // inv(A), computing the maximum column sum as we go.

                  AINVNM = ZERO;
                  for (I = 1; I <= N; I++) { // 40
                     for (J = 1; J <= N; J++) { // 30
                        X[J] = ZERO;
                     } // 30
                     X[I] = ONE;
                     zgttrs(TRANS, N, 1, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO );
                     AINVNM = max( AINVNM, DZASUM( N, X, 1 ) );
                  } // 40

                  // Compute RCONDC = 1 / (norm(A) * norm(inv(A))

                  if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                     RCONDC = ONE;
                  } else {
                     RCONDC = ( ONE / ANORM ) / AINVNM;
                  }
                  if ( ITRAN == 1 ) {
                     RCONDO = RCONDC;
                  } else {
                     RCONDI = RCONDC;
                  }
               } else {
                  RCONDC = ZERO;
               }

// +    TEST 7
               // Estimate the reciprocal of the condition number of the
               // matrix.

               SRNAMT = 'ZGTCON';
               zgtcon(NORM, N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, ANORM, RCOND, WORK, INFO );

               // Check error code from ZGTCON.

               if (INFO != 0) alaerh( PATH, 'ZGTCON', INFO, 0, NORM, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

               RESULT[7] = DGET06( RCOND, RCONDC );

               // Print the test ratio if it is >= THRESH.

               if ( RESULT( 7 ) >= THRESH ) {
                  if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                  WRITE( NOUT, FMT = 9997 )NORM, N, IMAT, 7, RESULT( 7 );
                  NFAIL = NFAIL + 1;
               }
               NRUN = NRUN + 1;
            } // 50

            // Skip the remaining tests if the matrix is singular.

            if (TRFCON) GO TO 100;

            for (IRHS = 1; IRHS <= NNS; IRHS++) { // 90
               NRHS = NSVAL( IRHS );

               // Generate NRHS random solution vectors.

               IX = 1;
               for (J = 1; J <= NRHS; J++) { // 60
                  zlarnv(2, ISEED, N, XACT( IX ) );
                  IX = IX + LDA;
               } // 60

               for (ITRAN = 1; ITRAN <= 3; ITRAN++) { // 80
                  TRANS = TRANSS( ITRAN );
                  if ( ITRAN == 1 ) {
                     RCONDC = RCONDO;
                  } else {
                     RCONDC = RCONDI;
                  }

                  // Set the right hand side.

                  zlagtm(TRANS, N, NRHS, ONE, A, A( M+1 ), A( N+M+1 ), XACT, LDA, ZERO, B, LDA );

// +    TEST 2
               // Solve op(A) * X = B and compute the residual.

                  zlacpy('Full', N, NRHS, B, LDA, X, LDA );
                  SRNAMT = 'ZGTTRS';
                  zgttrs(TRANS, N, NRHS, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO );

               // Check error code from ZGTTRS.

                  if (INFO != 0) alaerh( PATH, 'ZGTTRS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                  zgtt02(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), X, LDA, WORK, LDA, RESULT( 2 ) );

// +    TEST 3
               // Check solution from generated exact solution.

                  zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, and 6
               // Use iterative refinement to improve the solution.

                  SRNAMT = 'ZGTRFS';
                  zgtrfs(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

               // Check error code from ZGTRFS.

                  if (INFO != 0) alaerh( PATH, 'ZGTRFS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                  zgtt05(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

               // Print information about the tests that did not pass the
               // threshold.

                  for (K = 2; K <= 6; K++) { // 70
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 70
                  NRUN = NRUN + 5;
               } // 80
            } // 90
         } // 100
      } // 110

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 12X, 'N =', I5, ',', 10X, ' type ', I2, ', test(', I2, ') = ', G12.5 );
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') = ', G12.5 );
 9997 FORMAT( ' NORM =''', A1, ''', N =', I5, ',', 10X, ' type ', I2, ', test(', I2, ') = ', G12.5 );
      return;
      }
