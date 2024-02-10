      void schkgt(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A, AF, B, X, XACT, WORK, final Array<double> RWORK, final Array<int> IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NN, NNS, NOUT;
      double               THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NSVAL( * ), NVAL( * );
      double               A( * ), AF( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 12 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      bool               TRFCON, ZEROT;
      String             DIST, NORM, TRANS, TYPE;
      String             PATH;
      int                I, IMAT, IN, INFO, IRHS, ITRAN, IX, IZERO, J, K, KL, KOFF, KU, LDA, M, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      double               AINVNM, ANORM, COND, RCOND, RCONDC, RCONDI, RCONDO;
      String             TRANSS( 3 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS ), Z( 3 );
      // ..
      // .. External Functions ..
      //- REAL               SASUM, SGET06, SLANGT;
      // EXTERNAL SASUM, SGET06, SLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SCOPY, SERRGE, SGET04, SGTCON, SGTRFS, SGTT01, SGTT02, SGTT05, SGTTRF, SGTTRS, SLACPY, SLAGTM, SLARNV, SLATB4, SLATMS, SSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = 0, 0, 0, 1, TRANSS = 'N', 'T', 'C';

      PATH[1: 1] = 'Single precision';
      PATH[2: 3] = 'GT';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) serrge( PATH, NOUT );
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

            // Set up parameters with SLATB4.

            slatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST );

            ZEROT = IMAT >= 8 && IMAT <= 10;
            if ( IMAT <= 6 ) {

               // Types 1-6:  generate matrices of known condition number.

               KOFF = max( 2-KU, 3-max( 1, N ) );
              srnamc.SRNAMT = 'SLATMS';
               slatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'Z', AF( KOFF ), 3, WORK, INFO );

               // Check the error code from SLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'SLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100;
               }
               IZERO = 0;

               if ( N > 1 ) {
                  scopy(N-1, AF( 4 ), 3, A, 1 );
                  scopy(N-1, AF( 3 ), 3, A( N+M+1 ), 1 );
               }
               scopy(N, AF( 2 ), 3, A( M+1 ), 1 );
            } else {

               // Types 7-12:  generate tridiagonal matrices with
               // unknown condition numbers.

               if ( !ZEROT || !DOTYPE( 7 ) ) {

                  // Generate a matrix with elements from [-1,1].

                  slarnv(2, ISEED, N+2*M, A );
                  if (ANORM != ONE) sscal( N+2*M, ANORM, A, 1 );
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
            //    norm(L*U - A) / (n * norm(A) * EPS )

            scopy(N+2*M, A, 1, AF, 1 );
           srnamc.SRNAMT = 'SGTTRF';
            sgttrf(N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, INFO );

            // Check error code from SGTTRF.

            if (INFO != IZERO) alaerh( PATH, 'SGTTRF', INFO, IZERO, ' ', N, N, 1, 1, -1, IMAT, NFAIL, NERRS, NOUT );
            TRFCON = INFO != 0;

            sgtt01(N, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, WORK, LDA, RWORK, RESULT( 1 ) );

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
               ANORM = SLANGT( NORM, N, A, A( M+1 ), A( N+M+1 ) );

               if ( !TRFCON ) {

                  // Use SGTTRS to solve for one column at a time of inv(A)
                  // or inv(A^T), computing the maximum column sum as we
                  // go.

                  AINVNM = ZERO;
                  for (I = 1; I <= N; I++) { // 40
                     for (J = 1; J <= N; J++) { // 30
                        X[J] = ZERO;
                     } // 30
                     X[I] = ONE;
                     sgttrs(TRANS, N, 1, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO );
                     AINVNM = max( AINVNM, SASUM( N, X, 1 ) );
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

              srnamc.SRNAMT = 'SGTCON';
               sgtcon(NORM, N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, ANORM, RCOND, WORK, IWORK( N+1 ), INFO );

               // Check error code from SGTCON.

               if (INFO != 0) alaerh( PATH, 'SGTCON', INFO, 0, NORM, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

               RESULT[7] = SGET06( RCOND, RCONDC );

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
                  slarnv(2, ISEED, N, XACT( IX ) );
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

                  slagtm(TRANS, N, NRHS, ONE, A, A( M+1 ), A( N+M+1 ), XACT, LDA, ZERO, B, LDA );

// +    TEST 2
                  // Solve op(A) * X = B and compute the residual.

                  slacpy('Full', N, NRHS, B, LDA, X, LDA );
                 srnamc.SRNAMT = 'SGTTRS';
                  sgttrs(TRANS, N, NRHS, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO );

                  // Check error code from SGTTRS.

                  if (INFO != 0) alaerh( PATH, 'SGTTRS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  slacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                  sgtt02(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), X, LDA, WORK, LDA, RESULT( 2 ) );

// +    TEST 3
                  // Check solution from generated exact solution.

                  sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, and 6
                  // Use iterative refinement to improve the solution.

                 srnamc.SRNAMT = 'SGTRFS';
                  sgtrfs(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO );

                  // Check error code from SGTRFS.

                  if (INFO != 0) alaerh( PATH, 'SGTRFS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                  sgtt05(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

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

 9999 FORMAT('${' ' * 12}N =${.i5},${' ' * 10} type ${.i2}, test(${.i2}) = ${.g12_5};
 9998 FORMAT( ' TRANS=''${.a1}'', N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) = ${.g12_5};
 9997 FORMAT( ' NORM =''${.a1}'', N =${.i5},${' ' * 10} type ${.i2}, test(${.i2}) = ${.g12_5};
      }
