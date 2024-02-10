import 'common.dart';

      void dchkpt(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A, D, E, B, X, XACT, final Array<double> _WORK, final Array<double> RWORK, final int NOUT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NN, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                NSVAL( * ), NVAL( * );
      double             A( * ), B( * ), D( * ), E( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 12 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      bool               ZEROT;
      String             DIST, TYPE;
      String             PATH;
      int                I, IA, IMAT, IN, INFO, IRHS, IX, IZERO, J, K, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      double             AINVNM, ANORM, COND, DMAX, RCOND, RCONDC;
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS ), Z( 3 );
      // ..
      // .. External Functions ..
      //- int                idamax;
      //- double             DASUM, DGET06, DLANST;
      // EXTERNAL idamax, DASUM, DGET06, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DCOPY, DERRGT, DGET04, DLACPY, DLAPTM, DLARNV, DLATB4, DLATMS, DPTCON, DPTRFS, DPTT01, DPTT02, DPTT05, DPTTRF, DPTTRS, DSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
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
      const ISEEDY = [ 0, 0, 0, 1 ];

      PATH = '${'Double precision'[0]}';
      PATH[2: 3] = 'PT';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) derrgt( PATH, NOUT );
      infoc.INFOT = 0;

      for (IN = 1; IN <= NN; IN++) { // 110

         // Do for each value of N in NVAL.

         N = NVAL( IN );
         LDA = max( 1, N );
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( N > 0 && !DOTYPE( IMAT ) ) GO TO 100;

            // Set up parameters with DLATB4.

            dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST );

            ZEROT = IMAT >= 8 && IMAT <= 10;
            if ( IMAT <= 6 ) {

               // Type 1-6:  generate a symmetric tridiagonal matrix of
               // known condition number in lower triangular band storage.

               srnamc.SRNAMT = 'DLATMS';
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'B', A, 2, WORK, INFO );

               // Check the error code from DLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100;
               }
               IZERO = 0;

               // Copy the matrix to D and E.

               IA = 1;
               for (I = 1; I <= N - 1; I++) { // 20
                  D[I] = A( IA );
                  E[I] = A( IA+1 );
                  IA = IA + 2;
               } // 20
               if (N > 0) D( N ) = A( IA );
            } else {

               // Type 7-12:  generate a diagonally dominant matrix with
               // unknown condition number in the vectors D and E.

               if ( !ZEROT || !DOTYPE( 7 ) ) {

                  // Let D and E have values from [-1,1].

                  dlarnv(2, ISEED, N, D );
                  dlarnv(2, ISEED, N-1, E );

                  // Make the tridiagonal matrix diagonally dominant.

                  if ( N == 1 ) {
                     D[1] = ( D( 1 ) ).abs();
                  } else {
                     D[1] = ( D( 1 ) ).abs() + ( E( 1 ) ).abs();
                     D[N] = ( D( N ) ).abs() + ( E( N-1 ) ).abs();
                     for (I = 2; I <= N - 1; I++) { // 30
                        D[I] = ( D( I ) ).abs() + ( E( I ) ).abs() + ( E( I-1 ) ).abs();
                     } // 30
                  }

                  // Scale D and E so the maximum element is ANORM.

                  IX = idamax( N, D, 1 );
                  DMAX = D( IX );
                  dscal(N, ANORM / DMAX, D, 1 );
                  dscal(N-1, ANORM / DMAX, E, 1 );

               } else if ( IZERO > 0 ) {

                  // Reuse the last matrix by copying back the zeroed out
                  // elements.

                  if ( IZERO == 1 ) {
                     D[1] = Z( 2 );
                     if (N > 1) E( 1 ) = Z( 3 );
                  } else if ( IZERO == N ) {
                     E[N-1] = Z( 1 );
                     D[N] = Z( 2 );
                  } else {
                     E[IZERO-1] = Z( 1 );
                     D[IZERO] = Z( 2 );
                     E[IZERO] = Z( 3 );
                  }
               }

               // For types 8-10, set one row and column of the matrix to
               // zero.

               IZERO = 0;
               if ( IMAT == 8 ) {
                  IZERO = 1;
                  Z[2] = D( 1 );
                  D[1] = ZERO;
                  if ( N > 1 ) {
                     Z[3] = E( 1 );
                     E[1] = ZERO;
                  }
               } else if ( IMAT == 9 ) {
                  IZERO = N;
                  if ( N > 1 ) {
                     Z[1] = E( N-1 );
                     E[N-1] = ZERO;
                  }
                  Z[2] = D( N );
                  D[N] = ZERO;
               } else if ( IMAT == 10 ) {
                  IZERO = ( N+1 ) / 2;
                  if ( IZERO > 1 ) {
                     Z[1] = E( IZERO-1 );
                     E[IZERO-1] = ZERO;
                     Z[3] = E( IZERO );
                     E[IZERO] = ZERO;
                  }
                  Z[2] = D( IZERO );
                  D[IZERO] = ZERO;
               }
            }

            dcopy(N, D, 1, D( N+1 ), 1 );
            if (N > 1) dcopy( N-1, E, 1, E( N+1 ), 1 );

// +    TEST 1
            // Factor A as L*D*L' and compute the ratio
            //    norm(L*D*L' - A) / (n * norm(A) * EPS )

            dpttrf(N, D( N+1 ), E( N+1 ), INFO );

            // Check error code from DPTTRF.

            if ( INFO != IZERO ) {
               alaerh(PATH, 'DPTTRF', INFO, IZERO, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               GO TO 100;
            }

            if ( INFO > 0 ) {
               RCONDC = ZERO;
               GO TO 90;
            }

            dptt01(N, D, E, D( N+1 ), E( N+1 ), WORK, RESULT( 1 ) );

            // Print the test ratio if greater than or equal to THRESH.

            if ( RESULT( 1 ) >= THRESH ) {
               if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
               WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 );
               NFAIL = NFAIL + 1;
            }
            NRUN = NRUN + 1;

            // Compute RCONDC = 1 / (norm(A) * norm(inv(A))

            // Compute norm(A).

            ANORM = dlanst( '1', N, D, E );

            // Use DPTTRS to solve for one column at a time of inv(A),
            // computing the maximum column sum as we go.

            AINVNM = ZERO;
            for (I = 1; I <= N; I++) { // 50
               for (J = 1; J <= N; J++) { // 40
                  X[J] = ZERO;
               } // 40
               X[I] = ONE;
               dpttrs(N, 1, D( N+1 ), E( N+1 ), X, LDA, INFO );
               AINVNM = max( AINVNM, dasum( N, X, 1 ) );
            } // 50
            RCONDC = ONE / max( ONE, ANORM*AINVNM );

            for (IRHS = 1; IRHS <= NNS; IRHS++) { // 80
               NRHS = NSVAL( IRHS );

            // Generate NRHS random solution vectors.

               IX = 1;
               for (J = 1; J <= NRHS; J++) { // 60
                  dlarnv(2, ISEED, N, XACT( IX ) );
                  IX = IX + LDA;
               } // 60

            // Set the right hand side.

               dlaptm(N, NRHS, ONE, D, E, XACT, LDA, ZERO, B, LDA );

// +    TEST 2
            // Solve A*x = b and compute the residual.

               dlacpy('Full', N, NRHS, B, LDA, X, LDA );
               dpttrs(N, NRHS, D( N+1 ), E( N+1 ), X, LDA, INFO );

            // Check error code from DPTTRS.

               if (INFO != 0) alaerh( PATH, 'DPTTRS', INFO, 0, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

               dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
               dptt02(N, NRHS, D, E, X, LDA, WORK, LDA, RESULT( 2 ) );

// +    TEST 3
            // Check solution from generated exact solution.

               dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, and 6
            // Use iterative refinement to improve the solution.

               srnamc.SRNAMT = 'DPTRFS';
               dptrfs(N, NRHS, D, E, D( N+1 ), E( N+1 ), B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, INFO );

            // Check error code from DPTRFS.

               if (INFO != 0) alaerh( PATH, 'DPTRFS', INFO, 0, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

               dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
               dptt05(N, NRHS, D, E, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

            // Print information about the tests that did not pass the
            // threshold.

               for (K = 2; K <= 6; K++) { // 70
                  if ( RESULT( K ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9998 )N, NRHS, IMAT, K, RESULT( K );
                     NFAIL = NFAIL + 1;
                  }
               } // 70
               NRUN = NRUN + 5;
            } // 80

// +    TEST 7
            // Estimate the reciprocal of the condition number of the
            // matrix.

            } // 90
            srnamc.SRNAMT = 'DPTCON';
            dptcon(N, D( N+1 ), E( N+1 ), ANORM, RCOND, RWORK, INFO );

            // Check error code from DPTCON.

            if (INFO != 0) alaerh( PATH, 'DPTCON', INFO, 0, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

            RESULT[7] = DGET06( RCOND, RCONDC );

            // Print the test ratio if greater than or equal to THRESH.

            if ( RESULT( 7 ) >= THRESH ) {
               if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
               WRITE( NOUT, FMT = 9999 )N, IMAT, 7, RESULT( 7 );
               NFAIL = NFAIL + 1;
            }
            NRUN = NRUN + 1;
         } // 100
      } // 110

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' N =${.i5}, type ${.i2}, test ${.i2}, ratio = ${.g12_5};
 9998 FORMAT( ' N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) = ${.g12_5};
      }
