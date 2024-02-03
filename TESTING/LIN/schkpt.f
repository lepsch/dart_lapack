      SUBROUTINE SCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A, D, E, B, X, XACT, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NN, NNS, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NSVAL( * ), NVAL( * );
      REAL               A( * ), B( * ), D( * ), E( * ), RWORK( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      int                NTYPES;
      const              NTYPES = 12 ;
      int                NTESTS;
      const              NTESTS = 7 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, TYPE;
      String             PATH;
      int                I, IA, IMAT, IN, INFO, IRHS, IX, IZERO, J, K, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      REAL               AINVNM, ANORM, COND, DMAX, RCOND, RCONDC
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS ), Z( 3 )
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SASUM, SGET06, SLANST
      // EXTERNAL ISAMAX, SASUM, SGET06, SLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SCOPY, SERRGT, SGET04, SLACPY, SLAPTM, SLARNV, SLATB4, SLATMS, SPTCON, SPTRFS, SPTT01, SPTT02, SPTT05, SPTTRF, SPTTRS, SSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
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
      DATA               ISEEDY / 0, 0, 0, 1 /
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'PT'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      if (TSTERR) CALL SERRGT( PATH, NOUT );
      INFOT = 0

      for (IN = 1; IN <= NN; IN++) { // 110

         // Do for each value of N in NVAL.

         N = NVAL( IN )
         LDA = MAX( 1, N )
         NIMAT = NTYPES
         if (N.LE.0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( N.GT.0 .AND. .NOT.DOTYPE( IMAT ) ) GO TO 100

            // Set up parameters with SLATB4.

            slatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST );

            ZEROT = IMAT.GE.8 .AND. IMAT.LE.10
            if ( IMAT.LE.6 ) {

               // Type 1-6:  generate a symmetric tridiagonal matrix of
               // known condition number in lower triangular band storage.

               SRNAMT = 'SLATMS'
               slatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'B', A, 2, WORK, INFO );

               // Check the error code from SLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'SLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100
               }
               IZERO = 0

               // Copy the matrix to D and E.

               IA = 1
               for (I = 1; I <= N - 1; I++) { // 20
                  D( I ) = A( IA )
                  E( I ) = A( IA+1 )
                  IA = IA + 2
               } // 20
               if (N.GT.0) D( N ) = A( IA );
            } else {

               // Type 7-12:  generate a diagonally dominant matrix with
               // unknown condition number in the vectors D and E.

               if ( .NOT.ZEROT .OR. .NOT.DOTYPE( 7 ) ) {

                  // Let D and E have values from [-1,1].

                  slarnv(2, ISEED, N, D );
                  slarnv(2, ISEED, N-1, E );

                  // Make the tridiagonal matrix diagonally dominant.

                  if ( N == 1 ) {
                     D( 1 ) = ABS( D( 1 ) )
                  } else {
                     D( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                     D( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
                     for (I = 2; I <= N - 1; I++) { // 30
                        D( I ) = ABS( D( I ) ) + ABS( E( I ) ) + ABS( E( I-1 ) )
                     } // 30
                  }

                  // Scale D and E so the maximum element is ANORM.

                  IX = ISAMAX( N, D, 1 )
                  DMAX = D( IX )
                  sscal(N, ANORM / DMAX, D, 1 );
                  sscal(N-1, ANORM / DMAX, E, 1 );

               } else if ( IZERO.GT.0 ) {

                  // Reuse the last matrix by copying back the zeroed out
                  // elements.

                  if ( IZERO == 1 ) {
                     D( 1 ) = Z( 2 )
                     if (N.GT.1) E( 1 ) = Z( 3 );
                  } else if ( IZERO == N ) {
                     E( N-1 ) = Z( 1 )
                     D( N ) = Z( 2 )
                  } else {
                     E( IZERO-1 ) = Z( 1 )
                     D( IZERO ) = Z( 2 )
                     E( IZERO ) = Z( 3 )
                  }
               }

               // For types 8-10, set one row and column of the matrix to
               // zero.

               IZERO = 0
               if ( IMAT == 8 ) {
                  IZERO = 1
                  Z( 2 ) = D( 1 )
                  D( 1 ) = ZERO
                  if ( N.GT.1 ) {
                     Z( 3 ) = E( 1 )
                     E( 1 ) = ZERO
                  }
               } else if ( IMAT == 9 ) {
                  IZERO = N
                  if ( N.GT.1 ) {
                     Z( 1 ) = E( N-1 )
                     E( N-1 ) = ZERO
                  }
                  Z( 2 ) = D( N )
                  D( N ) = ZERO
               } else if ( IMAT == 10 ) {
                  IZERO = ( N+1 ) / 2
                  if ( IZERO.GT.1 ) {
                     Z( 1 ) = E( IZERO-1 )
                     E( IZERO-1 ) = ZERO
                     Z( 3 ) = E( IZERO )
                     E( IZERO ) = ZERO
                  }
                  Z( 2 ) = D( IZERO )
                  D( IZERO ) = ZERO
               }
            }

            scopy(N, D, 1, D( N+1 ), 1 );
            if (N.GT.1) CALL SCOPY( N-1, E, 1, E( N+1 ), 1 );

*+    TEST 1
            // Factor A as L*D*L' and compute the ratio
               // norm(L*D*L' - A) / (n * norm(A) * EPS )

            spttrf(N, D( N+1 ), E( N+1 ), INFO );

            // Check error code from SPTTRF.

            if ( INFO != IZERO ) {
               alaerh(PATH, 'SPTTRF', INFO, IZERO, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               GO TO 100
            }

            if ( INFO.GT.0 ) {
               RCONDC = ZERO
               GO TO 90
            }

            sptt01(N, D, E, D( N+1 ), E( N+1 ), WORK, RESULT( 1 ) );

            // Print the test ratio if greater than or equal to THRESH.

            if ( RESULT( 1 ).GE.THRESH ) {
               if (NFAIL == 0 .AND. NERRS == 0) CALL ALAHD( NOUT, PATH );
               WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 )
               NFAIL = NFAIL + 1
            }
            NRUN = NRUN + 1

            // Compute RCONDC = 1 / (norm(A) * norm(inv(A))

            // Compute norm(A).

            ANORM = SLANST( '1', N, D, E )

            // Use SPTTRS to solve for one column at a time of inv(A),
            // computing the maximum column sum as we go.

            AINVNM = ZERO
            for (I = 1; I <= N; I++) { // 50
               for (J = 1; J <= N; J++) { // 40
                  X( J ) = ZERO
               } // 40
               X( I ) = ONE
               spttrs(N, 1, D( N+1 ), E( N+1 ), X, LDA, INFO );
               AINVNM = MAX( AINVNM, SASUM( N, X, 1 ) )
            } // 50
            RCONDC = ONE / MAX( ONE, ANORM*AINVNM )

            for (IRHS = 1; IRHS <= NNS; IRHS++) { // 80
               NRHS = NSVAL( IRHS )

            // Generate NRHS random solution vectors.

               IX = 1
               for (J = 1; J <= NRHS; J++) { // 60
                  slarnv(2, ISEED, N, XACT( IX ) );
                  IX = IX + LDA
               } // 60

            // Set the right hand side.

               slaptm(N, NRHS, ONE, D, E, XACT, LDA, ZERO, B, LDA );

*+    TEST 2
            // Solve A*x = b and compute the residual.

               slacpy('Full', N, NRHS, B, LDA, X, LDA );
               spttrs(N, NRHS, D( N+1 ), E( N+1 ), X, LDA, INFO );

            // Check error code from SPTTRS.

               if (INFO != 0) CALL ALAERH( PATH, 'SPTTRS', INFO, 0, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

               slacpy('Full', N, NRHS, B, LDA, WORK, LDA );
               sptt02(N, NRHS, D, E, X, LDA, WORK, LDA, RESULT( 2 ) );

*+    TEST 3
            // Check solution from generated exact solution.

               sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

*+    TESTS 4, 5, and 6
            // Use iterative refinement to improve the solution.

               SRNAMT = 'SPTRFS'
               sptrfs(N, NRHS, D, E, D( N+1 ), E( N+1 ), B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, INFO );

            // Check error code from SPTRFS.

               if (INFO != 0) CALL ALAERH( PATH, 'SPTRFS', INFO, 0, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

               sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
               sptt05(N, NRHS, D, E, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

            // Print information about the tests that did not pass the
            // threshold.

               for (K = 2; K <= 6; K++) { // 70
                  if ( RESULT( K ).GE.THRESH ) {
                     if (NFAIL == 0 .AND. NERRS == 0) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9998 )N, NRHS, IMAT, K, RESULT( K );
                     NFAIL = NFAIL + 1
                  }
               } // 70
               NRUN = NRUN + 5
            } // 80

*+    TEST 7
            // Estimate the reciprocal of the condition number of the
            // matrix.

            } // 90
            SRNAMT = 'SPTCON'
            sptcon(N, D( N+1 ), E( N+1 ), ANORM, RCOND, RWORK, INFO );

            // Check error code from SPTCON.

            if (INFO != 0) CALL ALAERH( PATH, 'SPTCON', INFO, 0, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

            RESULT( 7 ) = SGET06( RCOND, RCONDC )

            // Print the test ratio if greater than or equal to THRESH.

            if ( RESULT( 7 ).GE.THRESH ) {
               if (NFAIL == 0 .AND. NERRS == 0) CALL ALAHD( NOUT, PATH );
               WRITE( NOUT, FMT = 9999 )N, IMAT, 7, RESULT( 7 )
               NFAIL = NFAIL + 1
            }
            NRUN = NRUN + 1
         } // 100
      } // 110

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' N =', I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
 9998 FORMAT( ' N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') = ', G12.5 )
      RETURN

      // End of SCHKPT

      }
