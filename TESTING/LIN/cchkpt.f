      SUBROUTINE CCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A, D, E, B, X, XACT, WORK, RWORK, NOUT )

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
      REAL               D( * ), RWORK( * )
      COMPLEX            A( * ), B( * ), E( * ), WORK( * ), X( * ), XACT( * )
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
      String             DIST, TYPE, UPLO;
      String             PATH;
      int                I, IA, IMAT, IN, INFO, IRHS, IUPLO, IX, IZERO, J, K, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      REAL               AINVNM, ANORM, COND, DMAX, RCOND, RCONDC
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS )
      COMPLEX            Z( 3 )
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               CLANHT, SCASUM, SGET06
      // EXTERNAL ISAMAX, CLANHT, SCASUM, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CCOPY, CERRGT, CGET04, CLACPY, CLAPTM, CLARNV, CLATB4, CLATMS, CPTCON, CPTRFS, CPTT01, CPTT02, CPTT05, CPTTRF, CPTTRS, CSSCAL, SCOPY, SLARNV, SSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL
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
      DATA               ISEEDY / 0, 0, 0, 1 / , UPLOS / 'U', 'L' /
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'PT'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      if (TSTERR) CALL CERRGT( PATH, NOUT );
      INFOT = 0

      for (IN = 1; IN <= NN; IN++) { // 120

         // Do for each value of N in NVAL.

         N = NVAL( IN )
         LDA = MAX( 1, N )
         NIMAT = NTYPES
         if (N.LE.0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 110

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( N > 0 && .NOT.DOTYPE( IMAT ) ) GO TO 110

            // Set up parameters with CLATB4.

            clatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST );

            ZEROT = IMAT >= 8 && IMAT.LE.10
            if ( IMAT.LE.6 ) {

               // Type 1-6:  generate a Hermitian tridiagonal matrix of
               // known condition number in lower triangular band storage.

               SRNAMT = 'CLATMS'
               clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'B', A, 2, WORK, INFO );

               // Check the error code from CLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'CLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 110
               }
               IZERO = 0

               // Copy the matrix to D and E.

               IA = 1
               for (I = 1; I <= N - 1; I++) { // 20
                  D( I ) = REAL( A( IA ) )
                  E( I ) = A( IA+1 )
                  IA = IA + 2
               } // 20
               if (N > 0) D( N ) = REAL( A( IA ) );
            } else {

               // Type 7-12:  generate a diagonally dominant matrix with
               // unknown condition number in the vectors D and E.

               if ( .NOT.ZEROT || .NOT.DOTYPE( 7 ) ) {

                  // Let E be complex, D real, with values from [-1,1].

                  slarnv(2, ISEED, N, D );
                  clarnv(2, ISEED, N-1, E );

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
                  csscal(N-1, ANORM / DMAX, E, 1 );

               } else if ( IZERO > 0 ) {

                  // Reuse the last matrix by copying back the zeroed out
                  // elements.

                  if ( IZERO == 1 ) {
                     D( 1 ) = REAL( Z( 2 ) )
                     if (N > 1) E( 1 ) = Z( 3 );
                  } else if ( IZERO == N ) {
                     E( N-1 ) = Z( 1 )
                     D( N ) = REAL( Z( 2 ) )
                  } else {
                     E( IZERO-1 ) = Z( 1 )
                     D( IZERO ) = REAL( Z( 2 ) )
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
                  if ( N > 1 ) {
                     Z( 3 ) = E( 1 )
                     E( 1 ) = ZERO
                  }
               } else if ( IMAT == 9 ) {
                  IZERO = N
                  if ( N > 1 ) {
                     Z( 1 ) = E( N-1 )
                     E( N-1 ) = ZERO
                  }
                  Z( 2 ) = D( N )
                  D( N ) = ZERO
               } else if ( IMAT == 10 ) {
                  IZERO = ( N+1 ) / 2
                  if ( IZERO > 1 ) {
                     Z( 1 ) = E( IZERO-1 )
                     Z( 3 ) = E( IZERO )
                     E( IZERO-1 ) = ZERO
                     E( IZERO ) = ZERO
                  }
                  Z( 2 ) = D( IZERO )
                  D( IZERO ) = ZERO
               }
            }

            scopy(N, D, 1, D( N+1 ), 1 );
            if (N > 1) CALL CCOPY( N-1, E, 1, E( N+1 ), 1 );

*+    TEST 1
            // Factor A as L*D*L' and compute the ratio
               // norm(L*D*L' - A) / (n * norm(A) * EPS )

            cpttrf(N, D( N+1 ), E( N+1 ), INFO );

            // Check error code from CPTTRF.

            if ( INFO != IZERO ) {
               alaerh(PATH, 'CPTTRF', INFO, IZERO, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               GO TO 110
            }

            if ( INFO > 0 ) {
               RCONDC = ZERO
               GO TO 100
            }

            cptt01(N, D, E, D( N+1 ), E( N+1 ), WORK, RESULT( 1 ) );

            // Print the test ratio if greater than or equal to THRESH.

            if ( RESULT( 1 ) >= THRESH ) {
               if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH );
               WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 )
               NFAIL = NFAIL + 1
            }
            NRUN = NRUN + 1

            // Compute RCONDC = 1 / (norm(A) * norm(inv(A))

            // Compute norm(A).

            ANORM = CLANHT( '1', N, D, E )

            // Use CPTTRS to solve for one column at a time of inv(A),
            // computing the maximum column sum as we go.

            AINVNM = ZERO
            for (I = 1; I <= N; I++) { // 50
               for (J = 1; J <= N; J++) { // 40
                  X( J ) = ZERO
               } // 40
               X( I ) = ONE
               cpttrs('Lower', N, 1, D( N+1 ), E( N+1 ), X, LDA, INFO );
               AINVNM = MAX( AINVNM, SCASUM( N, X, 1 ) )
            } // 50
            RCONDC = ONE / MAX( ONE, ANORM*AINVNM )

            for (IRHS = 1; IRHS <= NNS; IRHS++) { // 90
               NRHS = NSVAL( IRHS )

            // Generate NRHS random solution vectors.

               IX = 1
               for (J = 1; J <= NRHS; J++) { // 60
                  clarnv(2, ISEED, N, XACT( IX ) );
                  IX = IX + LDA
               } // 60

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 80

               // Do first for UPLO = 'U', then for UPLO = 'L'.

                  UPLO = UPLOS( IUPLO )

               // Set the right hand side.

                  claptm(UPLO, N, NRHS, ONE, D, E, XACT, LDA, ZERO, B, LDA );

*+    TEST 2
               // Solve A*x = b and compute the residual.

                  clacpy('Full', N, NRHS, B, LDA, X, LDA );
                  cpttrs(UPLO, N, NRHS, D( N+1 ), E( N+1 ), X, LDA, INFO );

               // Check error code from CPTTRS.

                  if (INFO != 0) CALL ALAERH( PATH, 'CPTTRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  clacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                  cptt02(UPLO, N, NRHS, D, E, X, LDA, WORK, LDA, RESULT( 2 ) );

*+    TEST 3
               // Check solution from generated exact solution.

                  cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

*+    TESTS 4, 5, and 6
               // Use iterative refinement to improve the solution.

                  SRNAMT = 'CPTRFS'
                  cptrfs(UPLO, N, NRHS, D, E, D( N+1 ), E( N+1 ), B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

               // Check error code from CPTRFS.

                  if (INFO != 0) CALL ALAERH( PATH, 'CPTRFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                  cptt05(N, NRHS, D, E, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

               // Print information about the tests that did not pass the
               // threshold.

                  for (K = 2; K <= 6; K++) { // 70
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1
                     }
                  } // 70
                  NRUN = NRUN + 5

               } // 80
            } // 90

*+    TEST 7
            // Estimate the reciprocal of the condition number of the
            // matrix.

            } // 100
            SRNAMT = 'CPTCON'
            cptcon(N, D( N+1 ), E( N+1 ), ANORM, RCOND, RWORK, INFO );

            // Check error code from CPTCON.

            if (INFO != 0) CALL ALAERH( PATH, 'CPTCON', INFO, 0, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

            RESULT( 7 ) = SGET06( RCOND, RCONDC )

            // Print the test ratio if greater than or equal to THRESH.

            if ( RESULT( 7 ) >= THRESH ) {
               if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH );
               WRITE( NOUT, FMT = 9999 )N, IMAT, 7, RESULT( 7 )
               NFAIL = NFAIL + 1
            }
            NRUN = NRUN + 1
         } // 110
      } // 120

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' N =', I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS =', I3, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
      RETURN

      // End of CCHKPT

      }
