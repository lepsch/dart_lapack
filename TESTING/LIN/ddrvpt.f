      SUBROUTINE DDRVPT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, D, E, B, X, XACT, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NVAL( * );
      double             A( * ), B( * ), D( * ), E( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      int                NTYPES;
      const              NTYPES = 12 ;
      int                NTESTS;
      const              NTESTS = 6 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, FACT, TYPE;
      String             PATH;
      int                I, IA, IFACT, IMAT, IN, INFO, IX, IZERO, J, K, K1, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NRUN, NT;
      double             AINVNM, ANORM, COND, DMAX, RCOND, RCONDC;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS ), Z( 3 );
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DASUM, DGET06, DLANST;
      // EXTERNAL IDAMAX, DASUM, DGET06, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DCOPY, DERRVX, DGET04, DLACPY, DLAPTM, DLARNV, DLASET, DLATB4, DLATMS, DPTSV, DPTSVX, DPTT01, DPTT02, DPTT05, DPTTRF, DPTTRS, DSCAL
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

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'PT'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      if (TSTERR) CALL DERRVX( PATH, NOUT );
      INFOT = 0

      for (IN = 1; IN <= NN; IN++) { // 120

         // Do for each value of N in NVAL.

         N = NVAL( IN )
         LDA = MAX( 1, N )
         NIMAT = NTYPES
         if (N.LE.0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 110

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( N.GT.0 && .NOT.DOTYPE( IMAT ) ) GO TO 110

            // Set up parameters with DLATB4.

            dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST );

            ZEROT = IMAT.GE.8 && IMAT.LE.10
            if ( IMAT.LE.6 ) {

               // Type 1-6:  generate a symmetric tridiagonal matrix of
               // known condition number in lower triangular band storage.

               SRNAMT = 'DLATMS'
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'B', A, 2, WORK, INFO );

               // Check the error code from DLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 110
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

                  dlarnv(2, ISEED, N, D );
                  dlarnv(2, ISEED, N-1, E );

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

                  IX = IDAMAX( N, D, 1 )
                  DMAX = D( IX )
                  dscal(N, ANORM / DMAX, D, 1 );
                  if (N.GT.1) CALL DSCAL( N-1, ANORM / DMAX, E, 1 );

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
                     Z( 3 ) = E( IZERO )
                     E( IZERO-1 ) = ZERO
                     E( IZERO ) = ZERO
                  }
                  Z( 2 ) = D( IZERO )
                  D( IZERO ) = ZERO
               }
            }

            // Generate NRHS random solution vectors.

            IX = 1
            for (J = 1; J <= NRHS; J++) { // 40
               dlarnv(2, ISEED, N, XACT( IX ) );
               IX = IX + LDA
            } // 40

            // Set the right hand side.

            dlaptm(N, NRHS, ONE, D, E, XACT, LDA, ZERO, B, LDA );

            for (IFACT = 1; IFACT <= 2; IFACT++) { // 100
               if ( IFACT == 1 ) {
                  FACT = 'F'
               } else {
                  FACT = 'N'
               }

               // Compute the condition number for comparison with
               // the value returned by DPTSVX.

               if ( ZEROT ) {
                  if (IFACT == 1) GO TO 100;
                  RCONDC = ZERO

               } else if ( IFACT == 1 ) {

                  // Compute the 1-norm of A.

                  ANORM = DLANST( '1', N, D, E )

                  dcopy(N, D, 1, D( N+1 ), 1 );
                  if (N.GT.1) CALL DCOPY( N-1, E, 1, E( N+1 ), 1 );

                  // Factor the matrix A.

                  dpttrf(N, D( N+1 ), E( N+1 ), INFO );

                  // Use DPTTRS to solve for one column at a time of
                  // inv(A), computing the maximum column sum as we go.

                  AINVNM = ZERO
                  for (I = 1; I <= N; I++) { // 60
                     for (J = 1; J <= N; J++) { // 50
                        X( J ) = ZERO
                     } // 50
                     X( I ) = ONE
                     dpttrs(N, 1, D( N+1 ), E( N+1 ), X, LDA, INFO );
                     AINVNM = MAX( AINVNM, DASUM( N, X, 1 ) )
                  } // 60

                  // Compute the 1-norm condition number of A.

                  if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                     RCONDC = ONE
                  } else {
                     RCONDC = ( ONE / ANORM ) / AINVNM
                  }
               }

               if ( IFACT == 2 ) {

                  // --- Test DPTSV --

                  dcopy(N, D, 1, D( N+1 ), 1 );
                  if (N.GT.1) CALL DCOPY( N-1, E, 1, E( N+1 ), 1 );
                  dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                  // Factor A as L*D*L' and solve the system A*X = B.

                  SRNAMT = 'DPTSV '
                  dptsv(N, NRHS, D( N+1 ), E( N+1 ), X, LDA, INFO );

                  // Check error code from DPTSV .

                  if (INFO != IZERO) CALL ALAERH( PATH, 'DPTSV ', INFO, IZERO, ' ', N, N, 1, 1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                  NT = 0
                  if ( IZERO == 0 ) {

                     // Check the factorization by computing the ratio
                        // norm(L*D*L' - A) / (n * norm(A) * EPS )

                     dptt01(N, D, E, D( N+1 ), E( N+1 ), WORK, RESULT( 1 ) );

                     // Compute the residual in the solution.

                     dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     dptt02(N, NRHS, D, E, X, LDA, WORK, LDA, RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                     NT = 3
                  }

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= NT; K++) { // 70
                     if ( RESULT( K ).GE.THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )'DPTSV ', N, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1
                     }
                  } // 70
                  NRUN = NRUN + NT
               }

               // --- Test DPTSVX ---

               if ( IFACT.GT.1 ) {

                  // Initialize D( N+1:2*N ) and E( N+1:2*N ) to zero.

                  for (I = 1; I <= N - 1; I++) { // 80
                     D( N+I ) = ZERO
                     E( N+I ) = ZERO
                  } // 80
                  if (N.GT.0) D( N+N ) = ZERO;
               }

               dlaset('Full', N, NRHS, ZERO, ZERO, X, LDA );

               // Solve the system and compute the condition number and
               // error bounds using DPTSVX.

               SRNAMT = 'DPTSVX'
               dptsvx(FACT, N, NRHS, D, E, D( N+1 ), E( N+1 ), B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, INFO );

               // Check the error code from DPTSVX.

               if (INFO != IZERO) CALL ALAERH( PATH, 'DPTSVX', INFO, IZERO, FACT, N, N, 1, 1, NRHS, IMAT, NFAIL, NERRS, NOUT );
               if ( IZERO == 0 ) {
                  if ( IFACT == 2 ) {

                     // Check the factorization by computing the ratio
                        // norm(L*D*L' - A) / (n * norm(A) * EPS )

                     K1 = 1
                     dptt01(N, D, E, D( N+1 ), E( N+1 ), WORK, RESULT( 1 ) );
                  } else {
                     K1 = 2
                  }

                  // Compute the residual in the solution.

                  dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                  dptt02(N, NRHS, D, E, X, LDA, WORK, LDA, RESULT( 2 ) );

                  // Check solution from generated exact solution.

                  dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

                  // Check error bounds from iterative refinement.

                  dptt05(N, NRHS, D, E, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
               } else {
                  K1 = 6
               }

               // Check the reciprocal of the condition number.

               RESULT( 6 ) = DGET06( RCOND, RCONDC )

               // Print information about the tests that did not pass
               // the threshold.

               for (K = K1; K <= 6; K++) { // 90
                  if ( RESULT( K ).GE.THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9998 )'DPTSVX', FACT, N, IMAT, K, RESULT( K );
                     NFAIL = NFAIL + 1
                  }
               } // 90
               NRUN = NRUN + 7 - K1
            } // 100
         } // 110
      } // 120

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
      RETURN

      // End of DDRVPT

      }
