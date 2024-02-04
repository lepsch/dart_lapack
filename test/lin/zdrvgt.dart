      void zdrvgt(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, AF, B, X, XACT, WORK, RWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
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
      const              NTESTS = 6 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, FACT, TRANS, TYPE;
      String             PATH;
      int                I, IFACT, IMAT, IN, INFO, ITRAN, IX, IZERO, J, K, K1, KL, KOFF, KU, LDA, M, MODE, N, NERRS, NFAIL, NIMAT, NRUN, NT;
      double             AINVNM, ANORM, ANORMI, ANORMO, COND, RCOND, RCONDC, RCONDI, RCONDO;
      // ..
      // .. Local Arrays ..
      String             TRANSS( 3 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS ), Z( 3 );
      // ..
      // .. External Functions ..
      //- double             DGET06, DZASUM, ZLANGT;
      // EXTERNAL DGET06, DZASUM, ZLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, ZCOPY, ZDSCAL, ZERRVX, ZGET04, ZGTSV, ZGTSVX, ZGTT01, ZGTT02, ZGTT05, ZGTTRF, ZGTTRS, ZLACPY, ZLAGTM, ZLARNV, ZLASET, ZLATB4, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX
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

      if (TSTERR) zerrvx( PATH, NOUT );
      INFOT = 0;

      for (IN = 1; IN <= NN; IN++) { // 140

         // Do for each value of N in NVAL.

         N = NVAL( IN );
         M = max( N-1, 0 );
         LDA = max( 1, N );
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 130

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 130;

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
                  GO TO 130;
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

                  // Generate a matrix with elements from [-1,1].

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
                  Z[2] = (A( N )).toDouble();
                  A[N] = ZERO;
                  if ( N > 1 ) {
                     Z[3] = (A( 1 )).toDouble();
                     A[1] = ZERO;
                  }
               } else if ( IMAT == 9 ) {
                  IZERO = N;
                  Z[1] = (A( 3*N-2 )).toDouble();
                  Z[2] = (A( 2*N-1 )).toDouble();
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

            for (IFACT = 1; IFACT <= 2; IFACT++) { // 120
               if ( IFACT == 1 ) {
                  FACT = 'F';
               } else {
                  FACT = 'N';
               }

               // Compute the condition number for comparison with
               // the value returned by ZGTSVX.

               if ( ZEROT ) {
                  if (IFACT == 1) GO TO 120;
                  RCONDO = ZERO;
                  RCONDI = ZERO;

               } else if ( IFACT == 1 ) {
                  zcopy(N+2*M, A, 1, AF, 1 );

                  // Compute the 1-norm and infinity-norm of A.

                  ANORMO = ZLANGT( '1', N, A, A( M+1 ), A( N+M+1 ) );
                  ANORMI = ZLANGT( 'I', N, A, A( M+1 ), A( N+M+1 ) );

                  // Factor the matrix A.

                  zgttrf(N, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, INFO );

                  // Use ZGTTRS to solve for one column at a time of
                  // inv(A), computing the maximum column sum as we go.

                  AINVNM = ZERO;
                  for (I = 1; I <= N; I++) { // 40
                     for (J = 1; J <= N; J++) { // 30
                        X[J] = ZERO;
                     } // 30
                     X[I] = ONE;
                     zgttrs('No transpose', N, 1, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO );
                     AINVNM = max( AINVNM, DZASUM( N, X, 1 ) );
                  } // 40

                  // Compute the 1-norm condition number of A.

                  if ( ANORMO <= ZERO || AINVNM <= ZERO ) {
                     RCONDO = ONE;
                  } else {
                     RCONDO = ( ONE / ANORMO ) / AINVNM;
                  }

                  // Use ZGTTRS to solve for one column at a time of
                  // inv(A'), computing the maximum column sum as we go.

                  AINVNM = ZERO;
                  for (I = 1; I <= N; I++) { // 60
                     for (J = 1; J <= N; J++) { // 50
                        X[J] = ZERO;
                     } // 50
                     X[I] = ONE;
                     zgttrs('Conjugate transpose', N, 1, AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, X, LDA, INFO );
                     AINVNM = max( AINVNM, DZASUM( N, X, 1 ) );
                  } // 60

                  // Compute the infinity-norm condition number of A.

                  if ( ANORMI <= ZERO || AINVNM <= ZERO ) {
                     RCONDI = ONE;
                  } else {
                     RCONDI = ( ONE / ANORMI ) / AINVNM;
                  }
               }

               for (ITRAN = 1; ITRAN <= 3; ITRAN++) { // 110
                  TRANS = TRANSS( ITRAN );
                  if ( ITRAN == 1 ) {
                     RCONDC = RCONDO;
                  } else {
                     RCONDC = RCONDI;
                  }

                  // Generate NRHS random solution vectors.

                  IX = 1;
                  for (J = 1; J <= NRHS; J++) { // 70
                     zlarnv(2, ISEED, N, XACT( IX ) );
                     IX = IX + LDA;
                  } // 70

                  // Set the right hand side.

                  zlagtm(TRANS, N, NRHS, ONE, A, A( M+1 ), A( N+M+1 ), XACT, LDA, ZERO, B, LDA );

                  if ( IFACT == 2 && ITRAN == 1 ) {

                     // --- Test ZGTSV  ---

                     // Solve the system using Gaussian elimination with
                     // partial pivoting.

                     zcopy(N+2*M, A, 1, AF, 1 );
                     zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                     SRNAMT = 'ZGTSV ';
                     zgtsv(N, NRHS, AF, AF( M+1 ), AF( N+M+1 ), X, LDA, INFO );

                     // Check error code from ZGTSV .

                     if (INFO != IZERO) alaerh( PATH, 'ZGTSV ', INFO, IZERO, ' ', N, N, 1, 1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                     NT = 1;
                     if ( IZERO == 0 ) {

                        // Check residual of computed solution.

                        zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                        zgtt02(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), X, LDA, WORK, LDA, RESULT( 2 ) );

                        // Check solution from generated exact solution.

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        NT = 3;
                     }

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 2; K <= NT; K++) { // 80
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                           WRITE( NOUT, FMT = 9999 )'ZGTSV ', N, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1;
                        }
                     } // 80
                     NRUN = NRUN + NT - 1;
                  }

                  // --- Test ZGTSVX ---

                  if ( IFACT > 1 ) {

                     // Initialize AF to zero.

                     for (I = 1; I <= 3*N - 2; I++) { // 90
                        AF[I] = ZERO;
                     } // 90
                  }
                  zlaset('Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDA );

                  // Solve the system and compute the condition number and
                  // error bounds using ZGTSVX.

                  SRNAMT = 'ZGTSVX';
                  zgtsvx(FACT, TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                  // Check the error code from ZGTSVX.

                  if (INFO != IZERO) alaerh( PATH, 'ZGTSVX', INFO, IZERO, FACT // TRANS, N, N, 1, 1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  if ( IFACT >= 2 ) {

                     // Reconstruct matrix from factors and compute
                     // residual.

                     zgtt01(N, A, A( M+1 ), A( N+M+1 ), AF, AF( M+1 ), AF( N+M+1 ), AF( N+2*M+1 ), IWORK, WORK, LDA, RWORK, RESULT( 1 ) );
                     K1 = 1;
                  } else {
                     K1 = 2;
                  }

                  if ( INFO == 0 ) {
                     TRFCON = false;

                     // Check residual of computed solution.

                     zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     zgtt02(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), X, LDA, WORK, LDA, RESULT( 2 ) );

                     // Check solution from generated exact solution.

                     zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

                     // Check the error bounds from iterative refinement.

                     zgtt05(TRANS, N, NRHS, A, A( M+1 ), A( N+M+1 ), B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                     NT = 5;
                  }

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = K1; K <= NT; K++) { // 100
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9998 )'ZGTSVX', FACT, TRANS, N, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 100

                  // Check the reciprocal of the condition number.

                  RESULT[6] = DGET06( RCOND, RCONDC );
                  if ( RESULT( 6 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) aladhd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9998 )'ZGTSVX', FACT, TRANS, N, IMAT, K, RESULT( K );
                     NFAIL = NFAIL + 1;
                  }
                  NRUN = NRUN + NT - K1 + 2;

               } // 110
            } // 120
         } // 130
      } // 140

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', N =', I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 );
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', TRANS=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 );
      return;
      }
