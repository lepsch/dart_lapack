      SUBROUTINE DCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A, D, E, B, X, XACT, WORK, RWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NN, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                NSVAL( * ), NVAL( * );
      double             A( * ), B( * ), D( * ), E( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
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
      // EXTERNAL ALAERH, ALAHD, ALASUM, DCOPY, DERRGT, DGET04, DLACPY, DLAPTM, DLARNV, DLATB4, DLATMS, DPTCON, DPTRFS, DPTT01, DPTT02, DPTT05, DPTTRF, DPTTRS, DSCAL
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
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
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
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL DERRGT( PATH, NOUT )
      INFOT = 0

      for (IN = 1; IN <= NN; IN++) { // 110

         // Do for each value of N in NVAL.

         N = NVAL( IN )
         LDA = MAX( 1, N )
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( N.GT.0 .AND. .NOT.DOTYPE( IMAT ) ) GO TO 100

            // Set up parameters with DLATB4.

            dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, COND, DIST );

            ZEROT = IMAT.GE.8 .AND. IMAT.LE.10
            if ( IMAT.LE.6 ) {

               // Type 1-6:  generate a symmetric tridiagonal matrix of
               // known condition number in lower triangular band storage.

               SRNAMT = 'DLATMS'
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, COND, ANORM, KL, KU, 'B', A, 2, WORK, INFO );

               // Check the error code from DLATMS.

               if ( INFO.NE.0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, ' ', N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100
               }
               IZERO = 0

               // Copy the matrix to D and E.

               IA = 1
               DO 20 I = 1, N - 1
                  D( I ) = A( IA )
                  E( I ) = A( IA+1 )
                  IA = IA + 2
   20          CONTINUE
               IF( N.GT.0 ) D( N ) = A( IA )
            } else {

               // Type 7-12:  generate a diagonally dominant matrix with
               // unknown condition number in the vectors D and E.

               if ( .NOT.ZEROT .OR. .NOT.DOTYPE( 7 ) ) {

                  // Let D and E have values from [-1,1].

                  dlarnv(2, ISEED, N, D );
                  dlarnv(2, ISEED, N-1, E );

                  // Make the tridiagonal matrix diagonally dominant.

                  if ( N.EQ.1 ) {
                     D( 1 ) = ABS( D( 1 ) )
                  } else {
                     D( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                     D( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
                     DO 30 I = 2, N - 1
                        D( I ) = ABS( D( I ) ) + ABS( E( I ) ) + ABS( E( I-1 ) )
   30                CONTINUE
                  }

                  // Scale D and E so the maximum element is ANORM.

                  IX = IDAMAX( N, D, 1 )
                  DMAX = D( IX )
                  dscal(N, ANORM / DMAX, D, 1 );
                  dscal(N-1, ANORM / DMAX, E, 1 );

               } else if ( IZERO.GT.0 ) {

                  // Reuse the last matrix by copying back the zeroed out
                  // elements.

                  if ( IZERO.EQ.1 ) {
                     D( 1 ) = Z( 2 )
                     IF( N.GT.1 ) E( 1 ) = Z( 3 )
                  } else if ( IZERO.EQ.N ) {
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
               if ( IMAT.EQ.8 ) {
                  IZERO = 1
                  Z( 2 ) = D( 1 )
                  D( 1 ) = ZERO
                  if ( N.GT.1 ) {
                     Z( 3 ) = E( 1 )
                     E( 1 ) = ZERO
                  }
               } else if ( IMAT.EQ.9 ) {
                  IZERO = N
                  if ( N.GT.1 ) {
                     Z( 1 ) = E( N-1 )
                     E( N-1 ) = ZERO
                  }
                  Z( 2 ) = D( N )
                  D( N ) = ZERO
               } else if ( IMAT.EQ.10 ) {
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

            dcopy(N, D, 1, D( N+1 ), 1 );
            IF( N.GT.1 ) CALL DCOPY( N-1, E, 1, E( N+1 ), 1 )

*+    TEST 1
            // Factor A as L*D*L' and compute the ratio
               // norm(L*D*L' - A) / (n * norm(A) * EPS )

            dpttrf(N, D( N+1 ), E( N+1 ), INFO );

            // Check error code from DPTTRF.

            if ( INFO.NE.IZERO ) {
               alaerh(PATH, 'DPTTRF', INFO, IZERO, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               GO TO 100
            }

            if ( INFO.GT.0 ) {
               RCONDC = ZERO
               GO TO 90
            }

            dptt01(N, D, E, D( N+1 ), E( N+1 ), WORK, RESULT( 1 ) );

            // Print the test ratio if greater than or equal to THRESH.

            if ( RESULT( 1 ).GE.THRESH ) {
               IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 )
               NFAIL = NFAIL + 1
            }
            NRUN = NRUN + 1

            // Compute RCONDC = 1 / (norm(A) * norm(inv(A))

            // Compute norm(A).

            ANORM = DLANST( '1', N, D, E )

            // Use DPTTRS to solve for one column at a time of inv(A),
            // computing the maximum column sum as we go.

            AINVNM = ZERO
            for (I = 1; I <= N; I++) { // 50
               for (J = 1; J <= N; J++) { // 40
                  X( J ) = ZERO
   40          CONTINUE
               X( I ) = ONE
               dpttrs(N, 1, D( N+1 ), E( N+1 ), X, LDA, INFO );
               AINVNM = MAX( AINVNM, DASUM( N, X, 1 ) )
   50       CONTINUE
            RCONDC = ONE / MAX( ONE, ANORM*AINVNM )

            for (IRHS = 1; IRHS <= NNS; IRHS++) { // 80
               NRHS = NSVAL( IRHS )

            // Generate NRHS random solution vectors.

               IX = 1
               for (J = 1; J <= NRHS; J++) { // 60
                  dlarnv(2, ISEED, N, XACT( IX ) );
                  IX = IX + LDA
   60          CONTINUE

            // Set the right hand side.

               dlaptm(N, NRHS, ONE, D, E, XACT, LDA, ZERO, B, LDA );

*+    TEST 2
            // Solve A*x = b and compute the residual.

               dlacpy('Full', N, NRHS, B, LDA, X, LDA );
               dpttrs(N, NRHS, D( N+1 ), E( N+1 ), X, LDA, INFO );

            // Check error code from DPTTRS.

               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DPTTRS', INFO, 0, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

               dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
               dptt02(N, NRHS, D, E, X, LDA, WORK, LDA, RESULT( 2 ) );

*+    TEST 3
            // Check solution from generated exact solution.

               dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

*+    TESTS 4, 5, and 6
            // Use iterative refinement to improve the solution.

               SRNAMT = 'DPTRFS'
               dptrfs(N, NRHS, D, E, D( N+1 ), E( N+1 ), B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, INFO );

            // Check error code from DPTRFS.

               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DPTRFS', INFO, 0, ' ', N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

               dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) )                CALL DPTT05( N, NRHS, D, E, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

            // Print information about the tests that did not pass the
            // threshold.

               for (K = 2; K <= 6; K++) { // 70
                  if ( RESULT( K ).GE.THRESH ) {
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9998 )N, NRHS, IMAT, K, RESULT( K )
                     NFAIL = NFAIL + 1
                  }
   70          CONTINUE
               NRUN = NRUN + 5
   80       CONTINUE

*+    TEST 7
            // Estimate the reciprocal of the condition number of the
            // matrix.

   90       CONTINUE
            SRNAMT = 'DPTCON'
            dptcon(N, D( N+1 ), E( N+1 ), ANORM, RCOND, RWORK, INFO );

            // Check error code from DPTCON.

            IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DPTCON', INFO, 0, ' ', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

            RESULT( 7 ) = DGET06( RCOND, RCONDC )

            // Print the test ratio if greater than or equal to THRESH.

            if ( RESULT( 7 ).GE.THRESH ) {
               IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )N, IMAT, 7, RESULT( 7 )
               NFAIL = NFAIL + 1
            }
            NRUN = NRUN + 1
  100    CONTINUE
  110 CONTINUE

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' N =', I5, ', type ', I2, ', test ', I2, ', ratio = ', G12.5 )
 9998 FORMAT( ' N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') = ', G12.5 )
      RETURN

      // End of DCHKPT

      }
