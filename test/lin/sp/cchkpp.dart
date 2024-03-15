      void cchkpp(final int DOTYPE, final int NN, final int NVAL, final int NNS, final int NSVAL, final int THRESH, final int TSTERR, final int NMAX, final int A, final int AFAC, final int AINV, final int B, final int X, final int XACT, final Array<double> _WORK_, final Array<double> RWORK_, final int NOUT,) {
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      double               THRESH;
      bool               DOTYPE( * );
      int                NSVAL( * ), NVAL( * );
      double               RWORK( * );
      Complex            A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      bool               ZEROT;
      String             DIST, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IMAT, IN, INFO, IOFF, IRHS, IUPLO, IZERO, K, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NPP, NRHS, NRUN;
      double               ANORM, CNDNUM, RCOND, RCONDC;
      String             PACKS( 2 ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- REAL               CLANHP, SGET06;
      // EXTERNAL CLANHP, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CCOPY, CERRPO, CGET04, CLACPY, CLAIPD, CLARHS, CLATB4, CLATMS, CPPCON, CPPRFS, CPPT01, CPPT02, CPPT03, CPPT05, CPPTRF, CPPTRI, CPPTRS
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
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', PACKS = 'C', 'R';

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Complex precision';
      PATH[2: 3] = 'PP';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) cerrpo( PATH, NOUT );
      INFOT = 0;

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 110
         N = NVAL( IN );
         LDA = max( N, 1 );
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 100;

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 5;
            if (ZEROT && N < IMAT-2) GO TO 100;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 90
               UPLO = UPLOS( IUPLO );
               PACKIT = PACKS( IUPLO );

               // Set up parameters with CLATB4 and generate a test matrix
               // with CLATMS.

               clatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

              srnamc.SRNAMT = 'CLATMS';
               clatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO );

               // Check error code from CLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'CLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 90;
               }

               // For types 3-5, zero one row and column of the matrix to
               // test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT == 3 ) {
                     IZERO = 1;
                  } else if ( IMAT == 4 ) {
                     IZERO = N;
                  } else {
                     IZERO = N / 2 + 1;
                  }

                  // Set row and column IZERO of A to 0.

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
                  IZERO = 0;
               }

               // Set the imaginary part of the diagonals.

               if ( IUPLO == 1 ) {
                  claipd(N, A, 2, 1 );
               } else {
                  claipd(N, A, N, -1 );
               }

               // Compute the L*L' or U'*U factorization of the matrix.

               NPP = N*( N+1 ) / 2;
               ccopy(NPP, A, 1, AFAC, 1 );
              srnamc.SRNAMT = 'CPPTRF';
               cpptrf(UPLO, N, AFAC, INFO );

               // Check error code from CPPTRF.

               if ( INFO != IZERO ) {
                  alaerh(PATH, 'CPPTRF', INFO, IZERO, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 90;
               }

               // Skip the tests if INFO is not 0.

               if (INFO != 0) GO TO 90;

// +    TEST 1
               // Reconstruct matrix from factors and compute residual.

               ccopy(NPP, AFAC, 1, AINV, 1 );
               cppt01(UPLO, N, A, AINV, RWORK, RESULT( 1 ) );

// +    TEST 2
               // Form the inverse and compute the residual.

               ccopy(NPP, AFAC, 1, AINV, 1 );
              srnamc.SRNAMT = 'CPPTRI';
               cpptri(UPLO, N, AINV, INFO );

               // Check error code from CPPTRI.

               if (INFO != 0) alaerh( PATH, 'CPPTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

               cppt03(UPLO, N, A, AINV, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );

               // Print information about the tests that did not pass
               // the threshold.

               for (K = 1; K <= 2; K++) { // 60
                  if ( RESULT( K ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, K, RESULT( K );
                     NFAIL = NFAIL + 1;
                  }
               } // 60
               NRUN = NRUN + 2;

               for (IRHS = 1; IRHS <= NNS; IRHS++) { // 80
                  NRHS = NSVAL( IRHS );

// +    TEST 3
               // Solve and compute residual for  A * X = B.

                 srnamc.SRNAMT = 'CLARHS';
                  clarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  clacpy('Full', N, NRHS, B, LDA, X, LDA );

                 srnamc.SRNAMT = 'CPPTRS';
                  cpptrs(UPLO, N, NRHS, AFAC, X, LDA, INFO );

               // Check error code from CPPTRS.

                  if (INFO != 0) alaerh( PATH, 'CPPTRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  clacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                  cppt02(UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

// +    TEST 4
               // Check solution from generated exact solution.

                  cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );

// +    TESTS 5, 6, and 7
               // Use iterative refinement to improve the solution.

                 srnamc.SRNAMT = 'CPPRFS';
                  cpprfs(UPLO, N, NRHS, A, AFAC, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

               // Check error code from CPPRFS.

                  if (INFO != 0) alaerh( PATH, 'CPPRFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                  cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) );
                  cppt05(UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 6 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 3; K <= 7; K++) { // 70
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 70
                  NRUN = NRUN + 5;
               } // 80

// +    TEST 8
               // Get an estimate of RCOND = 1/CNDNUM.

               ANORM = CLANHP( '1', UPLO, N, A, RWORK );
              srnamc.SRNAMT = 'CPPCON';
               cppcon(UPLO, N, AFAC, ANORM, RCOND, WORK, RWORK, INFO );

               // Check error code from CPPCON.

               if (INFO != 0) alaerh( PATH, 'CPPCON', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

               RESULT[8] = SGET06( RCOND, RCONDC );

               // Print the test ratio if greater than or equal to THRESH.

               if ( RESULT( 8 ) >= THRESH ) {
                  if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                  WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, 8, RESULT( 8 );
                  NFAIL = NFAIL + 1;
               }
               NRUN = NRUN + 1;

            } // 90
         } // 100
      } // 110

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = \'${.a1}\', N =${.i5}, type ${.i2}, test ${.i2}, ratio =${.g12_5};
 9998 FORMAT( ' UPLO = \'${.a1}\', N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) =${.g12_5};
      }