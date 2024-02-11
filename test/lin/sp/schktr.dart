      void schktr(final int DOTYPE, final int NN, final int NVAL, final int NNB, final int NBVAL, final int NNS, final int NSVAL, final int THRESH, final int TSTERR, final int NMAX, final int A, final int AINV, final int B, final int X, final int XACT, final Array<double> _WORK, final Array<double> RWORK, final Array<int> IWORK, final int NOUT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double               THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double               A( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

      int                NTYPE1, NTYPES;
      const              NTYPE1 = 10, NTYPES = 18 ;
      int                NTESTS;
      const              NTESTS = 10 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      String             DIAG, NORM, TRANS, UPLO, XTYPE;
      String             PATH;
      int                I, IDIAG, IMAT, IN, INB, INFO, IRHS, ITRAN, IUPLO, K, LDA, N, NB, NERRS, NFAIL, NRHS, NRUN;
      double               AINVNM, ANORM, BIGNUM, DUMMY, RCOND, RCONDC, RCONDI, RCONDO, RES, SCALE, SLAMCH;
      String             TRANSS( NTRAN ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS ), SCALE3( 2 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLANTR;
      // EXTERNAL lsame, SLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, SCOPY, SERRTR, SGET04, SLACPY, SLARHS, SLATRS, SLATRS3, SLATTR, SSCAL, STRCON, STRRFS, STRT01, STRT02, STRT03, STRT05, STRT06, STRTRI, STRTRS, XLAENV, SLAMCH
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, IOUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', TRANSS = 'N', 'T', 'C';

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Single precision';
      PATH[2: 3] = 'TR';
      BIGNUM = SLAMCH('Overflow') / SLAMCH('Precision');
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) serrtr( PATH, NOUT );
      INFOT = 0;
      xlaenv(2, 2 );

      for (IN = 1; IN <= NN; IN++) { // 120

         // Do for each value of N in NVAL

         N = NVAL( IN );
         LDA = max( 1, N );
         XTYPE = 'N';

         for (IMAT = 1; IMAT <= NTYPE1; IMAT++) { // 80

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 80;

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 70

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO );

               // Call SLATTR to generate a triangular test matrix.

              srnamc.SRNAMT = 'SLATTR';
               slattr(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, A, LDA, X, WORK, INFO );

               // Set IDIAG = 1 for non-unit matrices, 2 for unit.

               if ( lsame( DIAG, 'N' ) ) {
                  IDIAG = 1;
               } else {
                  IDIAG = 2;
               }

               for (INB = 1; INB <= NNB; INB++) { // 60

                  // Do for each blocksize in NBVAL

                  NB = NBVAL( INB );
                  xlaenv(1, NB );

// +    TEST 1
                  // Form the inverse of A.

                  slacpy(UPLO, N, N, A, LDA, AINV, LDA );
                 srnamc.SRNAMT = 'STRTRI';
                  strtri(UPLO, DIAG, N, AINV, LDA, INFO );

                  // Check error code from STRTRI.

                  if (INFO != 0) alaerh( PATH, 'STRTRI', INFO, 0, UPLO + DIAG, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                  // Compute the infinity-norm condition number of A.

                  ANORM = SLANTR( 'I', UPLO, DIAG, N, N, A, LDA, RWORK );
                  AINVNM = SLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, RWORK );
                  if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                     RCONDI = ONE;
                  } else {
                     RCONDI = ( ONE / ANORM ) / AINVNM;
                  }

                  // Compute the residual for the triangular matrix times
                  // its inverse.  Also compute the 1-norm condition number
                  // of A.

                  strt01(UPLO, DIAG, N, A, LDA, AINV, LDA, RCONDO, RWORK, RESULT( 1 ) );

                  // Print the test ratio if it is >= THRESH.

                  if ( RESULT( 1 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9999 )UPLO, DIAG, N, NB, IMAT, 1, RESULT( 1 );
                     NFAIL = NFAIL + 1;
                  }
                  NRUN = NRUN + 1;

                  // Skip remaining tests if not the first block size.

                  if (INB != 1) GO TO 60;

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 40
                     NRHS = NSVAL( IRHS );
                     XTYPE = 'N';

                     for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 30

                     // Do for op(A) = A, A**T, or A**H.

                        TRANS = TRANSS( ITRAN );
                        if ( ITRAN == 1 ) {
                           NORM = 'O';
                           RCONDC = RCONDO;
                        } else {
                           NORM = 'I';
                           RCONDC = RCONDI;
                        }

// +    TEST 2
                        // Solve and compute residual for op(A)*x = b.

                       srnamc.SRNAMT = 'SLARHS';
                        slarhs(PATH, XTYPE, UPLO, TRANS, N, N, 0, IDIAG, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                        XTYPE = 'C';
                        slacpy('Full', N, NRHS, B, LDA, X, LDA );

                       srnamc.SRNAMT = 'STRTRS';
                        strtrs(UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDA, INFO );

                        // Check error code from STRTRS.

                        if (INFO != 0) alaerh( PATH, 'STRTRS', INFO, 0, UPLO + TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        // This line is needed on a Sun SPARCstation.

                        if (N > 0) DUMMY = A( 1 );

                        strt02(UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDA, B, LDA, WORK, RESULT( 2 ) );

// +    TEST 3
                        // Check solution from generated exact solution.

                        sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, and 6
                        // Use iterative refinement to improve the solution
                        // and compute error bounds.

                       srnamc.SRNAMT = 'STRRFS';
                        strrfs(UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO );

                        // Check error code from STRRFS.

                        if (INFO != 0) alaerh( PATH, 'STRRFS', INFO, 0, UPLO + TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                        strt05(UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 2; K <= 6; K++) { // 20
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9998 )UPLO, TRANS, DIAG, N, NRHS, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 20
                        NRUN = NRUN + 5;
                     } // 30
                  } // 40

// +    TEST 7
                        // Get an estimate of RCOND = 1/CNDNUM.

                  for (ITRAN = 1; ITRAN <= 2; ITRAN++) { // 50
                     if ( ITRAN == 1 ) {
                        NORM = 'O';
                        RCONDC = RCONDO;
                     } else {
                        NORM = 'I';
                        RCONDC = RCONDI;
                     }
                    srnamc.SRNAMT = 'STRCON';
                     strcon(NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, IWORK, INFO );

                        // Check error code from STRCON.

                     if (INFO != 0) alaerh( PATH, 'STRCON', INFO, 0, NORM + UPLO // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     strt06(RCOND, RCONDC, UPLO, DIAG, N, A, LDA, RWORK, RESULT( 7 ) );

                     // Print the test ratio if it is >= THRESH.

                     if ( RESULT( 7 ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9997 )NORM, UPLO, N, IMAT, 7, RESULT( 7 );
                        NFAIL = NFAIL + 1;
                     }
                     NRUN = NRUN + 1;
                  } // 50
               } // 60
            } // 70
         } // 80

         // Use pathological test matrices to test SLATRS.

         for (IMAT = NTYPE1 + 1; IMAT <= NTYPES; IMAT++) { // 110

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 110;

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 100

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO );
               for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 90

                  // Do for op(A) = A, A**T, and A**H.

                  TRANS = TRANSS( ITRAN );

                  // Call SLATTR to generate a triangular test matrix.

                 srnamc.SRNAMT = 'SLATTR';
                  slattr(IMAT, UPLO, TRANS, DIAG, ISEED, N, A, LDA, X, WORK, INFO );

// +    TEST 8
                  // Solve the system op(A)*x = b.

                 srnamc.SRNAMT = 'SLATRS';
                  scopy(N, X, 1, B, 1 );
                  slatrs(UPLO, TRANS, DIAG, 'N', N, A, LDA, B, SCALE, RWORK, INFO );

                  // Check error code from SLATRS.

                  if (INFO != 0) alaerh( PATH, 'SLATRS', INFO, 0, UPLO + TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  strt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) );

// +    TEST 9
                  // Solve op(A)*X = b again with NORMIN = 'Y'.

                  scopy(N, X, 1, B( N+1 ), 1 );
                  slatrs(UPLO, TRANS, DIAG, 'Y', N, A, LDA, B( N+1 ), SCALE, RWORK, INFO );

                  // Check error code from SLATRS.

                  if (INFO != 0) alaerh( PATH, 'SLATRS', INFO, 0, UPLO + TRANS // DIAG // 'Y', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  strt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RESULT( 9 ) );

// +    TEST 10
                  // Solve op(A)*X = B

                 srnamc.SRNAMT = 'SLATRS3';
                  scopy(N, X, 1, B, 1 );
                  scopy(N, X, 1, B( N+1 ), 1 );
                  sscal(N, BIGNUM, B( N+1 ), 1 );
                  slatrs3(UPLO, TRANS, DIAG, 'N', N, 2, A, LDA, B, max(1, N), SCALE3, RWORK, WORK, NMAX, INFO );

                  // Check error code from SLATRS3.

                  if (INFO != 0) alaerh( PATH, 'SLATRS3', INFO, 0, UPLO + TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  strt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE3( 1 ), RWORK, ONE, B( 1 ), LDA, X, LDA, WORK, RESULT( 10 ) );
                  sscal(N, BIGNUM, X, 1 );
                  strt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE3( 2 ), RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RES );
                  RESULT[10] = max( RESULT( 10 ), RES );

                  // Print information about the tests that did not pass
                  // the threshold.

                  if ( RESULT( 8 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'SLATRS', UPLO, TRANS, DIAG, 'N', N, IMAT, 8, RESULT( 8 );
                     NFAIL = NFAIL + 1;
                  }
                  if ( RESULT( 9 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'SLATRS', UPLO, TRANS, DIAG, 'Y', N, IMAT, 9, RESULT( 9 );
                     NFAIL = NFAIL + 1;
                  }
                  if ( RESULT( 10 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'SLATRS3', UPLO, TRANS, DIAG, 'N', N, IMAT, 10, RESULT( 10 );
                     NFAIL = NFAIL + 1;
                  }
                  NRUN = NRUN + 3;
               } // 90
            } // 100
         } // 110
      } // 120

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO=''${.a1}'', DIAG=''${.a1}'', N=${.i5}, NB=${.i4}, type ${.i2}, test(${.i2})= ${.g12_5};
 9998 FORMAT( ' UPLO=''${.a1}'', TRANS=''${.a1}'', DIAG=''${.a1}'', N=${.i5}, NB=${.i4}, type ${.i2}, test(${.i2})= ${.g12_5};
 9997 FORMAT( ' NORM=''${.a1}'', UPLO =''${.a1}'', N=${.i5},${' ' * 11} type ${.i2}, test(${.i2})=${.g12_5};
 9996 FORMAT(' ${}( ''${.a1}''${.a1}''${.a1}''${.a1}'',${.i5}, ... ), type ${.i2}, test(${.i2})=${.g12_5};
      }
