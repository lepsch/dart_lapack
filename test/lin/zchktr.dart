      void zchktr(DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AINV, B, X, XACT, final Array<double> _WORK, final Array<double> RWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                NBVAL( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      Complex         A( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      int                NTYPE1, NTYPES;
      const              NTYPE1 = 10, NTYPES = 18 ;
      int                NTESTS;
      const              NTESTS = 10 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      String             DIAG, NORM, TRANS, UPLO, XTYPE;
      String             PATH;
      int                I, IDIAG, IMAT, IN, INB, INFO, IRHS, ITRAN, IUPLO, K, LDA, N, NB, NERRS, NFAIL, NRHS, NRUN;
      double             AINVNM, ANORM, BIGNUM, DUMMY, RCOND, RCONDC, RCONDI, RCONDO, RES, SCALE, DLAMCH;
      String             TRANSS( NTRAN ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS ), RWORK2( 2*NMAX ), SCALE3( 2 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             ZLANTR;
      // EXTERNAL lsame, ZLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DLAMCH, XLAENV, ZCOPY, ZDSCAL, ZERRTR, ZGET04, ZLACPY, ZLARHS, ZLATRS, ZLATRS3, ZLATTR, ZTRCON, ZTRRFS, ZTRT01, ZTRT02, ZTRT03, ZTRT05, ZTRT06, ZTRTRI, ZTRTRS
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, IOUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', TRANSS = 'N', 'T', 'C';

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'TR';
      BIGNUM = dlamch('Overflow') / dlamch('Precision');
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) zerrtr( PATH, NOUT );
      infoc.INFOT = 0;

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

               // Call ZLATTR to generate a triangular test matrix.

              srnamc.SRNAMT = 'ZLATTR';
               zlattr(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, A, LDA, X, WORK, RWORK, INFO );

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

                  zlacpy(UPLO, N, N, A, LDA, AINV, LDA );
                 srnamc.SRNAMT = 'ZTRTRI';
                  ztrtri(UPLO, DIAG, N, AINV, LDA, INFO );

                  // Check error code from ZTRTRI.

                  if (INFO != 0) alaerh( PATH, 'ZTRTRI', INFO, 0, UPLO + DIAG, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                  // Compute the infinity-norm condition number of A.

                  ANORM = ZLANTR( 'I', UPLO, DIAG, N, N, A, LDA, RWORK );
                  AINVNM = ZLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, RWORK );
                  if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                     RCONDI = ONE;
                  } else {
                     RCONDI = ( ONE / ANORM ) / AINVNM;
                  }

                  // Compute the residual for the triangular matrix times
                  // its inverse.  Also compute the 1-norm condition number
                  // of A.

                  ztrt01(UPLO, DIAG, N, A, LDA, AINV, LDA, RCONDO, RWORK, RESULT( 1 ) );
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

                       srnamc.SRNAMT = 'ZLARHS';
                        zlarhs(PATH, XTYPE, UPLO, TRANS, N, N, 0, IDIAG, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                        XTYPE = 'C';
                        zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                       srnamc.SRNAMT = 'ZTRTRS';
                        ztrtrs(UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDA, INFO );

                        // Check error code from ZTRTRS.

                        if (INFO != 0) alaerh( PATH, 'ZTRTRS', INFO, 0, UPLO + TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        // This line is needed on a Sun SPARCstation.

                        if (N > 0) DUMMY = (A( 1 )).toDouble();

                        ztrt02(UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDA, B, LDA, WORK, RWORK, RESULT( 2 ) );

// +    TEST 3
                        // Check solution from generated exact solution.

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, and 6
                        // Use iterative refinement to improve the solution
                        // and compute error bounds.

                       srnamc.SRNAMT = 'ZTRRFS';
                        ztrrfs(UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                        // Check error code from ZTRRFS.

                        if (INFO != 0) alaerh( PATH, 'ZTRRFS', INFO, 0, UPLO + TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                        ztrt05(UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

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
                    srnamc.SRNAMT = 'ZTRCON';
                     ztrcon(NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK, RWORK, INFO );

                        // Check error code from ZTRCON.

                     if (INFO != 0) alaerh( PATH, 'ZTRCON', INFO, 0, NORM + UPLO // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     ztrt06(RCOND, RCONDC, UPLO, DIAG, N, A, LDA, RWORK, RESULT( 7 ) );

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

         // Use pathological test matrices to test ZLATRS.

         for (IMAT = NTYPE1 + 1; IMAT <= NTYPES; IMAT++) { // 110

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 110;

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 100

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO );
               for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 90

                  // Do for op(A) = A, A**T, and A**H.

                  TRANS = TRANSS( ITRAN );

                  // Call ZLATTR to generate a triangular test matrix.

                 srnamc.SRNAMT = 'ZLATTR';
                  zlattr(IMAT, UPLO, TRANS, DIAG, ISEED, N, A, LDA, X, WORK, RWORK, INFO );

// +    TEST 8
                  // Solve the system op(A)*x = b.

                 srnamc.SRNAMT = 'ZLATRS';
                  zcopy(N, X, 1, B, 1 );
                  zlatrs(UPLO, TRANS, DIAG, 'N', N, A, LDA, B, SCALE, RWORK, INFO );

                  // Check error code from ZLATRS.

                  if (INFO != 0) alaerh( PATH, 'ZLATRS', INFO, 0, UPLO + TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  ztrt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) );

// +    TEST 9
                  // Solve op(A)*X = b again with NORMIN = 'Y'.

                  zcopy(N, X, 1, B( N+1 ), 1 );
                  zlatrs(UPLO, TRANS, DIAG, 'Y', N, A, LDA, B( N+1 ), SCALE, RWORK, INFO );

                  // Check error code from ZLATRS.

                  if (INFO != 0) alaerh( PATH, 'ZLATRS', INFO, 0, UPLO + TRANS // DIAG // 'Y', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  ztrt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE, RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RESULT( 9 ) );

// +    TEST 10
                  // Solve op(A)*X = B

                 srnamc.SRNAMT = 'ZLATRS3';
                  zcopy(N, X, 1, B, 1 );
                  zcopy(N, X, 1, B( N+1 ), 1 );
                  zdscal(N, BIGNUM, B( N+1 ), 1 );
                  zlatrs3(UPLO, TRANS, DIAG, 'N', N, 2, A, LDA, B, max(1, N), SCALE3, RWORK, RWORK2, 2*NMAX, INFO );

                  // Check error code from ZLATRS3.

                  if (INFO != 0) alaerh( PATH, 'ZLATRS3', INFO, 0, UPLO + TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  ztrt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE3( 1 ), RWORK, ONE, B( 1 ), LDA, X, LDA, WORK, RESULT( 10 ) );
                  zdscal(N, BIGNUM, X, 1 );
                  ztrt03(UPLO, TRANS, DIAG, N, 1, A, LDA, SCALE3( 2 ), RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RES );
                  RESULT[10] = max( RESULT( 10 ), RES );

                  // Print information about the tests that did not pass
                  // the threshold.

                  if ( RESULT( 8 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'ZLATRS', UPLO, TRANS, DIAG, 'N', N, IMAT, 8, RESULT( 8 );
                     NFAIL = NFAIL + 1;
                  }
                  if ( RESULT( 9 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'ZLATRS', UPLO, TRANS, DIAG, 'Y', N, IMAT, 9, RESULT( 9 );
                     NFAIL = NFAIL + 1;
                  }
                  if ( RESULT( 10 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'ZLATRS3', UPLO, TRANS, DIAG, 'N', N, IMAT, 10, RESULT( 10 );
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
