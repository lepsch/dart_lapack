      void zchktb(final Array<bool> DOTYPE_, final int NN, final Array<int> NVAL_, final int NNS, final Array<int> NSVAL_, final double THRESH, final bool TSTERR, final int NMAX, final int AB, final int AINV, final int B, final Array<double> X_, final Array<double> XACT_, final Array<double> WORK_, final Array<double> RWORK_, final Nout NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                NSVAL( * ), NVAL( * );
      double             RWORK( * );
      Complex         AB( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      int                NTYPE1, NTYPES;
      const              NTYPE1 = 9, NTYPES = 17 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      String             DIAG, NORM, TRANS, UPLO, XTYPE;
      String             PATH;
      int                I, IDIAG, IK, IMAT, IN, INFO, IRHS, ITRAN, IUPLO, J, K, KD, LDA, LDAB, N, NIMAT, NIMAT2, NK, NRHS;
      double             AINVNM, ANORM, RCOND, RCONDC, RCONDI, RCONDO, SCALE;
      String             TRANSS( NTRAN ), UPLOS( 2 );
      final                ISEED=Array<int>( 4 );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             ZLANTB, ZLANTR;
      // EXTERNAL lsame, ZLANTB, ZLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, ZCOPY, ZERRTR, ZGET04, ZLACPY, ZLARHS, ZLASET, ZLATBS, ZLATTB, ZTBCON, ZTBRFS, ZTBSV, ZTBT02, ZTBT03, ZTBT05, ZTBT06, ZTBTRS
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
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', TRANSS = 'N', 'T', 'C';

      // Initialize constants and the random number seed.

      final PATH = '${'Zomplex precision'[0]}TB';
      var NRUN = 0;
      var NFAIL = 0;
      var NERRS = Box(0);
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY[I - 1];
      } // 10

      // Test the error exits

      if (TSTERR) zerrtr( PATH, NOUT );
      infoc.INFOT = 0;

      for (IN = 1; IN <= NN; IN++) { // 140

         // Do for each value of N in NVAL

         final N = NVAL[IN];
         LDA = max( 1, N );
         XTYPE = 'N';
         NIMAT = NTYPE1;
         NIMAT2 = NTYPES;
         if ( N <= 0 ) {
            NIMAT = 1;
            NIMAT2 = NTYPE1 + 1;
         }

         NK = min( N+1, 4 );
         for (IK = 1; IK <= NK; IK++) { // 130

            // Do for KD = 0, N, (3N-1)/4, and (N+1)/4. This order makes
            // it easier to skip redundant values for small values of N.

            if ( IK == 1 ) {
               KD = 0;
            } else if ( IK == 2 ) {
               KD = max( N, 0 );
            } else if ( IK == 3 ) {
               KD = ( 3*N-1 ) / 4;
            } else if ( IK == 4 ) {
               KD = ( N+1 ) / 4;
            }
            LDAB = KD + 1;

            for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 90

               // Do the tests only if DOTYPE( IMAT ) is true.

               if( !DOTYPE[IMAT] ) GO TO 90;

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 80

                  // Do first for UPLO = 'U', then for UPLO = 'L'

                  final UPLO = UPLOS[IUPLO - 1];

                  // Call ZLATTB to generate a triangular test matrix.

                 srnamc.SRNAMT = 'ZLATTB';
                  zlattb(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, KD, AB, LDAB, X, WORK, RWORK, INFO );

                  // Set IDIAG = 1 for non-unit matrices, 2 for unit.

                  if ( lsame( DIAG, 'N' ) ) {
                     IDIAG = 1;
                  } else {
                     IDIAG = 2;
                  }

                  // Form the inverse of A so we can get a good estimate
                  // of RCONDC = 1/(norm(A) * norm(inv(A))).

                  zlaset('Full', N, N, Complex.zero, Complex.one, AINV, LDA );
                  if ( lsame( UPLO, 'U' ) ) {
                     for (J = 1; J <= N; J++) { // 20
                        ztbsv(UPLO, 'No transpose', DIAG, J, KD, AB, LDAB, AINV( ( J-1 )*LDA+1 ), 1 );
                     } // 20
                  } else {
                     for (J = 1; J <= N; J++) { // 30
                        ztbsv(UPLO, 'No transpose', DIAG, N-J+1, KD, AB( ( J-1 )*LDAB+1 ), LDAB, AINV( ( J-1 )*LDA+J ), 1 );
                     } // 30
                  }

                  // Compute the 1-norm condition number of A.

                  ANORM = ZLANTB( '1', UPLO, DIAG, N, KD, AB, LDAB, RWORK )                   AINVNM = ZLANTR( '1', UPLO, DIAG, N, N, AINV, LDA, RWORK );
                  if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                     RCONDO = ONE;
                  } else {
                     RCONDO = ( ONE / ANORM ) / AINVNM;
                  }

                  // Compute the infinity-norm condition number of A.

                  ANORM = ZLANTB( 'I', UPLO, DIAG, N, KD, AB, LDAB, RWORK )                   AINVNM = ZLANTR( 'I', UPLO, DIAG, N, N, AINV, LDA, RWORK );
                  if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                     RCONDI = ONE;
                  } else {
                     RCONDI = ( ONE / ANORM ) / AINVNM;
                  }

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 60
                     final NRHS = NSVAL[IRHS];
                     XTYPE = 'N';

                     for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 50

                     // Do for op(A) = A, A**T, or A**H.

                        final TRANS = TRANSS[ITRAN - 1];
                        if ( ITRAN == 1 ) {
                           NORM = 'O';
                           RCONDC = RCONDO;
                        } else {
                           NORM = 'I';
                           RCONDC = RCONDI;
                        }

// +    TEST 1
                     // Solve and compute residual for op(A)*x = b.

                       srnamc.SRNAMT = 'ZLARHS';
                        zlarhs(PATH, XTYPE, UPLO, TRANS, N, N, KD, IDIAG, NRHS, AB, LDAB, XACT, LDA, B, LDA, ISEED, INFO );
                        XTYPE = 'C';
                        zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                       srnamc.SRNAMT = 'ZTBTRS';
                        ztbtrs(UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, X, LDA, INFO );

                     // Check error code from ZTBTRS.

                        if (INFO != 0) alaerh( PATH, 'ZTBTRS', INFO, 0, UPLO + TRANS // DIAG, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        ztbt02(UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, X, LDA, B, LDA, WORK, RWORK, RESULT( 1 ) );

// +    TEST 2
                     // Check solution from generated exact solution.

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 2 ) );

// +    TESTS 3, 4, and 5
                     // Use iterative refinement to improve the solution
                     // and compute error bounds.

                       srnamc.SRNAMT = 'ZTBRFS';
                        ztbrfs(UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                     // Check error code from ZTBRFS.

                        if (INFO != 0) alaerh( PATH, 'ZTBRFS', INFO, 0, UPLO + TRANS // DIAG, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                        ztbt05(UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 1; K <= 5; K++) { // 40
                           if ( RESULT[K] >= THRESH ) {
                              if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                              NOUT.println( 9999 )UPLO, TRANS, DIAG, N, KD, NRHS, IMAT, K, RESULT[K];
                              NFAIL++;
                           }
                        } // 40
                        NRUN +=  5;
                     } // 50
                  } // 60

// +    TEST 6
                     // Get an estimate of RCOND = 1/CNDNUM.

                  for (ITRAN = 1; ITRAN <= 2; ITRAN++) { // 70
                     if ( ITRAN == 1 ) {
                        NORM = 'O';
                        RCONDC = RCONDO;
                     } else {
                        NORM = 'I';
                        RCONDC = RCONDI;
                     }
                    srnamc.SRNAMT = 'ZTBCON';
                     ztbcon(NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK, RWORK, INFO );

                     // Check error code from ZTBCON.

                     if (INFO != 0) alaerh( PATH, 'ZTBCON', INFO, 0, NORM + UPLO // DIAG, N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT );

                     ztbt06(RCOND, RCONDC, UPLO, DIAG, N, KD, AB, LDAB, RWORK, RESULT( 6 ) );

                     // Print the test ratio if it is >= THRESH.

                     if ( RESULT[6] >= THRESH ) {
                        if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                        NOUT.println( 9998 ) 'ZTBCON', NORM, UPLO, DIAG, N, KD, IMAT, 6, RESULT( 6 );
                        NFAIL++;
                     }
                     NRUN++;
                  } // 70
               } // 80
            } // 90

            // Use pathological test matrices to test ZLATBS.

            for (IMAT = NTYPE1 + 1; IMAT <= NIMAT2; IMAT++) { // 120

               // Do the tests only if DOTYPE( IMAT ) is true.

               if( !DOTYPE[IMAT] ) GO TO 120;

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 110

                  // Do first for UPLO = 'U', then for UPLO = 'L'

                  final UPLO = UPLOS[IUPLO - 1];
                  for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 100

                     // Do for op(A) = A, A**T, and A**H.

                     final TRANS = TRANSS[ITRAN - 1];

                     // Call ZLATTB to generate a triangular test matrix.

                    srnamc.SRNAMT = 'ZLATTB';
                     zlattb(IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, AB, LDAB, X, WORK, RWORK, INFO );

// +    TEST 7
                     // Solve the system op(A)*x = b

                    srnamc.SRNAMT = 'ZLATBS';
                     zcopy(N, X, 1, B, 1 );
                     zlatbs(UPLO, TRANS, DIAG, 'N', N, KD, AB, LDAB, B, SCALE, RWORK, INFO );

                     // Check error code from ZLATBS.

                     if (INFO != 0) alaerh( PATH, 'ZLATBS', INFO, 0, UPLO + TRANS // DIAG // 'N', N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT );

                     ztbt03(UPLO, TRANS, DIAG, N, KD, 1, AB, LDAB, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 7 ) );

// +    TEST 8
                     // Solve op(A)*x = b again with NORMIN = 'Y'.

                     zcopy(N, X, 1, B, 1 );
                     zlatbs(UPLO, TRANS, DIAG, 'Y', N, KD, AB, LDAB, B, SCALE, RWORK, INFO );

                     // Check error code from ZLATBS.

                     if (INFO != 0) alaerh( PATH, 'ZLATBS', INFO, 0, UPLO + TRANS // DIAG // 'Y', N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT );

                     ztbt03(UPLO, TRANS, DIAG, N, KD, 1, AB, LDAB, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) );

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( RESULT[7] >= THRESH ) {
                        if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                        NOUT.println( 9997 )'ZLATBS', UPLO, TRANS, DIAG, 'N', N, KD, IMAT, 7, RESULT( 7 );
                        NFAIL++;
                     }
                     if ( RESULT[8] >= THRESH ) {
                        if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                        NOUT.println( 9997 )'ZLATBS', UPLO, TRANS, DIAG, 'Y', N, KD, IMAT, 8, RESULT( 8 );
                        NFAIL++;
                     }
                     NRUN +=  2;
                  } // 100
               } // 110
            } // 120
         } // 130
      } // 140

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO=\'${.a1}\', TRANS=\'${TRANS.a1}\', DIAG=\'${.a1}\', N=${N.i5}, KD=${.i5}, NRHS=${.i5}, type ${IMAT.i2}, test(${.i2})=${.g12_5};
 9998 FORMAT(' ${}( \'${.a1}\'${.a1}\'${.a1}\',${.i5},${.i5},  ... ), type ${IMAT.i2}, test(${.i2})=${.g12_5};
 9997 FORMAT(' ${}( \'${.a1}\'${.a1}\'${.a1}\'${.a1}\',${.i5},${.i5}, ...  ),  type ${IMAT.i2}, test(${.i1})=${RESULT[].g12_5};
      }
