      void zdrvab(final Array<bool> DOTYPE_, final int NM, final Array<int> MVAL_, final int NNS, final Array<int> NSVAL_, final int THRESH, final int NMAX, final int A, final int AFAC, final int B, final int X, final Array<double> WORK_, final Array<double> RWORK_, final int SWORK, final Array<int> IWORK_, final Nout NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NM, NMAX, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                MVAL( * ), NSVAL( * ), IWORK( * );
      double             RWORK( * );
      Complex            SWORK( * );
      Complex         A( * ), AFAC( * ), B( * ), WORK( * ), X( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      bool               ZEROT;
      String             DIST, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IM, IMAT, INFO, IOFF, IRHS, IZERO, KL, KU, LDA, M, MODE, N, NIMAT, NRHS;
      double             ANORM, CNDNUM;
      final                ISEED=Array<int>( 4 );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. Local Variables ..
      int                ITER, KASE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ZGET08, ZLACPY, ZLARHS, ZLASET, ZLATB4, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, DBLE, MAX, MIN, SQRT
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 2006, 2007, 2008, 2009 ];

      // Initialize constants and the random number seed.

      KASE = 0;
      final PATH = '${'Zomplex precision'[0]}GE';
      var NRUN = 0;
      var NFAIL = 0;
      var NERRS = Box(0);
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY[I - 1];
      } // 10

      infoc.INFOT = 0;

      // Do for each value of M in MVAL

      for (IM = 1; IM <= NM; IM++) { // 120
         final M = MVAL[IM];
         final LDA = max( 1, M );

         N = M;
            final NIMAT = M <= 0 || N <= 0 ? 1 : NTYPES;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE[IMAT] ) GO TO 100;

            // Skip types 5, 6, or 7 if the matrix size is too small.

            final ZEROT = IMAT >= 5 && IMAT <= 7;
            if (ZEROT && N < IMAT-4) GO TO 100;

            // Set up parameters with ZLATB4 and generate a test matrix
            // with ZLATMS.

            final (:TYPE,:KL,:KU,:ANORM,:MODE,:CNDNUM,:DIST) = zlatb4(PATH, IMAT, M, N);

           srnamc.SRNAMT = 'ZLATMS';
            zlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

            // Check error code from ZLATMS.

            if ( INFO.value != 0 ) {
               alaerh(PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               GO TO 100;
            }

            // For types 5-7, zero one or more columns of the matrix to
            // test that INFO is returned correctly.

            final int IZERO;
            if ( ZEROT ) {
               if ( IMAT == 5 ) {
                  IZERO = 1;
               } else if ( IMAT == 6 ) {
                  IZERO = min( M, N );
               } else {
                  IZERO = min( M, N ) / 2 + 1;
               }
               IOFF = ( IZERO-1 )*LDA;
               if ( IMAT < 7 ) {
                  for (I = 1; I <= M; I++) { // 20
                     A[IOFF+I] = ZERO;
                  } // 20
               } else {
                  zlaset('Full', M, N-IZERO+1, Complex.zero, Complex.zero, A( IOFF+1 ), LDA );
               }
            } else {
               IZERO = 0;
            }

            for (IRHS = 1; IRHS <= NNS; IRHS++) { // 60
               final NRHS = NSVAL[IRHS];
               XTYPE = 'N';
               TRANS = 'N';

              srnamc.SRNAMT = 'ZLARHS';
               zlarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, X, LDA, B, LDA, ISEED, INFO );

              srnamc.SRNAMT = 'ZCGESV';

               KASE = KASE + 1;

               zlacpy('Full', M, N, A, LDA, AFAC, LDA );

               zcgesv(N, NRHS, A, LDA, IWORK, B, LDA, X, LDA, WORK, SWORK, RWORK, ITER, INFO);

               if (ITER < 0) {
                   zlacpy('Full', M, N, AFAC, LDA, A, LDA );
               }

               // Check error code from ZCGESV. This should be the same as
               // the one of DGETRF.

               if ( INFO.value != IZERO ) {

                  if (NFAIL == 0 && NERRS.value == 0) alahd( NOUT, PATH );
                  NERRS = NERRS + 1;

                  if ( INFO.value != IZERO && IZERO != 0 ) {
                     NOUT.println( 9988 )'ZCGESV',INFO, IZERO,M,IMAT;
                  } else {
                     NOUT.println( 9975 )'ZCGESV',INFO, M, IMAT;
                  }
               }

               // Skip the remaining test if the matrix is singular.

               if (INFO != 0) GO TO 100;

               // Check the quality of the solution

               zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );

               zget08(TRANS, N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 1 ) );

               // Check if the test passes the testing.
               // Print information about the tests that did not
               // pass the testing.

               // If iterative refinement has been used and claimed to
               // be successful (ITER>0), we want
               //   NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1

               // If double precision has been used (ITER<0), we want
               //   NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
               // (Cf. the linear solver testing routines)

               if ((THRESH <= 0.0e+00) || ((ITER >= 0) && (N > 0) && (RESULT(1) >= sqrt(N.toDouble()))) || ((ITER < 0) && (RESULT(1) >= THRESH))) {

                  if ( NFAIL == 0 && NERRS == 0 ) {
                     NOUT.println( 8999 )'DGE';
                     NOUT.println( ' Matrix types:' );
                     NOUT.println( 8979 );
                     NOUT.println( ' Test ratios:' );
                     NOUT.println( 8960 )1;
                     NOUT.println( ' Messages:' );
                  }

                  NOUT.println( 9998 )TRANS, N, NRHS, IMAT, 1, RESULT( 1 );
                  NFAIL++;
               }
               NRUN++;
            } // 60
         } // 100
      } // 120

      // Print a summary of the results.

      if ( NFAIL > 0 ) {
         NOUT.println( 9996 )'ZCGESV', NFAIL, NRUN;
      } else {
         NOUT.println( 9995 )'ZCGESV', NRUN;
      }
      if ( NERRS > 0 ) {
         NOUT.println( 9994 )NERRS;
      }

 9998 FORMAT( ' TRANS=\'${TRANS.a1}\', N =${N.i5}, NRHS=${NRHS.i3}, type ${IMAT.i2}, test(${.i2}) =${RESULT[].g12_5};
 9996 FORMAT(' ${.a6}: ${.i6} out of ${.i6} tests failed to pass the threshold' );
 9995 FORMAT('\n All tests for ${.a6} routines passed the threshold ( ${.i6} tests run)' );
 9994 FORMAT('${' ' * 6}${.i6} error messages recorded' );

      // SUBNAM, INFO, INFOE, M, IMAT

 9988 FORMAT( ' *** ${.a6} returned with INFO =${.i5} instead of ${.i5}\n ==> M =${M.i5}, type ${IMAT.i2}');

      // SUBNAM, INFO, M, IMAT

 9975 FORMAT( ' *** Error code from ${.a6}=${.i5} for M=${M.i5}, type ${IMAT.i2}');
 8999 FORMAT('\n ${.a3}:  General dense matrices' );
 8979 FORMAT('    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n    2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n    3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n    4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n    5. First column zero${' ' * 14}11. Scaled near overflow\n    6. Last column zero' );
 8960 FORMAT('   ${.i2}: norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS * sqrt(N) ) > 1 if ITERREF\n    or norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF' );
      }
