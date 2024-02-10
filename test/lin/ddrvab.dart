import 'common.dart';

      void ddrvab(DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX, A, AFAC, B, X, final Array<double> _WORK, final Array<double> RWORK, SWORK, final Array<int> IWORK, final int NOUT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NM, NMAX, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                MVAL( * ), NSVAL( * ), IWORK( * );
      double               SWORK(*);
      double             A( * ), AFAC( * ), B( * ), RWORK( * ), WORK( * ), X( * );
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
      int                I, IM, IMAT, INFO, IOFF, IRHS, IZERO, KL, KU, LDA, M, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      double             ANORM, CNDNUM;
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. Local Variables ..
      int                ITER, KASE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, DGET08, DLACPY, DLARHS, DLASET, DLATB4, DLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 2006, 2007, 2008, 2009 ];

      // Initialize constants and the random number seed.

      KASE = 0;
      PATH = '${'Double precision'[0]}';
      PATH[2: 3] = 'GE';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      infoc.INFOT = 0;

      // Do for each value of M in MVAL

      for (IM = 1; IM <= NM; IM++) { // 120
         M = MVAL( IM );
         LDA = max( 1, M );

         N = M;
         NIMAT = NTYPES;
         if (M <= 0 || N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 100;

            // Skip types 5, 6, or 7 if the matrix size is too small.

            ZEROT = IMAT >= 5 && IMAT <= 7;
            if (ZEROT && N < IMAT-4) GO TO 100;

            // Set up parameters with DLATB4 and generate a test matrix
            // with DLATMS.

            dlatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

            srnamc.SRNAMT = 'DLATMS';
            dlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

            // Check error code from DLATMS.

            if ( INFO != 0 ) {
               alaerh(PATH, 'DLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
               GO TO 100;
            }

            // For types 5-7, zero one or more columns of the matrix to
            // test that INFO is returned correctly.

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
                  dlaset('Full', M, N-IZERO+1, ZERO, ZERO, A( IOFF+1 ), LDA );
               }
            } else {
               IZERO = 0;
            }

            for (IRHS = 1; IRHS <= NNS; IRHS++) { // 60
               NRHS = NSVAL( IRHS );
               XTYPE = 'N';
               TRANS = 'N';

               srnamc.SRNAMT = 'DLARHS';
               dlarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, X, LDA, B, LDA, ISEED, INFO );

               srnamc.SRNAMT = 'DSGESV';

               KASE = KASE + 1;

               dlacpy('Full', M, N, A, LDA, AFAC, LDA );

               dsgesv(N, NRHS, A, LDA, IWORK, B, LDA, X, LDA, WORK, SWORK, ITER, INFO);

               if (ITER < 0) {
                   dlacpy('Full', M, N, AFAC, LDA, A, LDA );
               }

               // Check error code from DSGESV. This should be the same as
               // the one of DGETRF.

               if ( INFO != IZERO ) {

                  if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                  NERRS = NERRS + 1;

                  if ( INFO != IZERO && IZERO != 0 ) {
                     WRITE( NOUT, FMT = 9988 )'DSGESV',INFO, IZERO,M,IMAT;
                  } else {
                     WRITE( NOUT, FMT = 9975 )'DSGESV',INFO, M, IMAT;
                  }
               }

               // Skip the remaining test if the matrix is singular.

               if (INFO != 0) GO TO 100;

               // Check the quality of the solution

               dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );

               dget08(TRANS, N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 1 ) );

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
                     WRITE( NOUT, FMT = 8999 )'DGE';
                     WRITE( NOUT, FMT = '( '' Matrix types:'' )' );
                     WRITE( NOUT, FMT = 8979 );
                     WRITE( NOUT, FMT = '( '' Test ratios:'' )' );
                     WRITE( NOUT, FMT = 8960 )1;
                     WRITE( NOUT, FMT = '( '' Messages:'' )' );
                  }

                  WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, 1, RESULT( 1 );
                  NFAIL = NFAIL + 1;
               }
               NRUN = NRUN + 1;
            } // 60
         } // 100
      } // 120

      // Print a summary of the results.

      if ( NFAIL > 0 ) {
         WRITE( NOUT, FMT = 9996 )'DSGESV', NFAIL, NRUN;
      } else {
         WRITE( NOUT, FMT = 9995 )'DSGESV', NRUN;
      }
      if ( NERRS > 0 ) {
         WRITE( NOUT, FMT = 9994 )NERRS;
      }

 9998 FORMAT( ' TRANS=''${.a1}'', N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) =${.g12_5};
 9996 FORMAT(' ${.a6}: ${.i6} out of ${.i6} tests failed to pass the threshold' );
 9995 FORMAT('\n All tests for ${.a6} routines passed the threshold ( ${.i6} tests run)' );
 9994 FORMAT('${' ' * 6}${.i6} error messages recorded' );

      // SUBNAM, INFO, INFOE, M, IMAT

 9988 FORMAT( ' *** ${.a6} returned with INFO =${.i5} instead of ', I5, / ' ==> M =${.i5}, type ${.i2}');

      // SUBNAM, INFO, M, IMAT

 9975 FORMAT( ' *** Error code from ${.a6}=${.i5} for M=${.i5}, type ${.i2}');
 8999 FORMAT('\n ${.a3}:  General dense matrices' );
 8979 FORMAT('    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n${' ' * 4}2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n${' ' * 4}3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n${' ' * 4}4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n${' ' * 4}5. First column zero${' ' * 14}11. Scaled near overflow\n${' ' * 4}6. Last column zero' );
 8960 FORMAT('${' ' * 3}${.i2}: norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS * sqrt(N) ) > 1 if ITERREF\n${' ' * 4}or norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF' );
      }
