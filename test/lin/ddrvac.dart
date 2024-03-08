import 'common.dart';

      void ddrvac(final int DOTYPE, final int NM, final int MVAL, final int NNS, final int NSVAL, final int THRESH, final int NMAX, final int A, final int AFAC, final int B, final int X, final Array<double> WORK_, final Array<double> RWORK_, final int SWORK, final int NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NMAX, NM, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                MVAL( * ), NSVAL( * );
      double               SWORK(*);
      double             A( * ), AFAC( * ), B( * ), RWORK( * ), WORK( * ), X( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      bool               ZEROT;
      String             DIST, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IM, IMAT, INFO, IOFF, IRHS, IUPLO, IZERO, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      double             ANORM, CNDNUM;
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. Local Variables ..
      int                ITER, KASE;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, DLACPY, DLARHS, DLASET, DLATB4, DLATMS, DPOT06, DSPOSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, SQRT
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
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];

      // Initialize constants and the random number seed.

      KASE = 0;
      PATH = '${'Double precision'[0]}';
      PATH[2: 3] = 'PO';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      infoc.INFOT = 0;

      // Do for each value of N in MVAL

      for (IM = 1; IM <= NM; IM++) { // 120
         N = MVAL( IM );
         LDA = max( N, 1 );
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 110

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 110;

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 5;
            if (ZEROT && N < IMAT-2) GO TO 110;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 100
               UPLO = UPLOS( IUPLO );

               // Set up parameters with DLATB4 and generate a test matrix
               // with DLATMS.

               dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               srnamc.SRNAMT = 'DLATMS';
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from DLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100;
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
                  IOFF = ( IZERO-1 )*LDA;

                  // Set row and column IZERO of A to 0.

                  if ( IUPLO == 1 ) {
                     for (I = 1; I <= IZERO - 1; I++) { // 20
                        A[IOFF+I] = ZERO;
                     } // 20
                     IOFF = IOFF + IZERO;
                     for (I = IZERO; I <= N; I++) { // 30
                        A[IOFF] = ZERO;
                        IOFF = IOFF + LDA;
                     } // 30
                  } else {
                     IOFF = IZERO;
                     for (I = 1; I <= IZERO - 1; I++) { // 40
                        A[IOFF] = ZERO;
                        IOFF = IOFF + LDA;
                     } // 40
                     IOFF = IOFF - IZERO;
                     for (I = IZERO; I <= N; I++) { // 50
                        A[IOFF+I] = ZERO;
                     } // 50
                  }
               } else {
                  IZERO = 0;
               }

               for (IRHS = 1; IRHS <= NNS; IRHS++) { // 60
                  NRHS = NSVAL( IRHS );
                  XTYPE = 'N';

                  // Form an exact solution and set the right hand side.

                  srnamc.SRNAMT = 'DLARHS';
                  dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, X, LDA, B, LDA, ISEED, INFO );

                  // Compute the L*L' or U'*U factorization of the
                  // matrix and solve the system.

                  srnamc.SRNAMT = 'DSPOSV ';
                  KASE = KASE + 1;

                  dlacpy('All', N, N, A, LDA, AFAC, LDA);

                  dsposv(UPLO, N, NRHS, AFAC, LDA, B, LDA, X, LDA, WORK, SWORK, ITER, INFO );

                  if (ITER < 0) {
                     dlacpy('All', N, N, A, LDA, AFAC, LDA );
                  }

                  // Check error code from DSPOSV .

                  if ( INFO != IZERO ) {

                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     NERRS = NERRS + 1;

                     if ( INFO != IZERO && IZERO != 0 ) {
                        WRITE( NOUT, FMT = 9988 )'DSPOSV',INFO,IZERO,N, IMAT;
                     } else {
                        WRITE( NOUT, FMT = 9975 )'DSPOSV',INFO,N,IMAT;
                     }
                  }

                  // Skip the remaining test if the matrix is singular.

                  if (INFO != 0) GO TO 110;

                  // Check the quality of the solution

                  dlacpy('All', N, NRHS, B, LDA, WORK, LDA );

                  dpot06(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 1 ) );

                  // Check if the test passes the testing.
                  // Print information about the tests that did not
                  // pass the testing.

                  // If iterative refinement has been used and claimed to
                  // be successful (ITER>0), we want
                  // NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS*SRQT(N)) < 1

                  // If double precision has been used (ITER<0), we want
                  // NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS) < THRES
                  // (Cf. the linear solver testing routines)

                  if ((THRESH <= 0.0e+00) || ((ITER >= 0) && (N > 0) && (RESULT(1) >= sqrt(N.toDouble()))) || ((ITER < 0) && (RESULT(1) >= THRESH))) {

                     if ( NFAIL == 0 && NERRS == 0 ) {
                        WRITE( NOUT, FMT = 8999 )'DPO';
                        WRITE( NOUT, FMT = ' Matrix types:' );
                        WRITE( NOUT, FMT = 8979 );
                        WRITE( NOUT, FMT = ' Test ratios:' );
                        WRITE( NOUT, FMT = 8960 )1;
                        WRITE( NOUT, FMT = ' Messages:' );
                     }

                     WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, 1, RESULT( 1 );

                     NFAIL = NFAIL + 1;

                  }

                  NRUN = NRUN + 1;

               } // 60
            } // 100
         } // 110
      } // 120

      // Print a summary of the results.

      if ( NFAIL > 0 ) {
         WRITE( NOUT, FMT = 9996 )'DSPOSV', NFAIL, NRUN;
      } else {
         WRITE( NOUT, FMT = 9995 )'DSPOSV', NRUN;
      }
      if ( NERRS > 0 ) {
         WRITE( NOUT, FMT = 9994 )NERRS;
      }

 9998 FORMAT( ' UPLO=\'${.a1}\', N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) =${.g12_5};
 9996 FORMAT(' ${.a6}: ${.i6} out of ${.i6} tests failed to pass the threshold' );
 9995 FORMAT('\n All tests for ${.a6} routines passed the threshold ( ${.i6} tests run)' );
 9994 FORMAT('${' ' * 6}${.i6} error messages recorded' );

      // SUBNAM, INFO, INFOE, N, IMAT

 9988 FORMAT( ' *** ${.a6} returned with INFO =${.i5} instead of ', I5, / ' ==> N =${.i5}, type ${.i2}');

      // SUBNAM, INFO, N, IMAT

 9975 FORMAT( ' *** Error code from ${.a6}=${.i5} for M=${.i5}, type ${.i2}');
 8999 FORMAT('\n ${.a3}:  positive definite dense matrices' );
 8979 FORMAT('    1. Diagonal${' ' * 24}7. Last n/2 columns zero\n    2. Upper triangular${' ' * 16}8. Random, CNDNUM = sqrt(0.1/EPS)\n    3. Lower triangular${' ' * 16}9. Random, CNDNUM = 0.1/EPS\n    4. Random, CNDNUM = 2${' ' * 13}10. Scaled near underflow\n    5. First column zero${' ' * 14}11. Scaled near overflow\n    6. Last column zero' );
 8960 FORMAT('   ${.i2}: norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS * sqrt(N) ) > 1 if ITERREF\n    or norm_1( B - A * X )  / ( norm_1(A) * norm_1(X) * EPS ) > THRES if DPOTRF' );

      }
