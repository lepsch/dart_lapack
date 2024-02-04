      void zdrvab(DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX, A, AFAC, B, X, WORK, RWORK, SWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NM, NMAX, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                MVAL( * ), NSVAL( * ), IWORK( * );
      double             RWORK( * );
      Complex            SWORK( * );
      Complex         A( * ), AFAC( * ), B( * ), WORK( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IM, IMAT, INFO, IOFF, IRHS, IZERO, KL, KU, LDA, M, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      double             ANORM, CNDNUM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
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
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 2006, 2007, 2008, 2009 ];
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      KASE = 0;
      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'GE';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      INFOT = 0;

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

            // Set up parameters with ZLATB4 and generate a test matrix
            // with ZLATMS.

            zlatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

            SRNAMT = 'ZLATMS';
            zlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

            // Check error code from ZLATMS.

            if ( INFO != 0 ) {
               alaerh(PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
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
                  zlaset('Full', M, N-IZERO+1, Complex.zero, Complex.zero, A( IOFF+1 ), LDA );
               }
            } else {
               IZERO = 0;
            }

            for (IRHS = 1; IRHS <= NNS; IRHS++) { // 60
               NRHS = NSVAL( IRHS );
               XTYPE = 'N';
               TRANS = 'N';

               SRNAMT = 'ZLARHS';
               zlarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, X, LDA, B, LDA, ISEED, INFO );

               SRNAMT = 'ZCGESV';

               KASE = KASE + 1;

               zlacpy('Full', M, N, A, LDA, AFAC, LDA );

               zcgesv(N, NRHS, A, LDA, IWORK, B, LDA, X, LDA, WORK, SWORK, RWORK, ITER, INFO);

               if (ITER < 0) {
                   zlacpy('Full', M, N, AFAC, LDA, A, LDA );
               }

               // Check error code from ZCGESV. This should be the same as
               // the one of DGETRF.

               if ( INFO != IZERO ) {

                  if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                  NERRS = NERRS + 1;

                  if ( INFO != IZERO && IZERO != 0 ) {
                     WRITE( NOUT, FMT = 9988 )'ZCGESV',INFO, IZERO,M,IMAT;
                  } else {
                     WRITE( NOUT, FMT = 9975 )'ZCGESV',INFO, M, IMAT;
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
                 // NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1

               // If double precision has been used (ITER<0), we want
                 // NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
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
         WRITE( NOUT, FMT = 9996 )'ZCGESV', NFAIL, NRUN;
      } else {
         WRITE( NOUT, FMT = 9995 )'ZCGESV', NRUN;
      }
      if ( NERRS > 0 ) {
         WRITE( NOUT, FMT = 9994 )NERRS;
      }

 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') =', G12.5 );
 9996 FORMAT( 1X, A6, ': ', I6, ' out of ', I6, ' tests failed to pass the threshold' );
 9995 FORMAT( /1X, 'All tests for ', A6, ' routines passed the threshold ( ', I6, ' tests run)' );
 9994 FORMAT( 6X, I6, ' error messages recorded' );

      // SUBNAM, INFO, INFOE, M, IMAT

 9988 FORMAT( ' *** ', A6, ' returned with INFO =', I5, ' instead of ', I5, / ' ==> M =', I5, ', type ', I2 );

      // SUBNAM, INFO, M, IMAT

 9975 FORMAT( ' *** Error code from ', A6, '=', I5, ' for M=', I5, ', type ', I2 );
 8999 FORMAT( / 1X, A3, ':  General dense matrices' );
 8979 FORMAT( 4X, '1. Diagonal', 24X, '7. Last n/2 columns zero', / 4X, '2. Upper triangular', 16X, '8. Random, CNDNUM = sqrt(0.1/EPS)', / 4X, '3. Lower triangular', 16X, '9. Random, CNDNUM = 0.1/EPS', / 4X, '4. Random, CNDNUM = 2', 13X, '10. Scaled near underflow', / 4X, '5. First column zero', 14X, '11. Scaled near overflow', / 4X, '6. Last column zero' );
 8960 FORMAT( 3X, I2, ': norm_1( B - A * X )  / ', '( norm_1(A) * norm_1(X) * EPS * sqrt(N) ) > 1 if ITERREF', / 4x, 'or norm_1( B - A * X )  / ', '( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF' );
      return;
      }
