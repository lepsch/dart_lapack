      SUBROUTINE DDRVAB( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX, A, AFAC, B, X, WORK, RWORK, SWORK, IWORK, NOUT );

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
      REAL               SWORK(*);
      double             A( * ), AFAC( * ), B( * ), RWORK( * ), WORK( * ), X( * );
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
      // EXTERNAL ALAERH, ALAHD, DGET08, DLACPY, DLARHS, DLASET, DLATB4, DLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT
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
      DATA               ISEEDY / 2006, 2007, 2008, 2009 /;
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      KASE = 0;
      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'GE';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10

      INFOT = 0;

      // Do for each value of M in MVAL

      for (IM = 1; IM <= NM; IM++) { // 120
         M = MVAL( IM );
         LDA = MAX( 1, M );

         N = M;
         NIMAT = NTYPES;
         if (M <= 0 || N <= 0) NIMAT = 1;

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( !DOTYPE( IMAT ) ) GO TO 100;

            // Skip types 5, 6, or 7 if the matrix size is too small.

            ZEROT = IMAT >= 5 && IMAT <= 7;
            if (ZEROT && N < IMAT-4) GO TO 100;

            // Set up parameters with DLATB4 and generate a test matrix
            // with DLATMS.

            dlatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

            SRNAMT = 'DLATMS';
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
                  IZERO = MIN( M, N );
               } else {
                  IZERO = MIN( M, N ) / 2 + 1;
               }
               IOFF = ( IZERO-1 )*LDA;
               if ( IMAT < 7 ) {
                  for (I = 1; I <= M; I++) { // 20
                     A( IOFF+I ) = ZERO;
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

               SRNAMT = 'DLARHS';
               dlarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, X, LDA, B, LDA, ISEED, INFO );

               SRNAMT = 'DSGESV';

               KASE = KASE + 1;

               dlacpy('Full', M, N, A, LDA, AFAC, LDA );

               dsgesv(N, NRHS, A, LDA, IWORK, B, LDA, X, LDA, WORK, SWORK, ITER, INFO);

               if (ITER < 0) {
                   dlacpy('Full', M, N, AFAC, LDA, A, LDA );
               }

               // Check error code from DSGESV. This should be the same as
               // the one of DGETRF.

               if ( INFO != IZERO ) {

                  if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH );
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
                 // NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1

               // If double precision has been used (ITER<0), we want
                 // NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
               // (Cf. the linear solver testing routines)

               if ((THRESH <= 0.0e+00) || ((ITER >= 0) && (N > 0) && (RESULT(1) >= SQRT(DBLE(N)))) || ((ITER < 0) && (RESULT(1) >= THRESH))) {

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
 8960 FORMAT( 3X, I2, ': norm_1( B - A * X )  / ', '( norm_1(A) * norm_1(X) * EPS * SQRT(N) ) > 1 if ITERREF', / 4x, 'or norm_1( B - A * X )  / ', '( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF' );
      return;

      // End of DDRVAB

      }
