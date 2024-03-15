      void cerrpo(final int PATH, final int NUNIT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      String             C2;
      int                I, INFO, J;
      double               ANRM, RCOND;
      double               R( NMAX ), R1( NMAX ), R2( NMAX );
      Complex            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CPBCON, CPBEQU, CPBRFS, CPBTF2, CPBTRF, CPBTRS, CPOCON, CPOEQU, CPORFS, CPOTF2, CPOTRF, CPOTRI, CPOTRS, CPPCON, CPPEQU, CPPRFS, CPPTRF, CPPTRI, CPPTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF[I][J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
         } // 10
         B[J] = 0.;
         R1[J] = 0.;
         R2[J] = 0.;
         W[J] = 0.;
         X[J] = 0.;
      } // 20
      ANRM = 1.;
      OK = true;

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite matrix.

      if ( lsamen( 2, C2, 'PO' ) ) {

         // CPOTRF

        srnamc.SRNAMT = 'CPOTRF';
         INFOT = 1;
         cpotrf('/', 0, A, 1, INFO );
         chkxer('CPOTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpotrf('U', -1, A, 1, INFO );
         chkxer('CPOTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cpotrf('U', 2, A, 1, INFO );
         chkxer('CPOTRF', INFOT, NOUT, LERR, OK );

         // CPOTF2

        srnamc.SRNAMT = 'CPOTF2';
         INFOT = 1;
         cpotf2('/', 0, A, 1, INFO );
         chkxer('CPOTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpotf2('U', -1, A, 1, INFO );
         chkxer('CPOTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cpotf2('U', 2, A, 1, INFO );
         chkxer('CPOTF2', INFOT, NOUT, LERR, OK );

         // CPOTRI

        srnamc.SRNAMT = 'CPOTRI';
         INFOT = 1;
         cpotri('/', 0, A, 1, INFO );
         chkxer('CPOTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpotri('U', -1, A, 1, INFO );
         chkxer('CPOTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cpotri('U', 2, A, 1, INFO );
         chkxer('CPOTRI', INFOT, NOUT, LERR, OK );

         // CPOTRS

        srnamc.SRNAMT = 'CPOTRS';
         INFOT = 1;
         cpotrs('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('CPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpotrs('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('CPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpotrs('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('CPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cpotrs('U', 2, 1, A, 1, B, 2, INFO );
         chkxer('CPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cpotrs('U', 2, 1, A, 2, B, 1, INFO );
         chkxer('CPOTRS', INFOT, NOUT, LERR, OK );

         // CPORFS

        srnamc.SRNAMT = 'CPORFS';
         INFOT = 1;
         cporfs('/', 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cporfs('U', -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cporfs('U', 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cporfs('U', 2, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cporfs('U', 2, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cporfs('U', 2, 1, A, 2, AF, 2, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cporfs('U', 2, 1, A, 2, AF, 2, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CPORFS', INFOT, NOUT, LERR, OK );

         // CPOCON

        srnamc.SRNAMT = 'CPOCON';
         INFOT = 1;
         cpocon('/', 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpocon('U', -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cpocon('U', 2, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cpocon('U', 1, A, 1, -ANRM, RCOND, W, R, INFO );
         chkxer('CPOCON', INFOT, NOUT, LERR, OK );

         // CPOEQU

        srnamc.SRNAMT = 'CPOEQU';
         INFOT = 1;
         cpoequ(-1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPOEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpoequ(2, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPOEQU', INFOT, NOUT, LERR, OK );

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite packed matrix.

      } else if ( lsamen( 2, C2, 'PP' ) ) {

         // CPPTRF

        srnamc.SRNAMT = 'CPPTRF';
         INFOT = 1;
         cpptrf('/', 0, A, INFO );
         chkxer('CPPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpptrf('U', -1, A, INFO );
         chkxer('CPPTRF', INFOT, NOUT, LERR, OK );

         // CPPTRI

        srnamc.SRNAMT = 'CPPTRI';
         INFOT = 1;
         cpptri('/', 0, A, INFO );
         chkxer('CPPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpptri('U', -1, A, INFO );
         chkxer('CPPTRI', INFOT, NOUT, LERR, OK );

         // CPPTRS

        srnamc.SRNAMT = 'CPPTRS';
         INFOT = 1;
         cpptrs('/', 0, 0, A, B, 1, INFO );
         chkxer('CPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpptrs('U', -1, 0, A, B, 1, INFO );
         chkxer('CPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpptrs('U', 0, -1, A, B, 1, INFO );
         chkxer('CPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cpptrs('U', 2, 1, A, B, 1, INFO );
         chkxer('CPPTRS', INFOT, NOUT, LERR, OK );

         // CPPRFS

        srnamc.SRNAMT = 'CPPRFS';
         INFOT = 1;
         cpprfs('/', 0, 0, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpprfs('U', -1, 0, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpprfs('U', 0, -1, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cpprfs('U', 2, 1, A, AF, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cpprfs('U', 2, 1, A, AF, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CPPRFS', INFOT, NOUT, LERR, OK );

         // CPPCON

        srnamc.SRNAMT = 'CPPCON';
         INFOT = 1;
         cppcon('/', 0, A, ANRM, RCOND, W, R, INFO );
         chkxer('CPPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cppcon('U', -1, A, ANRM, RCOND, W, R, INFO );
         chkxer('CPPCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cppcon('U', 1, A, -ANRM, RCOND, W, R, INFO );
         chkxer('CPPCON', INFOT, NOUT, LERR, OK );

         // CPPEQU

        srnamc.SRNAMT = 'CPPEQU';
         INFOT = 1;
         cppequ('/', 0, A, R1, RCOND, ANRM, INFO );
         chkxer('CPPEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cppequ('U', -1, A, R1, RCOND, ANRM, INFO );
         chkxer('CPPEQU', INFOT, NOUT, LERR, OK );

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite band matrix.

      } else if ( lsamen( 2, C2, 'PB' ) ) {

         // CPBTRF

        srnamc.SRNAMT = 'CPBTRF';
         INFOT = 1;
         cpbtrf('/', 0, 0, A, 1, INFO );
         chkxer('CPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpbtrf('U', -1, 0, A, 1, INFO );
         chkxer('CPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpbtrf('U', 1, -1, A, 1, INFO );
         chkxer('CPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cpbtrf('U', 2, 1, A, 1, INFO );
         chkxer('CPBTRF', INFOT, NOUT, LERR, OK );

         // CPBTF2

        srnamc.SRNAMT = 'CPBTF2';
         INFOT = 1;
         cpbtf2('/', 0, 0, A, 1, INFO );
         chkxer('CPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpbtf2('U', -1, 0, A, 1, INFO );
         chkxer('CPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpbtf2('U', 1, -1, A, 1, INFO );
         chkxer('CPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cpbtf2('U', 2, 1, A, 1, INFO );
         chkxer('CPBTF2', INFOT, NOUT, LERR, OK );

         // CPBTRS

        srnamc.SRNAMT = 'CPBTRS';
         INFOT = 1;
         cpbtrs('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('CPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpbtrs('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('CPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpbtrs('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('CPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cpbtrs('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('CPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cpbtrs('U', 2, 1, 1, A, 1, B, 1, INFO );
         chkxer('CPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cpbtrs('U', 2, 0, 1, A, 1, B, 1, INFO );
         chkxer('CPBTRS', INFOT, NOUT, LERR, OK );

         // CPBRFS

        srnamc.SRNAMT = 'CPBRFS';
         INFOT = 1;
         cpbrfs('/', 0, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpbrfs('U', -1, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpbrfs('U', 1, -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cpbrfs('U', 0, 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cpbrfs('U', 2, 1, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cpbrfs('U', 2, 1, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cpbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cpbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CPBRFS', INFOT, NOUT, LERR, OK );

         // CPBCON

        srnamc.SRNAMT = 'CPBCON';
         INFOT = 1;
         cpbcon('/', 0, 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpbcon('U', -1, 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpbcon('U', 1, -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cpbcon('U', 2, 1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cpbcon('U', 1, 0, A, 1, -ANRM, RCOND, W, R, INFO );
         chkxer('CPBCON', INFOT, NOUT, LERR, OK );

         // CPBEQU

        srnamc.SRNAMT = 'CPBEQU';
         INFOT = 1;
         cpbequ('/', 0, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpbequ('U', -1, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpbequ('U', 1, -1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cpbequ('U', 2, 1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPBEQU', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }