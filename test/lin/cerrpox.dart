      void cerrpo(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      String             EQ;
      String             C2;
      int                I, INFO, J, N_ERR_BNDS, NPARAMS;
      double               ANRM, RCOND, BERR;
      // ..
      // .. Local Arrays ..
      double               S( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 );
      Complex            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CPBCON, CPBEQU, CPBRFS, CPBTF2, CPBTRF, CPBTRS, CPOCON, CPOEQU, CPORFS, CPOTF2, CPOTRF, CPOTRI, CPOTRS, CPPCON, CPPEQU, CPPRFS, CPPTRF, CPPTRI, CPPTRS, CPOEQUB, CPORFSX
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
         } // 10
         B[J] = 0.;
         R1[J] = 0.;
         R2[J] = 0.;
         W[J] = 0.;
         X[J] = 0.;
         S[J] = 0.;
      } // 20
      ANRM = 1.;
      OK = true;

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite matrix.

      if ( LSAMEN( 2, C2, 'PO' ) ) {

         // CPOTRF

         SRNAMT = 'CPOTRF';
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

         SRNAMT = 'CPOTF2';
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

         SRNAMT = 'CPOTRI';
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

         SRNAMT = 'CPOTRS';
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

         SRNAMT = 'CPORFS';
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

         // CPORFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
         SRNAMT = 'CPORFSX';
         INFOT = 1;
         cporfsx('/', EQ, 0, 0, A, 1, AF, 1, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cporfsx('U', '/', -1, 0, A, 1, AF, 1, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CPORFSX', INFOT, NOUT, LERR, OK );
         EQ = 'N';
         INFOT = 3;
         cporfsx('U', EQ, -1, 0, A, 1, AF, 1, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cporfsx('U', EQ, 0, -1, A, 1, AF, 1, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cporfsx('U', EQ, 2, 1, A, 1, AF, 2, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cporfsx('U', EQ, 2, 1, A, 2, AF, 1, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cporfsx('U', EQ, 2, 1, A, 2, AF, 2, S, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cporfsx('U', EQ, 2, 1, A, 2, AF, 2, S, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CPORFSX', INFOT, NOUT, LERR, OK );

         // CPOCON

         SRNAMT = 'CPOCON';
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

         SRNAMT = 'CPOEQU';
         INFOT = 1;
         cpoequ(-1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPOEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpoequ(2, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPOEQU', INFOT, NOUT, LERR, OK );

         // CPOEQUB

         SRNAMT = 'CPOEQUB';
         INFOT = 1;
         cpoequb(-1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPOEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpoequb(2, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('CPOEQUB', INFOT, NOUT, LERR, OK );

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite packed matrix.

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // CPPTRF

         SRNAMT = 'CPPTRF';
         INFOT = 1;
         cpptrf('/', 0, A, INFO );
         chkxer('CPPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpptrf('U', -1, A, INFO );
         chkxer('CPPTRF', INFOT, NOUT, LERR, OK );

         // CPPTRI

         SRNAMT = 'CPPTRI';
         INFOT = 1;
         cpptri('/', 0, A, INFO );
         chkxer('CPPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpptri('U', -1, A, INFO );
         chkxer('CPPTRI', INFOT, NOUT, LERR, OK );

         // CPPTRS

         SRNAMT = 'CPPTRS';
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

         SRNAMT = 'CPPRFS';
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

         SRNAMT = 'CPPCON';
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

         SRNAMT = 'CPPEQU';
         INFOT = 1;
         cppequ('/', 0, A, R1, RCOND, ANRM, INFO );
         chkxer('CPPEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cppequ('U', -1, A, R1, RCOND, ANRM, INFO );
         chkxer('CPPEQU', INFOT, NOUT, LERR, OK );

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite band matrix.

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // CPBTRF

         SRNAMT = 'CPBTRF';
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

         SRNAMT = 'CPBTF2';
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

         SRNAMT = 'CPBTRS';
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

         SRNAMT = 'CPBRFS';
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

         SRNAMT = 'CPBCON';
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

         SRNAMT = 'CPBEQU';
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

      return;
      }
