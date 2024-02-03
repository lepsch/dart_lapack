      SUBROUTINE DERRPO( PATH, NUNIT );

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
      double             ANRM, RCOND, BERR;
      // ..
      // .. Local Arrays ..
      int                IW( NMAX );
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), R1( NMAX ), R2( NMAX ), W( 3*NMAX ), X( NMAX ), S( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3), PARAMS( 1 );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DPBCON, DPBEQU, DPBRFS, DPBTF2, DPBTRF, DPBTRS, DPOCON, DPOEQU, DPORFS, DPOTF2, DPOTRF, DPOTRI, DPOTRS, DPPCON, DPPEQU, DPPRFS, DPPTRF, DPPTRI, DPPTRS, DPOEQUB, DPORFSX
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
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1.0 / DBLE( I+J );
            AF( I, J ) = 1.0 / DBLE( I+J );
         } // 10
         B( J ) = 0.0;
         R1( J ) = 0.0;
         R2( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
         S( J ) = 0.0;
         IW( J ) = J;
      } // 20
      OK = true;

      if ( LSAMEN( 2, C2, 'PO' ) ) {

         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive definite matrix.

         // DPOTRF

         SRNAMT = 'DPOTRF';
         INFOT = 1;
         dpotrf('/', 0, A, 1, INFO );
         chkxer('DPOTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpotrf('U', -1, A, 1, INFO );
         chkxer('DPOTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dpotrf('U', 2, A, 1, INFO );
         chkxer('DPOTRF', INFOT, NOUT, LERR, OK );

         // DPOTF2

         SRNAMT = 'DPOTF2';
         INFOT = 1;
         dpotf2('/', 0, A, 1, INFO );
         chkxer('DPOTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpotf2('U', -1, A, 1, INFO );
         chkxer('DPOTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dpotf2('U', 2, A, 1, INFO );
         chkxer('DPOTF2', INFOT, NOUT, LERR, OK );

         // DPOTRI

         SRNAMT = 'DPOTRI';
         INFOT = 1;
         dpotri('/', 0, A, 1, INFO );
         chkxer('DPOTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpotri('U', -1, A, 1, INFO );
         chkxer('DPOTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dpotri('U', 2, A, 1, INFO );
         chkxer('DPOTRI', INFOT, NOUT, LERR, OK );

         // DPOTRS

         SRNAMT = 'DPOTRS';
         INFOT = 1;
         dpotrs('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('DPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpotrs('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('DPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpotrs('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('DPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dpotrs('U', 2, 1, A, 1, B, 2, INFO );
         chkxer('DPOTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dpotrs('U', 2, 1, A, 2, B, 1, INFO );
         chkxer('DPOTRS', INFOT, NOUT, LERR, OK );

         // DPORFS

         SRNAMT = 'DPORFS';
         INFOT = 1;
         dporfs('/', 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dporfs('U', -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dporfs('U', 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dporfs('U', 2, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('DPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dporfs('U', 2, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('DPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dporfs('U', 2, 1, A, 2, AF, 2, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DPORFS', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dporfs('U', 2, 1, A, 2, AF, 2, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPORFS', INFOT, NOUT, LERR, OK );

         // DPORFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
         SRNAMT = 'DPORFSX';
         INFOT = 1;
         dporfsx('/', EQ, 0, 0, A, 1, AF, 1, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dporfsx('U', "/", -1, 0, A, 1, AF, 1, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPORFSX', INFOT, NOUT, LERR, OK );
         EQ = 'N';
         INFOT = 3;
         dporfsx('U', EQ, -1, 0, A, 1, AF, 1, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dporfsx('U', EQ, 0, -1, A, 1, AF, 1, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dporfsx('U', EQ, 2, 1, A, 1, AF, 2, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dporfsx('U', EQ, 2, 1, A, 2, AF, 1, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dporfsx('U', EQ, 2, 1, A, 2, AF, 2, S, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPORFSX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dporfsx('U', EQ, 2, 1, A, 2, AF, 2, S, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPORFSX', INFOT, NOUT, LERR, OK );

         // DPOCON

         SRNAMT = 'DPOCON';
         INFOT = 1;
         dpocon('/', 0, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('DPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpocon('U', -1, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('DPOCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dpocon('U', 2, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('DPOCON', INFOT, NOUT, LERR, OK );

         // DPOEQU

         SRNAMT = 'DPOEQU';
         INFOT = 1;
         dpoequ(-1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('DPOEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpoequ(2, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('DPOEQU', INFOT, NOUT, LERR, OK );

         // DPOEQUB

         SRNAMT = 'DPOEQUB';
         INFOT = 1;
         dpoequb(-1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('DPOEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpoequb(2, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('DPOEQUB', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive definite packed matrix.

         // DPPTRF

         SRNAMT = 'DPPTRF';
         INFOT = 1;
         dpptrf('/', 0, A, INFO );
         chkxer('DPPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpptrf('U', -1, A, INFO );
         chkxer('DPPTRF', INFOT, NOUT, LERR, OK );

         // DPPTRI

         SRNAMT = 'DPPTRI';
         INFOT = 1;
         dpptri('/', 0, A, INFO );
         chkxer('DPPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpptri('U', -1, A, INFO );
         chkxer('DPPTRI', INFOT, NOUT, LERR, OK );

         // DPPTRS

         SRNAMT = 'DPPTRS';
         INFOT = 1;
         dpptrs('/', 0, 0, A, B, 1, INFO );
         chkxer('DPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpptrs('U', -1, 0, A, B, 1, INFO );
         chkxer('DPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpptrs('U', 0, -1, A, B, 1, INFO );
         chkxer('DPPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dpptrs('U', 2, 1, A, B, 1, INFO );
         chkxer('DPPTRS', INFOT, NOUT, LERR, OK );

         // DPPRFS

         SRNAMT = 'DPPRFS';
         INFOT = 1;
         dpprfs('/', 0, 0, A, AF, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpprfs('U', -1, 0, A, AF, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpprfs('U', 0, -1, A, AF, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dpprfs('U', 2, 1, A, AF, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DPPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dpprfs('U', 2, 1, A, AF, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPPRFS', INFOT, NOUT, LERR, OK );

         // DPPCON

         SRNAMT = 'DPPCON';
         INFOT = 1;
         dppcon('/', 0, A, ANRM, RCOND, W, IW, INFO );
         chkxer('DPPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dppcon('U', -1, A, ANRM, RCOND, W, IW, INFO );
         chkxer('DPPCON', INFOT, NOUT, LERR, OK );

         // DPPEQU

         SRNAMT = 'DPPEQU';
         INFOT = 1;
         dppequ('/', 0, A, R1, RCOND, ANRM, INFO );
         chkxer('DPPEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dppequ('U', -1, A, R1, RCOND, ANRM, INFO );
         chkxer('DPPEQU', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive definite band matrix.

         // DPBTRF

         SRNAMT = 'DPBTRF';
         INFOT = 1;
         dpbtrf('/', 0, 0, A, 1, INFO );
         chkxer('DPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpbtrf('U', -1, 0, A, 1, INFO );
         chkxer('DPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpbtrf('U', 1, -1, A, 1, INFO );
         chkxer('DPBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dpbtrf('U', 2, 1, A, 1, INFO );
         chkxer('DPBTRF', INFOT, NOUT, LERR, OK );

         // DPBTF2

         SRNAMT = 'DPBTF2';
         INFOT = 1;
         dpbtf2('/', 0, 0, A, 1, INFO );
         chkxer('DPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpbtf2('U', -1, 0, A, 1, INFO );
         chkxer('DPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpbtf2('U', 1, -1, A, 1, INFO );
         chkxer('DPBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dpbtf2('U', 2, 1, A, 1, INFO );
         chkxer('DPBTF2', INFOT, NOUT, LERR, OK );

         // DPBTRS

         SRNAMT = 'DPBTRS';
         INFOT = 1;
         dpbtrs('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('DPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpbtrs('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('DPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpbtrs('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('DPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dpbtrs('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('DPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dpbtrs('U', 2, 1, 1, A, 1, B, 1, INFO );
         chkxer('DPBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dpbtrs('U', 2, 0, 1, A, 1, B, 1, INFO );
         chkxer('DPBTRS', INFOT, NOUT, LERR, OK );

         // DPBRFS

         SRNAMT = 'DPBRFS';
         INFOT = 1;
         dpbrfs('/', 0, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpbrfs('U', -1, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpbrfs('U', 1, -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dpbrfs('U', 0, 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dpbrfs('U', 2, 1, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('DPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dpbrfs('U', 2, 1, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('DPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dpbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DPBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dpbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DPBRFS', INFOT, NOUT, LERR, OK );

         // DPBCON

         SRNAMT = 'DPBCON';
         INFOT = 1;
         dpbcon('/', 0, 0, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('DPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpbcon('U', -1, 0, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('DPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpbcon('U', 1, -1, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('DPBCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dpbcon('U', 2, 1, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('DPBCON', INFOT, NOUT, LERR, OK );

         // DPBEQU

         SRNAMT = 'DPBEQU';
         INFOT = 1;
         dpbequ('/', 0, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('DPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpbequ('U', -1, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('DPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpbequ('U', 1, -1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('DPBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dpbequ('U', 2, 1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('DPBEQU', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of DERRPOX

      }
