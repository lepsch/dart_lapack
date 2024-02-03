      SUBROUTINE ZERRGE( PATH, NUNIT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      String             EQ;
      String             C2;
      int                I, INFO, J, N_ERR_BNDS, NPARAMS;
      double             ANRM, CCOND, RCOND, BERR;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             R( NMAX ), R1( NMAX ), R2( NMAX ), CS( NMAX ), RS( NMAX );
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( 2*NMAX ), X( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS;
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGBCON, ZGBEQU, ZGBRFS, ZGBTF2, ZGBTRF, ZGBTRS, ZGECON, ZGEEQU, ZGERFS, ZGETF2, ZGETRF, ZGETRI, ZGETRS, ZGEEQUB, ZGERFSX, ZGBEQUB, ZGBRFSX
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
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) );
         } // 10
         B( J ) = 0.0;
         R1( J ) = 0.0;
         R2( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
         CS( J ) = 0.0;
         RS( J ) = 0.0;
         IP( J ) = J;
      } // 20
      OK = true;

      // Test error exits of the routines that use the LU decomposition
      // of a general matrix.

      if ( LSAMEN( 2, C2, 'GE' ) ) {

         // ZGETRF

         SRNAMT = 'ZGETRF';
         INFOT = 1;
         zgetrf(-1, 0, A, 1, IP, INFO );
         chkxer('ZGETRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgetrf(0, -1, A, 1, IP, INFO );
         chkxer('ZGETRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgetrf(2, 1, A, 1, IP, INFO );
         chkxer('ZGETRF', INFOT, NOUT, LERR, OK );

         // ZGETF2

         SRNAMT = 'ZGETF2';
         INFOT = 1;
         zgetf2(-1, 0, A, 1, IP, INFO );
         chkxer('ZGETF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgetf2(0, -1, A, 1, IP, INFO );
         chkxer('ZGETF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgetf2(2, 1, A, 1, IP, INFO );
         chkxer('ZGETF2', INFOT, NOUT, LERR, OK );

         // ZGETRI

         SRNAMT = 'ZGETRI';
         INFOT = 1;
         zgetri(-1, A, 1, IP, W, 1, INFO );
         chkxer('ZGETRI', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgetri(2, A, 1, IP, W, 2, INFO );
         chkxer('ZGETRI', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgetri(2, A, 2, IP, W, 1, INFO );
         chkxer('ZGETRI', INFOT, NOUT, LERR, OK );

         // ZGETRS

         SRNAMT = 'ZGETRS';
         INFOT = 1;
         zgetrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgetrs('N', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgetrs('N', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgetrs('N', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zgetrs('N', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZGETRS', INFOT, NOUT, LERR, OK );

         // ZGERFS

         SRNAMT = 'ZGERFS';
         INFOT = 1;
         zgerfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgerfs('N', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgerfs('N', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgerfs('N', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgerfs('N', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgerfs('N', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zgerfs('N', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGERFS', INFOT, NOUT, LERR, OK );

         // ZGERFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
         SRNAMT = 'ZGERFSX';
         INFOT = 1;
         zgerfsx('/', EQ, 0, 0, A, 1, AF, 1, IP, RS, CS, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         EQ = '/';
         zgerfsx('N', EQ, 2, 1, A, 1, AF, 2, IP, RS, CS, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         EQ = 'R';
         zgerfsx('N', EQ, -1, 0, A, 1, AF, 1, IP, RS, CS, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgerfsx('N', EQ, 0, -1, A, 1, AF, 1, IP, RS, CS, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgerfsx('N', EQ, 2, 1, A, 1, AF, 2, IP, RS, CS, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zgerfsx('N', EQ, 2, 1, A, 2, AF, 1, IP, RS, CS, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'C';
         zgerfsx('N', EQ, 2, 1, A, 2, AF, 2, IP, RS, CS, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zgerfsx('N', EQ, 2, 1, A, 2, AF, 2, IP, RS, CS, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGERFSX', INFOT, NOUT, LERR, OK );

         // ZGECON

         SRNAMT = 'ZGECON';
         INFOT = 1;
         zgecon('/', 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZGECON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgecon('1', -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZGECON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgecon('1', 2, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZGECON', INFOT, NOUT, LERR, OK );

         // ZGEEQU

         SRNAMT = 'ZGEEQU';
         INFOT = 1;
         zgeequ(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgeequ(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQU', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgeequ(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQU', INFOT, NOUT, LERR, OK );

         // ZGEEQUB

         SRNAMT = 'ZGEEQUB';
         INFOT = 1;
         zgeequb(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgeequb(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgeequb(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQUB', INFOT, NOUT, LERR, OK );

      // Test error exits of the routines that use the LU decomposition
      // of a general band matrix.

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // ZGBTRF

         SRNAMT = 'ZGBTRF';
         INFOT = 1;
         zgbtrf(-1, 0, 0, 0, A, 1, IP, INFO );
         chkxer('ZGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgbtrf(0, -1, 0, 0, A, 1, IP, INFO );
         chkxer('ZGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgbtrf(1, 1, -1, 0, A, 1, IP, INFO );
         chkxer('ZGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgbtrf(1, 1, 0, -1, A, 1, IP, INFO );
         chkxer('ZGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgbtrf(2, 2, 1, 1, A, 3, IP, INFO );
         chkxer('ZGBTRF', INFOT, NOUT, LERR, OK );

         // ZGBTF2

         SRNAMT = 'ZGBTF2';
         INFOT = 1;
         zgbtf2(-1, 0, 0, 0, A, 1, IP, INFO );
         chkxer('ZGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgbtf2(0, -1, 0, 0, A, 1, IP, INFO );
         chkxer('ZGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgbtf2(1, 1, -1, 0, A, 1, IP, INFO );
         chkxer('ZGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgbtf2(1, 1, 0, -1, A, 1, IP, INFO );
         chkxer('ZGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgbtf2(2, 2, 1, 1, A, 3, IP, INFO );
         chkxer('ZGBTF2', INFOT, NOUT, LERR, OK );

         // ZGBTRS

         SRNAMT = 'ZGBTRS';
         INFOT = 1;
         zgbtrs('/', 0, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgbtrs('N', -1, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgbtrs('N', 1, -1, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgbtrs('N', 1, 0, -1, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgbtrs('N', 1, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgbtrs('N', 2, 1, 1, 1, A, 3, IP, B, 2, INFO );
         chkxer('ZGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgbtrs('N', 2, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', INFOT, NOUT, LERR, OK );

         // ZGBRFS

         SRNAMT = 'ZGBRFS';
         INFOT = 1;
         zgbrfs('/', 0, 0, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgbrfs('N', -1, 0, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgbrfs('N', 1, -1, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgbrfs('N', 1, 0, -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgbrfs('N', 1, 0, 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgbrfs('N', 2, 1, 1, 1, A, 2, AF, 4, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zgbrfs('N', 2, 1, 1, 1, A, 3, AF, 3, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         zgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', INFOT, NOUT, LERR, OK );

         // ZGBRFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
         SRNAMT = 'ZGBRFSX';
         INFOT = 1;
         zgbrfsx('/', EQ, 0, 0, 0, 0, A, 1, AF, 1, IP, RS, CS, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS,  W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         EQ = '/';
         zgbrfsx('N', EQ, 2, 1, 1, 1, A, 1, AF, 2, IP, RS, CS, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         EQ = 'R';
         zgbrfsx('N', EQ, -1, 1, 1, 0, A, 1, AF, 1, IP, RS, CS, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         EQ = 'R';
         zgbrfsx('N', EQ, 2, -1, 1, 1, A, 3, AF, 4, IP, RS, CS, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         EQ = 'R';
         zgbrfsx('N', EQ, 2, 1, -1, 1, A, 3, AF, 4, IP, RS, CS, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgbrfsx('N', EQ, 0, 0, 0, -1, A, 1, AF, 1, IP, RS, CS, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zgbrfsx('N', EQ, 2, 1, 1, 1, A, 1, AF, 2, IP, RS, CS, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgbrfsx('N', EQ, 2, 1, 1, 1, A, 3, AF, 3, IP, RS, CS, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'C';
         zgbrfsx('N', EQ, 2, 1, 1, 1, A, 3, AF, 5, IP, RS, CS, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zgbrfsx('N', EQ, 2, 1, 1, 1, A, 3, AF, 5, IP, RS, CS, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('ZGBRFSX', INFOT, NOUT, LERR, OK );

         // ZGBCON

         SRNAMT = 'ZGBCON';
         INFOT = 1;
         zgbcon('/', 0, 0, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgbcon('1', -1, 0, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgbcon('1', 1, -1, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgbcon('1', 1, 0, -1, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgbcon('1', 2, 1, 1, A, 3, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', INFOT, NOUT, LERR, OK );

         // ZGBEQU

         SRNAMT = 'ZGBEQU';
         INFOT = 1;
         zgbequ(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgbequ(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgbequ(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgbequ(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgbequ(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', INFOT, NOUT, LERR, OK );

         // ZGBEQUB

         SRNAMT = 'ZGBEQUB';
         INFOT = 1;
         zgbequb(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgbequb(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgbequb(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgbequb(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgbequb(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQUB', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of ZERRGEX

      }
