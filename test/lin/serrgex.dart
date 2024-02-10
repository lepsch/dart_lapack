      void serrge(PATH, final int NUNIT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX, LW;
      const              NMAX = 4, LW = 3*NMAX ;
      String             EQ;
      String             C2;
      int                I, INFO, J, N_ERR_BNDS, NPARAMS;
      double               ANRM, CCOND, RCOND, BERR;
      int                IP( NMAX ), IW( NMAX );
      double               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), W( LW ), X( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGBCON, SGBEQU, SGBRFS, SGBTF2, SGBTRF, SGBTRS, SGECON, SGEEQU, SGERFS, SGETF2, SGETRF, SGETRI, SGETRS, SGEEQUB, SGERFSX, SGBEQUB, SGBRFSX
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
      // INTRINSIC REAL

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = 1. / REAL( I+J );
            AF[I][J] = 1. / REAL( I+J );
         } // 10
         B[J] = 0.;
         R1[J] = 0.;
         R2[J] = 0.;
         W[J] = 0.;
         X[J] = 0.;
         C[J] = 0.;
         R[J] = 0.;
         IP[J] = J;
         IW[J] = J;
      } // 20
      OK = true;

      if ( lsamen( 2, C2, 'GE' ) ) {

         // Test error exits of the routines that use the LU decomposition
         // of a general matrix.

         // SGETRF

        srnamc.SRNAMT = 'SGETRF';
         INFOT = 1;
         sgetrf(-1, 0, A, 1, IP, INFO );
         chkxer('SGETRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgetrf(0, -1, A, 1, IP, INFO );
         chkxer('SGETRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgetrf(2, 1, A, 1, IP, INFO );
         chkxer('SGETRF', INFOT, NOUT, LERR, OK );

         // SGETF2

        srnamc.SRNAMT = 'SGETF2';
         INFOT = 1;
         sgetf2(-1, 0, A, 1, IP, INFO );
         chkxer('SGETF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgetf2(0, -1, A, 1, IP, INFO );
         chkxer('SGETF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgetf2(2, 1, A, 1, IP, INFO );
         chkxer('SGETF2', INFOT, NOUT, LERR, OK );

         // SGETRI

        srnamc.SRNAMT = 'SGETRI';
         INFOT = 1;
         sgetri(-1, A, 1, IP, W, LW, INFO );
         chkxer('SGETRI', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgetri(2, A, 1, IP, W, LW, INFO );
         chkxer('SGETRI', INFOT, NOUT, LERR, OK );

         // SGETRS

        srnamc.SRNAMT = 'SGETRS';
         INFOT = 1;
         sgetrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('SGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgetrs('N', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('SGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgetrs('N', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('SGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgetrs('N', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('SGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgetrs('N', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('SGETRS', INFOT, NOUT, LERR, OK );

         // SGERFS

        srnamc.SRNAMT = 'SGERFS';
         INFOT = 1;
         sgerfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgerfs('N', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgerfs('N', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgerfs('N', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgerfs('N', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgerfs('N', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('SGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sgerfs('N', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGERFS', INFOT, NOUT, LERR, OK );

         // SGERFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
        srnamc.SRNAMT = 'SGERFSX';
         INFOT = 1;
         sgerfsx('/', EQ, 0, 0, A, 1, AF, 1, IP, R, C, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS,  W, IW, INFO );
         chkxer('SGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         EQ = '/';
         sgerfsx('N', EQ, 2, 1, A, 1, AF, 2, IP, R, C, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         EQ = 'R';
         sgerfsx('N', EQ, -1, 0, A, 1, AF, 1, IP, R, C, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgerfsx('N', EQ, 0, -1, A, 1, AF, 1, IP, R, C, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgerfsx('N', EQ, 2, 1, A, 1, AF, 2, IP, R, C, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgerfsx('N', EQ, 2, 1, A, 2, AF, 1, IP, R, C, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'C';
         sgerfsx('N', EQ, 2, 1, A, 2, AF, 2, IP, R, C, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgerfsx('N', EQ, 2, 1, A, 2, AF, 2, IP, R, C, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGERFSX', INFOT, NOUT, LERR, OK );

         // SGECON

        srnamc.SRNAMT = 'SGECON';
         INFOT = 1;
         sgecon('/', 0, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SGECON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgecon('1', -1, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SGECON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgecon('1', 2, A, 1, ANRM, RCOND, W, IW, INFO );
         chkxer('SGECON', INFOT, NOUT, LERR, OK );

         // SGEEQU

        srnamc.SRNAMT = 'SGEEQU';
         INFOT = 1;
         sgeequ(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGEEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgeequ(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGEEQU', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgeequ(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGEEQU', INFOT, NOUT, LERR, OK );

         // SGEEQUB

        srnamc.SRNAMT = 'SGEEQUB';
         INFOT = 1;
         sgeequb(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGEEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgeequb(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGEEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgeequb(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGEEQUB', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'GB' ) ) {

         // Test error exits of the routines that use the LU decomposition
         // of a general band matrix.

         // SGBTRF

        srnamc.SRNAMT = 'SGBTRF';
         INFOT = 1;
         sgbtrf(-1, 0, 0, 0, A, 1, IP, INFO );
         chkxer('SGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgbtrf(0, -1, 0, 0, A, 1, IP, INFO );
         chkxer('SGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgbtrf(1, 1, -1, 0, A, 1, IP, INFO );
         chkxer('SGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgbtrf(1, 1, 0, -1, A, 1, IP, INFO );
         chkxer('SGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgbtrf(2, 2, 1, 1, A, 3, IP, INFO );
         chkxer('SGBTRF', INFOT, NOUT, LERR, OK );

         // SGBTF2

        srnamc.SRNAMT = 'SGBTF2';
         INFOT = 1;
         sgbtf2(-1, 0, 0, 0, A, 1, IP, INFO );
         chkxer('SGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgbtf2(0, -1, 0, 0, A, 1, IP, INFO );
         chkxer('SGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgbtf2(1, 1, -1, 0, A, 1, IP, INFO );
         chkxer('SGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgbtf2(1, 1, 0, -1, A, 1, IP, INFO );
         chkxer('SGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgbtf2(2, 2, 1, 1, A, 3, IP, INFO );
         chkxer('SGBTF2', INFOT, NOUT, LERR, OK );

         // SGBTRS

        srnamc.SRNAMT = 'SGBTRS';
         INFOT = 1;
         sgbtrs('/', 0, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('SGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgbtrs('N', -1, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('SGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgbtrs('N', 1, -1, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('SGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgbtrs('N', 1, 0, -1, 1, A, 1, IP, B, 1, INFO );
         chkxer('SGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgbtrs('N', 1, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('SGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgbtrs('N', 2, 1, 1, 1, A, 3, IP, B, 2, INFO );
         chkxer('SGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgbtrs('N', 2, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('SGBTRS', INFOT, NOUT, LERR, OK );

         // SGBRFS

        srnamc.SRNAMT = 'SGBRFS';
         INFOT = 1;
         sgbrfs('/', 0, 0, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgbrfs('N', -1, 0, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgbrfs('N', 1, -1, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgbrfs('N', 1, 0, -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgbrfs('N', 1, 0, 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgbrfs('N', 2, 1, 1, 1, A, 2, AF, 4, IP, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgbrfs('N', 2, 1, 1, 1, A, 3, AF, 3, IP, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGBRFS', INFOT, NOUT, LERR, OK );

         // SGBRFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
        srnamc.SRNAMT = 'SGBRFSX';
         INFOT = 1;
         sgbrfsx('/', EQ, 0, 0, 0, 0, A, 1, AF, 1, IP, R, C, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS,  W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         EQ = '/';
         sgbrfsx('N', EQ, 2, 1, 1, 1, A, 1, AF, 2, IP, R, C, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         EQ = 'R';
         sgbrfsx('N', EQ, -1, 1, 1, 0, A, 1, AF, 1, IP, R, C, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         EQ = 'R';
         sgbrfsx('N', EQ, 2, -1, 1, 1, A, 3, AF, 4, IP, R, C, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         EQ = 'R';
         sgbrfsx('N', EQ, 2, 1, -1, 1, A, 3, AF, 4, IP, R, C, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgbrfsx('N', EQ, 0, 0, 0, -1, A, 1, AF, 1, IP, R, C, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgbrfsx('N', EQ, 2, 1, 1, 1, A, 1, AF, 2, IP, R, C, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgbrfsx('N', EQ, 2, 1, 1, 1, A, 3, AF, 3, IP, R, C, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'C';
         sgbrfsx('N', EQ, 2, 1, 1, 1, A, 3, AF, 5, IP, R, C, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgbrfsx('N', EQ, 2, 1, 1, 1, A, 3, AF, 5, IP, R, C, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SGBRFSX', INFOT, NOUT, LERR, OK );

         // SGBCON

        srnamc.SRNAMT = 'SGBCON';
         INFOT = 1;
         sgbcon('/', 0, 0, 0, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgbcon('1', -1, 0, 0, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgbcon('1', 1, -1, 0, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgbcon('1', 1, 0, -1, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgbcon('1', 2, 1, 1, A, 3, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SGBCON', INFOT, NOUT, LERR, OK );

         // SGBEQU

        srnamc.SRNAMT = 'SGBEQU';
         INFOT = 1;
         sgbequ(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgbequ(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgbequ(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgbequ(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgbequ(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQU', INFOT, NOUT, LERR, OK );

         // SGBEQUB

        srnamc.SRNAMT = 'SGBEQUB';
         INFOT = 1;
         sgbequb(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgbequb(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgbequb(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgbequb(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQUB', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgbequb(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('SGBEQUB', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }
