      void cerrvx(PATH, NUNIT ) {

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
      double               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      String             EQ;
      String             C2;
      int                I, INFO, J, N_ERR_BNDS, NPARAMS;
      double               RCOND, RPVGRW, BERR;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double               C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), RF( NMAX ), RW( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 );
      Complex            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGBSV, CGBSVX, CGESV, CGESVX, CGTSV, CGTSVX, CHESV, CHESV_RK, CHESV_ROOK, CHESVX, CHKXER, CHPSV, CHPSVX, CPBSV, CPBSVX, CPOSV, CPOSVX, CPPSV, CPPSVX, CPTSV, CPTSVX, CSPSV, CSPSVX, CSYSV, CSYSV_RK, CSYSV_ROOK, CSYSVX, CGESVXX, CPOSVXX, CSYSVXX, CHESVXX, CGBSVXX
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
         B[J] = 0.0;
         E[J] = 0E+0;
         R1[J] = 0.0;
         R2[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
         C[J] = 0.0;
         R[J] = 0.0;
         IP[J] = J;
      } // 20
      EQ = ' ';
      OK = true;

      if ( lsamen( 2, C2, 'GE' ) ) {

         // CGESV

        srnamc.SRNAMT = 'CGESV ';
         INFOT = 1;
         cgesv(-1, 0, A, 1, IP, B, 1, INFO );
         chkxer('CGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgesv(0, -1, A, 1, IP, B, 1, INFO );
         chkxer('CGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgesv(2, 1, A, 1, IP, B, 2, INFO );
         chkxer('CGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgesv(2, 1, A, 2, IP, B, 1, INFO );
         chkxer('CGESV ', INFOT, NOUT, LERR, OK );

         // CGESVX

        srnamc.SRNAMT = 'CGESVX';
         INFOT = 1;
         cgesvx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgesvx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgesvx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgesvx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgesvx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgesvx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = '/';
         cgesvx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ = 'R';
         cgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         EQ = 'C';
         cgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGESVX', INFOT, NOUT, LERR, OK );

         // CGESVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
        srnamc.SRNAMT = 'CGESVXX';
         INFOT = 1;
         cgesvxx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgesvxx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgesvxx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgesvxx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgesvxx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgesvxx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = '/';
         cgesvxx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ = 'R';
         cgesvxx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         EQ = 'C';
         cgesvxx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgesvxx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgesvxx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGESVXX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'GB' ) ) {

         // CGBSV

        srnamc.SRNAMT = 'CGBSV ';
         INFOT = 1;
         cgbsv(-1, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('CGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbsv(1, -1, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('CGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbsv(1, 0, -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('CGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbsv(0, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('CGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgbsv(1, 1, 1, 0, A, 3, IP, B, 1, INFO );
         chkxer('CGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgbsv(2, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('CGBSV ', INFOT, NOUT, LERR, OK );

         // CGBSVX

        srnamc.SRNAMT = 'CGBSVX';
         INFOT = 1;
         cgbsvx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbsvx('N', '/', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbsvx('N', 'N', -1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbsvx('N', 'N', 1, -1, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgbsvx('N', 'N', 1, 0, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgbsvx('N', 'N', 0, 0, 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgbsvx('N', 'N', 1, 1, 1, 0, A, 2, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgbsvx('N', 'N', 1, 1, 1, 0, A, 3, AF, 3, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         EQ = '/';
         cgbsvx('F', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'R';
         cgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         EQ = 'C';
         cgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         cgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGBSVX', INFOT, NOUT, LERR, OK );

         // CGBSVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
        srnamc.SRNAMT = 'CGBSVXX';
         INFOT = 1;
         cgbsvxx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbsvxx('N', '/', 0, 1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbsvxx('N', 'N', -1, 1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbsvxx('N', 'N', 2, -1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgbsvxx('N', 'N', 2, 1, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgbsvxx('N', 'N', 0, 1, 1, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgbsvxx('N', 'N', 2, 1, 1, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 3, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         EQ = '/';
         cgbsvxx('F', 'N', 0, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'R';
         cgbsvxx('F', 'N', 1, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         EQ = 'C';
         cgbsvxx('F', 'N', 1, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 4, IP, EQ, R, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CGBSVXX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'GT' ) ) {

         // CGTSV

        srnamc.SRNAMT = 'CGTSV ';
         INFOT = 1;
         cgtsv(-1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('CGTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgtsv(0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('CGTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgtsv(2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('CGTSV ', INFOT, NOUT, LERR, OK );

         // CGTSVX

        srnamc.SRNAMT = 'CGTSVX';
         INFOT = 1;
         cgtsvx('/', 'N', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgtsvx('N', '/', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgtsvx('N', 'N', -1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgtsvx('N', 'N', 0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CGTSVX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'PO' ) ) {

         // CPOSV

        srnamc.SRNAMT = 'CPOSV ';
         INFOT = 1;
         cposv('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('CPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cposv('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('CPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cposv('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('CPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cposv('U', 2, 0, A, 1, B, 2, INFO );
         chkxer('CPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cposv('U', 2, 0, A, 2, B, 1, INFO );
         chkxer('CPOSV ', INFOT, NOUT, LERR, OK );

         // CPOSVX

        srnamc.SRNAMT = 'CPOSVX';
         INFOT = 1;
         cposvx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cposvx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cposvx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cposvx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cposvx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cposvx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         EQ = '/';
         cposvx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = 'Y';
         cposvx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPOSVX', INFOT, NOUT, LERR, OK );

         // CPOSVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
        srnamc.SRNAMT = 'CPOSVXX';
         INFOT = 1;
         cposvxx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cposvxx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cposvxx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cposvxx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cposvxx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cposvxx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         EQ = '/';
         cposvxx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = 'Y';
         cposvxx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cposvxx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cposvxx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CPOSVXX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'PP' ) ) {

         // CPPSV

        srnamc.SRNAMT = 'CPPSV ';
         INFOT = 1;
         cppsv('/', 0, 0, A, B, 1, INFO );
         chkxer('CPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cppsv('U', -1, 0, A, B, 1, INFO );
         chkxer('CPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cppsv('U', 0, -1, A, B, 1, INFO );
         chkxer('CPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cppsv('U', 2, 0, A, B, 1, INFO );
         chkxer('CPPSV ', INFOT, NOUT, LERR, OK );

         // CPPSVX

        srnamc.SRNAMT = 'CPPSVX';
         INFOT = 1;
         cppsvx('/', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cppsvx('N', '/', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cppsvx('N', 'U', -1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cppsvx('N', 'U', 0, -1, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         EQ = '/';
         cppsvx('F', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         EQ = 'Y';
         cppsvx('F', 'U', 1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPPSVX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'PB' ) ) {

         // CPBSV

        srnamc.SRNAMT = 'CPBSV ';
         INFOT = 1;
         cpbsv('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('CPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpbsv('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('CPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpbsv('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('CPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cpbsv('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('CPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cpbsv('U', 1, 1, 0, A, 1, B, 2, INFO );
         chkxer('CPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cpbsv('U', 2, 0, 0, A, 1, B, 1, INFO );
         chkxer('CPBSV ', INFOT, NOUT, LERR, OK );

         // CPBSVX

        srnamc.SRNAMT = 'CPBSVX';
         INFOT = 1;
         cpbsvx('/', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpbsvx('N', '/', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cpbsvx('N', 'U', -1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cpbsvx('N', 'U', 1, -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cpbsvx('N', 'U', 0, 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cpbsvx('N', 'U', 1, 1, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cpbsvx('N', 'U', 1, 1, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = '/';
         cpbsvx('F', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ = 'Y';
         cpbsvx('F', 'U', 1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPBSVX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'PT' ) ) {

         // CPTSV

        srnamc.SRNAMT = 'CPTSV ';
         INFOT = 1;
         cptsv(-1, 0, R, A( 1, 1 ), B, 1, INFO );
         chkxer('CPTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cptsv(0, -1, R, A( 1, 1 ), B, 1, INFO );
         chkxer('CPTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cptsv(2, 0, R, A( 1, 1 ), B, 1, INFO );
         chkxer('CPTSV ', INFOT, NOUT, LERR, OK );

         // CPTSVX

        srnamc.SRNAMT = 'CPTSVX';
         INFOT = 1;
         cptsvx('/', 0, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cptsvx('N', -1, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cptsvx('N', 0, -1, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cptsvx('N', 2, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cptsvx('N', 2, 0, R, A( 1, 1 ), RF, AF( 1, 1 ), B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CPTSVX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'HE' ) ) {

         // CHESV

        srnamc.SRNAMT = 'CHESV ';
         INFOT = 1;
         chesv('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chesv('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chesv('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chesv('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('CHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chesv('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('CHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chesv('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('CHESV ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chesv('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('CHESV ', INFOT, NOUT, LERR, OK );

         // CHESVX

        srnamc.SRNAMT = 'CHESVX';
         INFOT = 1;
         chesvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chesvx('N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chesvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chesvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         chesvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chesvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         chesvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, RW, INFO );
         chkxer('CHESVX', INFOT, NOUT, LERR, OK );

         // CHESVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
        srnamc.SRNAMT = 'CHESVXX';
         INFOT = 1;
         chesvxx('/', 'U', 0, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chesvxx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chesvxx('N', 'U', -1, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chesvxx('N', 'U', 0, -1, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         chesvxx('N', 'U', 2, 0, A, 1, AF, 2, IP, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chesvxx('N', 'U', 2, 0, A, 2, AF, 1, IP, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         EQ = '/';
         chesvxx('F', 'U', 0, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = 'Y';
         chesvxx('F', 'U', 1, 0, A, 1, AF, 1, IP, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         chesvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         chesvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CHESVXX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'HR' ) ) {

         // CHESV_ROOK

        srnamc.SRNAMT = 'CHESV_ROOK';
         INFOT = 1;
         chesv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chesv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chesv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chesv_rook('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chesv_rook('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chesv_rook('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'HK' ) ) {

         // CHESV_RK

         // Test error exits of the driver that uses factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

        srnamc.SRNAMT = 'CHESV_RK';
         INFOT = 1;
         chesv_rk('/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chesv_rk('U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chesv_rk('U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chesv_rk('U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO );
         chkxer('CHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chesv_rk('U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chesv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO );
         chkxer('CHESV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chesv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO );
         chkxer('CHESV_RK', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'HP' ) ) {

         // CHPSV

        srnamc.SRNAMT = 'CHPSV ';
         INFOT = 1;
         chpsv('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('CHPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chpsv('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('CHPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chpsv('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('CHPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chpsv('U', 2, 0, A, IP, B, 1, INFO );
         chkxer('CHPSV ', INFOT, NOUT, LERR, OK );

         // CHPSVX

        srnamc.SRNAMT = 'CHPSVX';
         INFOT = 1;
         chpsvx('/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chpsvx('N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chpsvx('N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chpsvx('N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chpsvx('N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CHPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chpsvx('N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CHPSVX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'SY' ) ) {

         // CSYSV

        srnamc.SRNAMT = 'CSYSV ';
         INFOT = 1;
         csysv('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysv('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysv('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csysv('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csysv('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csysv('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );

         // CSYSVX

        srnamc.SRNAMT = 'CSYSVX';
         INFOT = 1;
         csysvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysvx('N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csysvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         csysvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csysvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         csysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         csysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         csysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, RW, INFO );
         chkxer('CSYSVX', INFOT, NOUT, LERR, OK );

         // CSYSVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
        srnamc.SRNAMT = 'CSYSVXX';
         INFOT = 1;
         EQ = 'N';
         csysvxx('/', 'U', 0, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysvxx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C,  NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysvxx('N', 'U', -1, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         EQ = '/';
         csysvxx('N', 'U', 0, -1, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         EQ = 'Y';
         INFOT = 6;
         csysvxx('N', 'U', 2, 0, A, 1, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csysvxx('N', 'U', 2, 0, A, 2, AF, 1, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, 'A', R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ='Y';
         csysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ='Y';
         R[1] = -ONE;
         csysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'N';
         csysvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         csysvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, RW, INFO );
         chkxer('CSYSVXX', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'SR' ) ) {

         // CSYSV_ROOK

        srnamc.SRNAMT = 'CSYSV_ROOK';
         INFOT = 1;
         csysv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csysv_rook('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csysv_rook('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csysv_rook('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'SK' ) ) {

         // CSYSV_RK

         // Test error exits of the driver that uses factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

        srnamc.SRNAMT = 'CSYSV_RK';
         INFOT = 1;
         csysv_rk('/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysv_rk('U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysv_rk('U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csysv_rk('U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO );
         chkxer('CSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         csysv_rk('U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         csysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO );
         chkxer('CSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         csysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO );
         chkxer('CSYSV_RK', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'SP' ) ) {

         // CSPSV

        srnamc.SRNAMT = 'CSPSV ';
         INFOT = 1;
         cspsv('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('CSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cspsv('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('CSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cspsv('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('CSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cspsv('U', 2, 0, A, IP, B, 1, INFO );
         chkxer('CSPSV ', INFOT, NOUT, LERR, OK );

         // CSPSVX

        srnamc.SRNAMT = 'CSPSVX';
         INFOT = 1;
         cspsvx('/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cspsvx('N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cspsvx('N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cspsvx('N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cspsvx('N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, RW, INFO );
         chkxer('CSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cspsvx('N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, RW, INFO );
         chkxer('CSPSVX', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH;
      } else {
         WRITE( NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT( 1X, A3, ' drivers passed the tests of the error exits' );
 9998 FORMAT( ' *** ', A3, ' drivers failed the tests of the error ', 'exits ***' );

      return;
      }
