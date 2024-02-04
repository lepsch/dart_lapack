      void derrvx(PATH, NUNIT ) {

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
      double             RCOND, RPVGRW, BERR;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX ), IW( NMAX );
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), C( NMAX ), E( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), W( 2*NMAX ), X( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DGBSV, DGBSVX, DGESV, DGESVX, DGTSV, DGTSVX, DPBSV, DPBSVX, DPOSV, DPOSVX, DPPSV, DPPSVX, DPTSV, DPTSVX, DSPSV, DSPSVX, DSYSV, DSYSV_RK, DSYSV_ROOK, DSYSVX, DGESVXX, DSYSVXX, DPOSVXX, DGBSVXX
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
            A[I, J] = 1.0 / (I+J).toDouble();
            AF[I, J] = 1.0 / (I+J).toDouble();
         } // 10
         B[J] = 0.0;
         E[J] = 0.0;
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

      if ( LSAMEN( 2, C2, 'GE' ) ) {

         // DGESV

         SRNAMT = 'DGESV ';
         INFOT = 1;
         dgesv(-1, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgesv(0, -1, A, 1, IP, B, 1, INFO );
         chkxer('DGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgesv(2, 1, A, 1, IP, B, 2, INFO );
         chkxer('DGESV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgesv(2, 1, A, 2, IP, B, 1, INFO );
         chkxer('DGESV ', INFOT, NOUT, LERR, OK );

         // DGESVX

         SRNAMT = 'DGESVX';
         INFOT = 1;
         dgesvx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgesvx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgesvx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgesvx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dgesvx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgesvx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = '/';
         dgesvx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ = 'R';
         dgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         EQ = 'C';
         dgesvx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dgesvx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGESVX', INFOT, NOUT, LERR, OK );

         // DGESVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
         SRNAMT = 'DGESVXX';
         INFOT = 1;
         dgesvxx('/', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgesvxx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgesvxx('N', 'N', -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgesvxx('N', 'N', 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dgesvxx('N', 'N', 2, 1, A, 1, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgesvxx('N', 'N', 2, 1, A, 2, AF, 1, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = '/';
         dgesvxx('F', 'N', 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ = 'R';
         dgesvxx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         EQ = 'C';
         dgesvxx('F', 'N', 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dgesvxx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dgesvxx('N', 'N', 2, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGESVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // DGBSV

         SRNAMT = 'DGBSV ';
         INFOT = 1;
         dgbsv(-1, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgbsv(1, -1, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgbsv(1, 0, -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgbsv(0, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dgbsv(1, 1, 1, 0, A, 3, IP, B, 1, INFO );
         chkxer('DGBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dgbsv(2, 0, 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('DGBSV ', INFOT, NOUT, LERR, OK );

         // DGBSVX

         SRNAMT = 'DGBSVX';
         INFOT = 1;
         dgbsvx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgbsvx('N', '/', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgbsvx('N', 'N', -1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgbsvx('N', 'N', 1, -1, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgbsvx('N', 'N', 1, 0, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dgbsvx('N', 'N', 0, 0, 0, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgbsvx('N', 'N', 1, 1, 1, 0, A, 2, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dgbsvx('N', 'N', 1, 1, 1, 0, A, 3, AF, 3, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         EQ = '/';
         dgbsvx('F', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'R';
         dgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         EQ = 'C';
         dgbsvx('F', 'N', 1, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dgbsvx('N', 'N', 2, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGBSVX', INFOT, NOUT, LERR, OK );

         // DGBSVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
         SRNAMT = 'DGBSVXX';
         INFOT = 1;
         dgbsvxx('/', 'N', 0, 0, 0, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgbsvxx('N', '/', 0, 1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgbsvxx('N', 'N', -1, 1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgbsvxx('N', 'N', 2, -1, 1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgbsvxx('N', 'N', 2, 1, -1, 0, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dgbsvxx('N', 'N', 0, 1, 1, -1, A, 1, AF, 1, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgbsvxx('N', 'N', 2, 1, 1, 1, A, 2, AF, 2, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 3, IP, EQ, R, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         EQ = '/';
         dgbsvxx('F', 'N', 0, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'R';
         dgbsvxx('F', 'N', 1, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         EQ = 'C';
         dgbsvxx('F', 'N', 1, 1, 1, 0, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 4, IP, EQ, R, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dgbsvxx('N', 'N', 2, 1, 1, 1, A, 3, AF, 4, IP, EQ, R, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DGBSVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'GT' ) ) {

         // DGTSV

         SRNAMT = 'DGTSV ';
         INFOT = 1;
         dgtsv(-1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('DGTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgtsv(0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('DGTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgtsv(2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), B, 1, INFO );
         chkxer('DGTSV ', INFOT, NOUT, LERR, OK );

         // DGTSVX

         SRNAMT = 'DGTSVX';
         INFOT = 1;
         dgtsvx('/', 'N', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgtsvx('N', '/', 0, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgtsvx('N', 'N', -1, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgtsvx('N', 'N', 0, -1, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dgtsvx('N', 'N', 2, 0, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), AF( 1, 1 ), AF( 1, 2 ), AF( 1, 3 ), AF( 1, 4 ), IP, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DGTSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PO' ) ) {

         // DPOSV

         SRNAMT = 'DPOSV ';
         INFOT = 1;
         dposv('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('DPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dposv('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('DPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dposv('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('DPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dposv('U', 2, 0, A, 1, B, 2, INFO );
         chkxer('DPOSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dposv('U', 2, 0, A, 2, B, 1, INFO );
         chkxer('DPOSV ', INFOT, NOUT, LERR, OK );

         // DPOSVX

         SRNAMT = 'DPOSVX';
         INFOT = 1;
         dposvx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dposvx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dposvx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dposvx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dposvx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dposvx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         EQ = '/';
         dposvx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = 'Y';
         dposvx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dposvx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPOSVX', INFOT, NOUT, LERR, OK );

         // DPOSVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
         SRNAMT = 'DPOSVXX';
         INFOT = 1;
         dposvxx('/', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dposvxx('N', '/', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dposvxx('N', 'U', -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dposvxx('N', 'U', 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dposvxx('N', 'U', 2, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dposvxx('N', 'U', 2, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         EQ = '/';
         dposvxx('F', 'U', 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = 'Y';
         dposvxx('F', 'U', 1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dposvxx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dposvxx('N', 'U', 2, 0, A, 2, AF, 2, EQ, C, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DPOSVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // DPPSV

         SRNAMT = 'DPPSV ';
         INFOT = 1;
         dppsv('/', 0, 0, A, B, 1, INFO );
         chkxer('DPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dppsv('U', -1, 0, A, B, 1, INFO );
         chkxer('DPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dppsv('U', 0, -1, A, B, 1, INFO );
         chkxer('DPPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dppsv('U', 2, 0, A, B, 1, INFO );
         chkxer('DPPSV ', INFOT, NOUT, LERR, OK );

         // DPPSVX

         SRNAMT = 'DPPSVX';
         INFOT = 1;
         dppsvx('/', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dppsvx('N', '/', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dppsvx('N', 'U', -1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dppsvx('N', 'U', 0, -1, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         EQ = '/';
         dppsvx('F', 'U', 0, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         EQ = 'Y';
         dppsvx('F', 'U', 1, 0, A, AF, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dppsvx('N', 'U', 2, 0, A, AF, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPPSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // DPBSV

         SRNAMT = 'DPBSV ';
         INFOT = 1;
         dpbsv('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('DPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpbsv('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('DPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpbsv('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('DPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dpbsv('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('DPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dpbsv('U', 1, 1, 0, A, 1, B, 2, INFO );
         chkxer('DPBSV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dpbsv('U', 2, 0, 0, A, 1, B, 1, INFO );
         chkxer('DPBSV ', INFOT, NOUT, LERR, OK );

         // DPBSVX

         SRNAMT = 'DPBSVX';
         INFOT = 1;
         dpbsvx('/', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpbsvx('N', '/', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dpbsvx('N', 'U', -1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dpbsvx('N', 'U', 1, -1, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dpbsvx('N', 'U', 0, 0, -1, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dpbsvx('N', 'U', 1, 1, 0, A, 1, AF, 2, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dpbsvx('N', 'U', 1, 1, 0, A, 2, AF, 1, EQ, C, B, 2, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         EQ = '/';
         dpbsvx('F', 'U', 0, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ = 'Y';
         dpbsvx('F', 'U', 1, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dpbsvx('N', 'U', 2, 0, 0, A, 1, AF, 1, EQ, C, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DPBSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // DPTSV

         SRNAMT = 'DPTSV ';
         INFOT = 1;
         dptsv(-1, 0, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('DPTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dptsv(0, -1, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('DPTSV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dptsv(2, 0, A( 1, 1 ), A( 1, 2 ), B, 1, INFO );
         chkxer('DPTSV ', INFOT, NOUT, LERR, OK );

         // DPTSVX

         SRNAMT = 'DPTSVX';
         INFOT = 1;
         dptsvx('/', 0, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dptsvx('N', -1, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dptsvx('N', 0, -1, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dptsvx('N', 2, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 1, X, 2, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dptsvx('N', 2, 0, A( 1, 1 ), A( 1, 2 ), AF( 1, 1 ), AF( 1, 2 ), B, 2, X, 1, RCOND, R1, R2, W, INFO );
         chkxer('DPTSVX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SY' ) ) {

         // DSYSV

         SRNAMT = 'DSYSV ';
         INFOT = 1;
         dsysv('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsysv('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsysv('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsysv('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsysv('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsysv('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('DSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsysv('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('DSYSV ', INFOT, NOUT, LERR, OK );

         // DSYSVX

         SRNAMT = 'DSYSVX';
         INFOT = 1;
         dsysvx('/', 'U', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsysvx('N', '/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsysvx('N', 'U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsysvx('N', 'U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, RCOND, R1, R2, W, 1, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsysvx('N', 'U', 2, 0, A, 1, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsysvx('N', 'U', 2, 0, A, 2, AF, 1, IP, B, 2, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 1, X, 2, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 1, RCOND, R1, R2, W, 4, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dsysvx('N', 'U', 2, 0, A, 2, AF, 2, IP, B, 2, X, 2, RCOND, R1, R2, W, 3, IW, INFO );
         chkxer('DSYSVX', INFOT, NOUT, LERR, OK );

         // DSYSVXX

         N_ERR_BNDS = 3;
         NPARAMS = 1;
         SRNAMT = 'DSYSVXX';
         INFOT = 1;
         EQ = 'N';
         dsysvxx('/', 'U', 0, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C,  NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsysvxx('N', '/', 0, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsysvxx('N', 'U', -1, 0, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C,  NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         EQ = '/';
         dsysvxx('N', 'U', 0, -1, A, 1, AF, 1, IP, EQ, R, B, 1, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         EQ = 'Y';
         INFOT = 6;
         dsysvxx('N', 'U', 2, 0, A, 1, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsysvxx('N', 'U', 2, 0, A, 2, AF, 1, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, 'A', R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ='Y';
         dsysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         EQ='Y';
         R[1] = -ONE;
         dsysvxx('F', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         EQ = 'N';
         dsysvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 1, X, 2, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dsysvxx('N', 'U', 2, 0, A, 2, AF, 2, IP, EQ, R, B, 2, X, 1, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYSVXX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // DSYSV_ROOK

         SRNAMT = 'DSYSV_ROOK';
         INFOT = 1;
         dsysv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsysv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsysv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsysv_rook('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsysv_rook('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsysv_rook('U', 0, 0, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('DSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsysv_rook('U', 0, 0, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('DSYSV_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // DSYSV_RK

         // Test error exits of the driver that uses factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         SRNAMT = 'DSYSV_RK';
         INFOT = 1;
         dsysv_rk('/', 0, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsysv_rk('U', -1, 0, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsysv_rk('U', 0, -1, A, 1, E, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsysv_rk('U', 2, 0, A, 1, E, IP, B, 2, W, 1, INFO );
         chkxer('DSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsysv_rk('U', 2, 0, A, 2, E, IP, B, 1, W, 1, INFO );
         chkxer('DSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, 0, INFO );
         chkxer('DSYSV_RK', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsysv_rk('U', 0, 0, A, 1, E, IP, B, 1, W, -2, INFO );
         chkxer('DSYSV_RK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // DSPSV

         SRNAMT = 'DSPSV ';
         INFOT = 1;
         dspsv('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('DSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dspsv('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('DSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dspsv('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('DSPSV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dspsv('U', 2, 0, A, IP, B, 1, INFO );
         chkxer('DSPSV ', INFOT, NOUT, LERR, OK );

         // DSPSVX

         SRNAMT = 'DSPSVX';
         INFOT = 1;
         dspsvx('/', 'U', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dspsvx('N', '/', 0, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dspsvx('N', 'U', -1, 0, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dspsvx('N', 'U', 0, -1, A, AF, IP, B, 1, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dspsvx('N', 'U', 2, 0, A, AF, IP, B, 1, X, 2, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dspsvx('N', 'U', 2, 0, A, AF, IP, B, 2, X, 1, RCOND, R1, R2, W, IW, INFO );
         chkxer('DSPSVX', INFOT, NOUT, LERR, OK );
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
