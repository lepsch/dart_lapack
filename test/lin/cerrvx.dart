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
      // ..
      // .. Local Scalars ..
      String             EQ;
      String             C2;
      int                I, INFO, J;
      double               RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double               C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), RF( NMAX ), RW( NMAX )       Complex            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGBSV, CGBSVX, CGESV, CGESVX, CGTSV, CGTSVX, CHESV, CHESV_RK, CHESV_ROOK, CHESVX, CHKXER, CHPSV, CHPSVX, CPBSV, CPBSVX, CPOSV, CPOSVX, CPPSV, CPPSVX, CPTSV, CPTSVX, CSPSV, CSPSVX, CSYSV, CSYSV_AA, CSYSV_RK, CSYSV_ROOK, CSYSVX, CSYSV_AA_2STAGE
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

         // CGESV

         SRNAMT = 'CGESV ';
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

         SRNAMT = 'CGESVX';
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

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // CGBSV

         SRNAMT = 'CGBSV ';
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

         SRNAMT = 'CGBSVX';
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

      } else if ( LSAMEN( 2, C2, 'GT' ) ) {

         // CGTSV

         SRNAMT = 'CGTSV ';
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

         SRNAMT = 'CGTSVX';
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

      } else if ( LSAMEN( 2, C2, 'PO' ) ) {

         // CPOSV

         SRNAMT = 'CPOSV ';
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

         SRNAMT = 'CPOSVX';
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

      } else if ( LSAMEN( 2, C2, 'PP' ) ) {

         // CPPSV

         SRNAMT = 'CPPSV ';
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

         SRNAMT = 'CPPSVX';
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

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // CPBSV

         SRNAMT = 'CPBSV ';
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

         SRNAMT = 'CPBSVX';
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

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // CPTSV

         SRNAMT = 'CPTSV ';
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

         SRNAMT = 'CPTSVX';
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

      } else if ( LSAMEN( 2, C2, 'HE' ) ) {

         // CHESV

         SRNAMT = 'CHESV ';
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

         SRNAMT = 'CHESVX';
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

      } else if ( LSAMEN( 2, C2, 'HR' ) ) {

         // CHESV_ROOK

         SRNAMT = 'CHESV_ROOK';
         INFOT = 1;
         chesv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chesv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chesv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chesv_rook('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
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

      } else if ( LSAMEN( 2, C2, 'HK' ) ) {

         // CHESV_RK

         // Test error exits of the driver that uses factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         SRNAMT = 'CHESV_RK';
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

      } else if ( LSAMEN( 2, C2, 'HA' ) ) {

         // CHESV_AASEN

         SRNAMT = 'CHESV_AA';
         INFOT = 1;
         chesv_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chesv_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chesv_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chesv_aa('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('CHESV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chesv_aa('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chesv_aa('U', 3, 1, A, 3, IP, B, 3, W, 6, INFO );
         chkxer('CHESV_AA', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'H2' ) ) {

         // CHESV_AASEN_2STAGE

         SRNAMT = 'CHESV_AA_2STAGE';
         INFOT = 1;
         chesv_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chesv_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chesv_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chesv_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chesv_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 2, W, 1, INFO );
         chkxer('CHESV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chesv_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, W, 1, INFO );
         chkxer('CHESV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chesv_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 2, W, 1, INFO );
         chkxer('CHESV_AA_2STAGE', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SA' ) ) {

         // CSYSV_AASEN

         SRNAMT = 'CSYSV_AA';
         INFOT = 1;
         csysv_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysv_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysv_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csysv_aa('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('CSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csysv_aa('U', 2, 0, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csysv_aa('U', 3, 1, A, 3, IP, B, 3, W, 6, INFO );
         chkxer('CSYSV_AA', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'S2' ) ) {

         // CSYSV_AASEN_2STAGE

         SRNAMT = 'CSYSV_AA_2STAGE';
         INFOT = 1;
         csysv_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysv_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysv_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csysv_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csysv_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 2, W, 1, INFO );
         chkxer('CSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         csysv_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         csysv_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 2, W, 1, INFO );
         chkxer('CSYSV_AA_2STAGE', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HP' ) ) {

         // CHPSV

         SRNAMT = 'CHPSV ';
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

         SRNAMT = 'CHPSVX';
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

      } else if ( LSAMEN( 2, C2, 'SY' ) ) {

         // CSYSV

         SRNAMT = 'CSYSV ';
         INFOT = 1;
         csysv('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysv('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysv('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csysv('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
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

         SRNAMT = 'CSYSVX';
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

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // CSYSV_ROOK

         SRNAMT = 'CSYSV_ROOK';
         INFOT = 1;
         csysv_rook('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csysv_rook('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csysv_rook('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYSV_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csysv_rook('U', 2, 0, A, 1, IP, B, 2, W, 1, INFO );
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

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // CSYSV_RK

         // Test error exits of the driver that uses factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         SRNAMT = 'CSYSV_RK';
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

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // CSPSV

         SRNAMT = 'CSPSV ';
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

         SRNAMT = 'CSPSVX';
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
