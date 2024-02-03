      SUBROUTINE CERRGE( PATH, NUNIT );

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
      String             C2;
      int                I, INFO, J;
      REAL               ANRM, CCOND, RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      REAL               R( NMAX ), R1( NMAX ), R2( NMAX );
      COMPLEX            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CGBCON, CGBEQU, CGBRFS, CGBTF2, CGBTRF, CGBTRS, CGECON, CGEEQU, CGERFS, CGETF2, CGETRF, CGETRI, CGETRS, CHKXER
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
            A( I, J ) = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF( I, J ) = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
         } // 10
         B( J ) = 0.;
         R1( J ) = 0.;
         R2( J ) = 0.;
         W( J ) = 0.;
         X( J ) = 0.;
         IP( J ) = J;
      } // 20
      OK = true;

      // Test error exits of the routines that use the LU decomposition
      // of a general matrix.

      if ( LSAMEN( 2, C2, 'GE' ) ) {

         // CGETRF

         SRNAMT = 'CGETRF';
         INFOT = 1;
         cgetrf(-1, 0, A, 1, IP, INFO );
         chkxer('CGETRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgetrf(0, -1, A, 1, IP, INFO );
         chkxer('CGETRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgetrf(2, 1, A, 1, IP, INFO );
         chkxer('CGETRF', INFOT, NOUT, LERR, OK );

         // CGETF2

         SRNAMT = 'CGETF2';
         INFOT = 1;
         cgetf2(-1, 0, A, 1, IP, INFO );
         chkxer('CGETF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgetf2(0, -1, A, 1, IP, INFO );
         chkxer('CGETF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgetf2(2, 1, A, 1, IP, INFO );
         chkxer('CGETF2', INFOT, NOUT, LERR, OK );

         // CGETRI

         SRNAMT = 'CGETRI';
         INFOT = 1;
         cgetri(-1, A, 1, IP, W, 1, INFO );
         chkxer('CGETRI', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgetri(2, A, 1, IP, W, 2, INFO );
         chkxer('CGETRI', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgetri(2, A, 2, IP, W, 1, INFO );
         chkxer('CGETRI', INFOT, NOUT, LERR, OK );

         // CGETRS

         SRNAMT = 'CGETRS';
         INFOT = 1;
         cgetrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('CGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgetrs('N', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('CGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgetrs('N', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('CGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgetrs('N', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('CGETRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgetrs('N', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('CGETRS', INFOT, NOUT, LERR, OK );

         // CGERFS

         SRNAMT = 'CGERFS';
         INFOT = 1;
         cgerfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgerfs('N', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgerfs('N', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgerfs('N', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgerfs('N', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgerfs('N', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CGERFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgerfs('N', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CGERFS', INFOT, NOUT, LERR, OK );

         // CGECON

         SRNAMT = 'CGECON';
         INFOT = 1;
         cgecon('/', 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CGECON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgecon('1', -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CGECON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgecon('1', 2, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('CGECON', INFOT, NOUT, LERR, OK );

         // CGEEQU

         SRNAMT = 'CGEEQU';
         INFOT = 1;
         cgeequ(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('CGEEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgeequ(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('CGEEQU', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgeequ(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('CGEEQU', INFOT, NOUT, LERR, OK );

      // Test error exits of the routines that use the LU decomposition
      // of a general band matrix.

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // CGBTRF

         SRNAMT = 'CGBTRF';
         INFOT = 1;
         cgbtrf(-1, 0, 0, 0, A, 1, IP, INFO );
         chkxer('CGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbtrf(0, -1, 0, 0, A, 1, IP, INFO );
         chkxer('CGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbtrf(1, 1, -1, 0, A, 1, IP, INFO );
         chkxer('CGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbtrf(1, 1, 0, -1, A, 1, IP, INFO );
         chkxer('CGBTRF', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgbtrf(2, 2, 1, 1, A, 3, IP, INFO );
         chkxer('CGBTRF', INFOT, NOUT, LERR, OK );

         // CGBTF2

         SRNAMT = 'CGBTF2';
         INFOT = 1;
         cgbtf2(-1, 0, 0, 0, A, 1, IP, INFO );
         chkxer('CGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbtf2(0, -1, 0, 0, A, 1, IP, INFO );
         chkxer('CGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbtf2(1, 1, -1, 0, A, 1, IP, INFO );
         chkxer('CGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbtf2(1, 1, 0, -1, A, 1, IP, INFO );
         chkxer('CGBTF2', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgbtf2(2, 2, 1, 1, A, 3, IP, INFO );
         chkxer('CGBTF2', INFOT, NOUT, LERR, OK );

         // CGBTRS

         SRNAMT = 'CGBTRS';
         INFOT = 1;
         cgbtrs('/', 0, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('CGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbtrs('N', -1, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('CGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbtrs('N', 1, -1, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('CGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbtrs('N', 1, 0, -1, 1, A, 1, IP, B, 1, INFO );
         chkxer('CGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgbtrs('N', 1, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('CGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgbtrs('N', 2, 1, 1, 1, A, 3, IP, B, 2, INFO );
         chkxer('CGBTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgbtrs('N', 2, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('CGBTRS', INFOT, NOUT, LERR, OK );

         // CGBRFS

         SRNAMT = 'CGBRFS';
         INFOT = 1;
         cgbrfs('/', 0, 0, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbrfs('N', -1, 0, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbrfs('N', 1, -1, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbrfs('N', 1, 0, -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgbrfs('N', 1, 0, 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgbrfs('N', 2, 1, 1, 1, A, 2, AF, 4, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgbrfs('N', 2, 1, 1, 1, A, 3, AF, 3, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CGBRFS', INFOT, NOUT, LERR, OK );

         // CGBCON

         SRNAMT = 'CGBCON';
         INFOT = 1;
         cgbcon('/', 0, 0, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('CGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbcon('1', -1, 0, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('CGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbcon('1', 1, -1, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('CGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbcon('1', 1, 0, -1, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('CGBCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgbcon('1', 2, 1, 1, A, 3, IP, ANRM, RCOND, W, R, INFO );
         chkxer('CGBCON', INFOT, NOUT, LERR, OK );

         // CGBEQU

         SRNAMT = 'CGBEQU';
         INFOT = 1;
         cgbequ(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('CGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgbequ(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('CGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgbequ(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('CGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgbequ(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('CGBEQU', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgbequ(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('CGBEQU', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN;

      // End of CERRGE

      }
