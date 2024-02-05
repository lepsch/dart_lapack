      void zerrge(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                infoc.NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, INFO, J;
      double             ANRM, CCOND, RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             R( NMAX ), R1( NMAX ), R2( NMAX );
      Complex         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGBCON, ZGBEQU, ZGBRFS, ZGBTF2, ZGBTRF, ZGBTRS, ZGECON, ZGEEQU, ZGERFS, ZGETF2, ZGETRF, ZGETRI, ZGETRS
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / (I+J).toDouble() );
         } // 10
         B[J] = 0.0;
         R1[J] = 0.0;
         R2[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
         IP[J] = J;
      } // 20
      infoc.OK = true;

      // Test error exits of the routines that use the LU decomposition
      // of a general matrix.

      if ( lsamen( 2, C2, 'GE' ) ) {

         // ZGETRF

        srnamc.SRNAMT = 'ZGETRF';
         infoc.INFOT = 1;
         zgetrf(-1, 0, A, 1, IP, INFO );
         chkxer('ZGETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgetrf(0, -1, A, 1, IP, INFO );
         chkxer('ZGETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgetrf(2, 1, A, 1, IP, INFO );
         chkxer('ZGETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGETF2

        srnamc.SRNAMT = 'ZGETF2';
         infoc.INFOT = 1;
         zgetf2(-1, 0, A, 1, IP, INFO );
         chkxer('ZGETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgetf2(0, -1, A, 1, IP, INFO );
         chkxer('ZGETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgetf2(2, 1, A, 1, IP, INFO );
         chkxer('ZGETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGETRI

        srnamc.SRNAMT = 'ZGETRI';
         infoc.INFOT = 1;
         zgetri(-1, A, 1, IP, W, 1, INFO );
         chkxer('ZGETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgetri(2, A, 1, IP, W, 2, INFO );
         chkxer('ZGETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zgetri(2, A, 2, IP, W, 1, INFO );
         chkxer('ZGETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGETRS

        srnamc.SRNAMT = 'ZGETRS';
         infoc.INFOT = 1;
         zgetrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgetrs('N', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgetrs('N', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zgetrs('N', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgetrs('N', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZGETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGERFS

        srnamc.SRNAMT = 'ZGERFS';
         infoc.INFOT = 1;
         zgerfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgerfs('N', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgerfs('N', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zgerfs('N', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zgerfs('N', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zgerfs('N', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         zgerfs('N', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGECON

        srnamc.SRNAMT = 'ZGECON';
         infoc.INFOT = 1;
         zgecon('/', 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZGECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgecon('1', -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZGECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgecon('1', 2, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZGECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGEEQU

        srnamc.SRNAMT = 'ZGEEQU';
         infoc.INFOT = 1;
         zgeequ(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgeequ(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgeequ(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGEEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Test error exits of the routines that use the LU decomposition
      // of a general band matrix.

      } else if ( lsamen( 2, C2, 'GB' ) ) {

         // ZGBTRF

        srnamc.SRNAMT = 'ZGBTRF';
         infoc.INFOT = 1;
         zgbtrf(-1, 0, 0, 0, A, 1, IP, INFO );
         chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgbtrf(0, -1, 0, 0, A, 1, IP, INFO );
         chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgbtrf(1, 1, -1, 0, A, 1, IP, INFO );
         chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgbtrf(1, 1, 0, -1, A, 1, IP, INFO );
         chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zgbtrf(2, 2, 1, 1, A, 3, IP, INFO );
         chkxer('ZGBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGBTF2

        srnamc.SRNAMT = 'ZGBTF2';
         infoc.INFOT = 1;
         zgbtf2(-1, 0, 0, 0, A, 1, IP, INFO );
         chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgbtf2(0, -1, 0, 0, A, 1, IP, INFO );
         chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgbtf2(1, 1, -1, 0, A, 1, IP, INFO );
         chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgbtf2(1, 1, 0, -1, A, 1, IP, INFO );
         chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zgbtf2(2, 2, 1, 1, A, 3, IP, INFO );
         chkxer('ZGBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGBTRS

        srnamc.SRNAMT = 'ZGBTRS';
         infoc.INFOT = 1;
         zgbtrs('/', 0, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgbtrs('N', -1, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgbtrs('N', 1, -1, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgbtrs('N', 1, 0, -1, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zgbtrs('N', 1, 0, 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zgbtrs('N', 2, 1, 1, 1, A, 3, IP, B, 2, INFO );
         chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zgbtrs('N', 2, 0, 0, 1, A, 1, IP, B, 1, INFO );
         chkxer('ZGBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGBRFS

        srnamc.SRNAMT = 'ZGBRFS';
         infoc.INFOT = 1;
         zgbrfs('/', 0, 0, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgbrfs('N', -1, 0, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgbrfs('N', 1, -1, 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgbrfs('N', 1, 0, -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zgbrfs('N', 1, 0, 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zgbrfs('N', 2, 1, 1, 1, A, 2, AF, 4, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         zgbrfs('N', 2, 1, 1, 1, A, 3, AF, 3, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         zgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 14;
         zgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZGBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGBCON

        srnamc.SRNAMT = 'ZGBCON';
         infoc.INFOT = 1;
         zgbcon('/', 0, 0, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgbcon('1', -1, 0, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgbcon('1', 1, -1, 0, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgbcon('1', 1, 0, -1, A, 1, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zgbcon('1', 2, 1, 1, A, 3, IP, ANRM, RCOND, W, R, INFO );
         chkxer('ZGBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGBEQU

        srnamc.SRNAMT = 'ZGBEQU';
         infoc.INFOT = 1;
         zgbequ(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgbequ(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgbequ(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgbequ(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zgbequ(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO );
         chkxer('ZGBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
