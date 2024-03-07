      void zerrpo(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      String             C2;
      int                I, INFO, J;
      double             ANRM, RCOND;
      double             R( NMAX ), R1( NMAX ), R2( NMAX );
      Complex         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZPBCON, ZPBEQU, ZPBRFS, ZPBTF2, ZPBTRF, ZPBTRS, ZPOCON, ZPOEQU, ZPORFS, ZPOTF2, ZPOTRF, ZPOTRI, ZPOTRS, ZPPCON, ZPPEQU, ZPPRFS, ZPPTRF, ZPPTRI, ZPPTRS
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

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH.substring( 1, 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / (I+J).toDouble() );
         } // 10
         B[J] = 0.0;
         R1[J] = 0.0;
         R2[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
      } // 20
      ANRM = 1.0;
      infoc.OK.value = true;

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite matrix.

      if ( lsamen( 2, C2, 'PO' ) ) {

         // ZPOTRF

        srnamc.SRNAMT = 'ZPOTRF';
         infoc.INFOT = 1;
         zpotrf('/', 0, A, 1, INFO );
         chkxer('ZPOTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpotrf('U', -1, A, 1, INFO );
         chkxer('ZPOTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zpotrf('U', 2, A, 1, INFO );
         chkxer('ZPOTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPOTF2

        srnamc.SRNAMT = 'ZPOTF2';
         infoc.INFOT = 1;
         zpotf2('/', 0, A, 1, INFO );
         chkxer('ZPOTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpotf2('U', -1, A, 1, INFO );
         chkxer('ZPOTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zpotf2('U', 2, A, 1, INFO );
         chkxer('ZPOTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPOTRI

        srnamc.SRNAMT = 'ZPOTRI';
         infoc.INFOT = 1;
         zpotri('/', 0, A, 1, INFO );
         chkxer('ZPOTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpotri('U', -1, A, 1, INFO );
         chkxer('ZPOTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zpotri('U', 2, A, 1, INFO );
         chkxer('ZPOTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPOTRS

        srnamc.SRNAMT = 'ZPOTRS';
         infoc.INFOT = 1;
         zpotrs('/', 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPOTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpotrs('U', -1, 0, A, 1, B, 1, INFO );
         chkxer('ZPOTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpotrs('U', 0, -1, A, 1, B, 1, INFO );
         chkxer('ZPOTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zpotrs('U', 2, 1, A, 1, B, 2, INFO );
         chkxer('ZPOTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zpotrs('U', 2, 1, A, 2, B, 1, INFO );
         chkxer('ZPOTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPORFS

        srnamc.SRNAMT = 'ZPORFS';
         infoc.INFOT = 1;
         zporfs('/', 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPORFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zporfs('U', -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPORFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zporfs('U', 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPORFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zporfs('U', 2, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPORFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zporfs('U', 2, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPORFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         zporfs('U', 2, 1, A, 2, AF, 2, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPORFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         zporfs('U', 2, 1, A, 2, AF, 2, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPORFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPOCON

        srnamc.SRNAMT = 'ZPOCON';
         infoc.INFOT = 1;
         zpocon('/', 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPOCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpocon('U', -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPOCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zpocon('U', 2, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPOCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zpocon('U', 1, A, 1, -ANRM, RCOND, W, R, INFO );
         chkxer('ZPOCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPOEQU

        srnamc.SRNAMT = 'ZPOEQU';
         infoc.INFOT = 1;
         zpoequ(-1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPOEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpoequ(2, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPOEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite packed matrix.

      } else if ( lsamen( 2, C2, 'PP' ) ) {

         // ZPPTRF

        srnamc.SRNAMT = 'ZPPTRF';
         infoc.INFOT = 1;
         zpptrf('/', 0, A, INFO );
         chkxer('ZPPTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpptrf('U', -1, A, INFO );
         chkxer('ZPPTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPPTRI

        srnamc.SRNAMT = 'ZPPTRI';
         infoc.INFOT = 1;
         zpptri('/', 0, A, INFO );
         chkxer('ZPPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpptri('U', -1, A, INFO );
         chkxer('ZPPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPPTRS

        srnamc.SRNAMT = 'ZPPTRS';
         infoc.INFOT = 1;
         zpptrs('/', 0, 0, A, B, 1, INFO );
         chkxer('ZPPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpptrs('U', -1, 0, A, B, 1, INFO );
         chkxer('ZPPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpptrs('U', 0, -1, A, B, 1, INFO );
         chkxer('ZPPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zpptrs('U', 2, 1, A, B, 1, INFO );
         chkxer('ZPPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPPRFS

        srnamc.SRNAMT = 'ZPPRFS';
         infoc.INFOT = 1;
         zpprfs('/', 0, 0, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpprfs('U', -1, 0, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpprfs('U', 0, -1, A, AF, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zpprfs('U', 2, 1, A, AF, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         zpprfs('U', 2, 1, A, AF, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPPCON

        srnamc.SRNAMT = 'ZPPCON';
         infoc.INFOT = 1;
         zppcon('/', 0, A, ANRM, RCOND, W, R, INFO );
         chkxer('ZPPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zppcon('U', -1, A, ANRM, RCOND, W, R, INFO );
         chkxer('ZPPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zppcon('U', 1, A, -ANRM, RCOND, W, R, INFO );
         chkxer('ZPPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPPEQU

        srnamc.SRNAMT = 'ZPPEQU';
         infoc.INFOT = 1;
         zppequ('/', 0, A, R1, RCOND, ANRM, INFO );
         chkxer('ZPPEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zppequ('U', -1, A, R1, RCOND, ANRM, INFO );
         chkxer('ZPPEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Test error exits of the routines that use the Cholesky
      // decomposition of a Hermitian positive definite band matrix.

      } else if ( lsamen( 2, C2, 'PB' ) ) {

         // ZPBTRF

        srnamc.SRNAMT = 'ZPBTRF';
         infoc.INFOT = 1;
         zpbtrf('/', 0, 0, A, 1, INFO );
         chkxer('ZPBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpbtrf('U', -1, 0, A, 1, INFO );
         chkxer('ZPBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpbtrf('U', 1, -1, A, 1, INFO );
         chkxer('ZPBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zpbtrf('U', 2, 1, A, 1, INFO );
         chkxer('ZPBTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPBTF2

        srnamc.SRNAMT = 'ZPBTF2';
         infoc.INFOT = 1;
         zpbtf2('/', 0, 0, A, 1, INFO );
         chkxer('ZPBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpbtf2('U', -1, 0, A, 1, INFO );
         chkxer('ZPBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpbtf2('U', 1, -1, A, 1, INFO );
         chkxer('ZPBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zpbtf2('U', 2, 1, A, 1, INFO );
         chkxer('ZPBTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPBTRS

        srnamc.SRNAMT = 'ZPBTRS';
         infoc.INFOT = 1;
         zpbtrs('/', 0, 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpbtrs('U', -1, 0, 0, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpbtrs('U', 1, -1, 0, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zpbtrs('U', 0, 0, -1, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zpbtrs('U', 2, 1, 1, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zpbtrs('U', 2, 0, 1, A, 1, B, 1, INFO );
         chkxer('ZPBTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPBRFS

        srnamc.SRNAMT = 'ZPBRFS';
         infoc.INFOT = 1;
         zpbrfs('/', 0, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpbrfs('U', -1, 0, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpbrfs('U', 1, -1, 0, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zpbrfs('U', 0, 0, -1, A, 1, AF, 1, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zpbrfs('U', 2, 1, 1, A, 1, AF, 2, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zpbrfs('U', 2, 1, 1, A, 2, AF, 1, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zpbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         zpbrfs('U', 2, 0, 1, A, 1, AF, 1, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZPBRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPBCON

        srnamc.SRNAMT = 'ZPBCON';
         infoc.INFOT = 1;
         zpbcon('/', 0, 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpbcon('U', -1, 0, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpbcon('U', 1, -1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zpbcon('U', 2, 1, A, 1, ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zpbcon('U', 1, 0, A, 1, -ANRM, RCOND, W, R, INFO );
         chkxer('ZPBCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPBEQU

        srnamc.SRNAMT = 'ZPBEQU';
         infoc.INFOT = 1;
         zpbequ('/', 0, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpbequ('U', -1, 0, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpbequ('U', 1, -1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zpbequ('U', 2, 1, A, 1, R1, RCOND, ANRM, INFO );
         chkxer('ZPBEQU', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      }
