      void zerrqp(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 3 ;
      String             C2;
      int                INFO, LW;
      int                IP( NMAX );
      double             RW( 2*NMAX );
      Complex         A( NMAX, NMAX ), TAU( NMAX ), W( 2*NMAX+3*NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQP3
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
      // INTRINSIC DCMPLX

      NOUT = infoc.NUNIT;
      C2 = PATH.substring( 1, 3 );
      LW = NMAX + 1;
      A[1][1] = DCMPLX( 1.0, -1.0 );
      A[1][2] = DCMPLX( 2.0, -2.0 );
      A[2][2] = DCMPLX( 3.0, -3.0 );
      A[2][1] = DCMPLX( 4.0, -4.0 );
      infoc.OK.value = true;
      NOUT.println( * );

      // Test error exits for QR factorization with pivoting

      if ( lsamen( 2, C2, 'QP' ) ) {

         // ZGEQP3

        srnamc.SRNAMT = 'ZGEQP3';
         infoc.INFOT = 1;
         zgeqp3(-1, 0, A, 1, IP, TAU, W, LW, RW, INFO );
         chkxer('ZGEQP3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgeqp3(1, -1, A, 1, IP, TAU, W, LW, RW, INFO );
         chkxer('ZGEQP3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgeqp3(2, 3, A, 1, IP, TAU, W, LW, RW, INFO );
         chkxer('ZGEQP3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgeqp3(2, 2, A, 2, IP, TAU, W, LW-10, RW, INFO );
         chkxer('ZGEQP3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      }
