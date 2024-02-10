      void cerrqp(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 3 ;
      String             C2;
      int                INFO, LW;
      int                IP( NMAX );
      double               RW( 2*NMAX );
      Complex            A( NMAX, NMAX ), TAU( NMAX ), W( 2*NMAX+3*NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CGEQP3, CHKXER
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
      // INTRINSIC CMPLX

      NOUT = NUNIT;
      C2 = PATH( 2: 3 );
      LW = NMAX + 1;
      A[1][1] = CMPLX( 1.0, -1.0 );
      A[1][2] = CMPLX( 2.0, -2.0 );
      A[2][2] = CMPLX( 3.0, -3.0 );
      A[2][1] = CMPLX( 4.0, -4.0 );
      OK = true;
      WRITE( NOUT, FMT = * );

      // Test error exits for QR factorization with pivoting

      if ( lsamen( 2, C2, 'QP' ) ) {

         // CGEQP3

        srnamc.SRNAMT = 'CGEQP3';
         INFOT = 1;
         cgeqp3(-1, 0, A, 1, IP, TAU, W, LW, RW, INFO );
         chkxer('CGEQP3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgeqp3(1, -1, A, 1, IP, TAU, W, LW, RW, INFO );
         chkxer('CGEQP3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgeqp3(2, 3, A, 1, IP, TAU, W, LW, RW, INFO );
         chkxer('CGEQP3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgeqp3(2, 2, A, 2, IP, TAU, W, LW-10, RW, INFO );
         chkxer('CGEQP3', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }
