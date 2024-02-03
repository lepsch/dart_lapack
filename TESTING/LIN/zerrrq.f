      SUBROUTINE ZERRRQ( PATH, NUNIT )

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
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGERQ2, ZGERQF, ZGERQS, ZUNGR2, ZUNGRQ, ZUNMR2, ZUNMRQ
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

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )
         } // 10
         B( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
      } // 20
      OK = true;

      // Error exits for RQ factorization

      // ZGERQF

      SRNAMT = 'ZGERQF'
      INFOT = 1
      zgerqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('ZGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgerqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('ZGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zgerqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('ZGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 7
      zgerqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('ZGERQF', INFOT, NOUT, LERR, OK );

      // ZGERQ2

      SRNAMT = 'ZGERQ2'
      INFOT = 1
      zgerq2(-1, 0, A, 1, B, W, INFO );
      chkxer('ZGERQ2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgerq2(0, -1, A, 1, B, W, INFO );
      chkxer('ZGERQ2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zgerq2(2, 1, A, 1, B, W, INFO );
      chkxer('ZGERQ2', INFOT, NOUT, LERR, OK );

      // ZGERQS

      SRNAMT = 'ZGERQS'
      INFOT = 1
      zgerqs(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgerqs(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgerqs(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('ZGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zgerqs(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zgerqs(2, 2, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('ZGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 8
      zgerqs(2, 2, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('ZGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 10
      zgerqs(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGERQS', INFOT, NOUT, LERR, OK );

      // ZUNGRQ

      SRNAMT = 'ZUNGRQ'
      INFOT = 1
      zungrq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zungrq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zungrq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('ZUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zungrq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('ZUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zungrq(1, 2, 2, A, 1, X, W, 1, INFO );
      chkxer('ZUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zungrq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('ZUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      zungrq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('ZUNGRQ', INFOT, NOUT, LERR, OK );

      // ZUNGR2

      SRNAMT = 'ZUNGR2'
      INFOT = 1
      zungr2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('ZUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zungr2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('ZUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zungr2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('ZUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zungr2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('ZUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zungr2(1, 2, 2, A, 2, X, W, INFO );
      chkxer('ZUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zungr2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('ZUNGR2', INFOT, NOUT, LERR, OK );

      // ZUNMRQ

      SRNAMT = 'ZUNMRQ'
      INFOT = 1
      zunmrq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zunmrq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zunmrq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zunmrq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zunmrq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zunmrq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zunmrq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      zunmrq('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      zunmrq('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 10
      zunmrq('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      zunmrq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      zunmrq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMRQ', INFOT, NOUT, LERR, OK );

      // ZUNMR2

      SRNAMT = 'ZUNMR2'
      INFOT = 1
      zunmr2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zunmr2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zunmr2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zunmr2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zunmr2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zunmr2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zunmr2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      zunmr2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      zunmr2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 10
      zunmr2('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNMR2', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of ZERRRQ

      }
