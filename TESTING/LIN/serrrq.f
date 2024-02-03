      SUBROUTINE SERRRQ( PATH, NUNIT )

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
      REAL               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGERQ2, SGERQF, SGERQS, SORGR2, SORGRQ, SORMR2, SORMRQ
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
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1. / REAL( I+J )
            AF( I, J ) = 1. / REAL( I+J )
         } // 10
         B( J ) = 0.
         W( J ) = 0.
         X( J ) = 0.
      } // 20
      OK = true;

      // Error exits for RQ factorization

      // SGERQF

      SRNAMT = 'SGERQF'
      INFOT = 1
      sgerqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('SGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgerqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('SGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sgerqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('SGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sgerqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('SGERQF', INFOT, NOUT, LERR, OK );

      // SGERQ2

      SRNAMT = 'SGERQ2'
      INFOT = 1
      sgerq2(-1, 0, A, 1, B, W, INFO );
      chkxer('SGERQ2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgerq2(0, -1, A, 1, B, W, INFO );
      chkxer('SGERQ2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sgerq2(2, 1, A, 1, B, W, INFO );
      chkxer('SGERQ2', INFOT, NOUT, LERR, OK );

      // SGERQS

      SRNAMT = 'SGERQS'
      INFOT = 1
      sgerqs(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgerqs(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgerqs(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('SGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sgerqs(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sgerqs(2, 2, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('SGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 8
      sgerqs(2, 2, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('SGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 10
      sgerqs(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGERQS', INFOT, NOUT, LERR, OK );

      // SORGRQ

      SRNAMT = 'SORGRQ'
      INFOT = 1
      sorgrq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorgrq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorgrq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('SORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorgrq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('SORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorgrq(1, 2, 2, A, 1, X, W, 1, INFO );
      chkxer('SORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sorgrq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('SORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      sorgrq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('SORGRQ', INFOT, NOUT, LERR, OK );

      // SORGR2

      SRNAMT = 'SORGR2'
      INFOT = 1
      sorgr2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('SORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorgr2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('SORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorgr2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('SORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorgr2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('SORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorgr2(1, 2, 2, A, 2, X, W, INFO );
      chkxer('SORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sorgr2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('SORGR2', INFOT, NOUT, LERR, OK );

      // SORMRQ

      SRNAMT = 'SORMRQ'
      INFOT = 1
      sormrq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sormrq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sormrq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sormrq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormrq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormrq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormrq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sormrq('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sormrq('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 10
      sormrq('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      sormrq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      sormrq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('SORMRQ', INFOT, NOUT, LERR, OK );

      // SORMR2

      SRNAMT = 'SORMR2'
      INFOT = 1
      sormr2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sormr2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sormr2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sormr2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormr2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormr2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormr2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sormr2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sormr2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 10
      sormr2('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORMR2', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of SERRRQ

      }
