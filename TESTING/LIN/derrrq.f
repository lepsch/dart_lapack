      SUBROUTINE DERRRQ( PATH, NUNIT )

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
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGERQ2, DGERQF, DGERQS, DORGR2, DORGRQ, DORMR2, DORMRQ
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

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1.0 / DBLE( I+J )
            AF( I, J ) = 1.0 / DBLE( I+J )
         } // 10
         B( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
      } // 20
      OK = true;

      // Error exits for RQ factorization

      // DGERQF

      SRNAMT = 'DGERQF'
      INFOT = 1
      dgerqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('DGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgerqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('DGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgerqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('DGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dgerqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('DGERQF', INFOT, NOUT, LERR, OK );

      // DGERQ2

      SRNAMT = 'DGERQ2'
      INFOT = 1
      dgerq2(-1, 0, A, 1, B, W, INFO );
      chkxer('DGERQ2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgerq2(0, -1, A, 1, B, W, INFO );
      chkxer('DGERQ2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgerq2(2, 1, A, 1, B, W, INFO );
      chkxer('DGERQ2', INFOT, NOUT, LERR, OK );

      // DGERQS

      SRNAMT = 'DGERQS'
      INFOT = 1
      dgerqs(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgerqs(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgerqs(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dgerqs(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgerqs(2, 2, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('DGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dgerqs(2, 2, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dgerqs(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', INFOT, NOUT, LERR, OK );

      // DORGRQ

      SRNAMT = 'DORGRQ'
      INFOT = 1
      dorgrq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgrq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgrq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('DORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgrq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('DORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgrq(1, 2, 2, A, 1, X, W, 1, INFO );
      chkxer('DORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorgrq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('DORGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dorgrq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('DORGRQ', INFOT, NOUT, LERR, OK );

      // DORGR2

      SRNAMT = 'DORGR2'
      INFOT = 1
      dorgr2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('DORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgr2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('DORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgr2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('DORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgr2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('DORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgr2(1, 2, 2, A, 2, X, W, INFO );
      chkxer('DORGR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorgr2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('DORGR2', INFOT, NOUT, LERR, OK );

      // DORMRQ

      SRNAMT = 'DORMRQ'
      INFOT = 1
      dormrq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dormrq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dormrq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dormrq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormrq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormrq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormrq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormrq('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormrq('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dormrq('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dormrq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dormrq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMRQ', INFOT, NOUT, LERR, OK );

      // DORMR2

      SRNAMT = 'DORMR2'
      INFOT = 1
      dormr2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dormr2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dormr2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dormr2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormr2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormr2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormr2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormr2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormr2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dormr2('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of DERRRQ

      }
