      SUBROUTINE DERRLQ( PATH, NUNIT )

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
      // EXTERNAL ALAESM, CHKXER, DGELQ2, DGELQF, DORGL2, DORGLQ, DORML2, DORMLQ
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
            A( I, J ) = 1.D0 / DBLE( I+J )
            AF( I, J ) = 1.D0 / DBLE( I+J )
         } // 10
         B( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
      } // 20
      OK = true;

      // Error exits for LQ factorization

      // DGELQF

      SRNAMT = 'DGELQF'
      INFOT = 1
      dgelqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('DGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgelqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('DGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgelqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('DGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dgelqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('DGELQF', INFOT, NOUT, LERR, OK );

      // DGELQ2

      SRNAMT = 'DGELQ2'
      INFOT = 1
      dgelq2(-1, 0, A, 1, B, W, INFO );
      chkxer('DGELQ2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgelq2(0, -1, A, 1, B, W, INFO );
      chkxer('DGELQ2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgelq2(2, 1, A, 1, B, W, INFO );
      chkxer('DGELQ2', INFOT, NOUT, LERR, OK );

      // DORGLQ

      SRNAMT = 'DORGLQ'
      INFOT = 1
      dorglq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorglq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorglq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('DORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorglq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('DORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorglq(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('DORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorglq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('DORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dorglq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('DORGLQ', INFOT, NOUT, LERR, OK );

      // DORGL2

      SRNAMT = 'DORGL2'
      INFOT = 1
      dorgl2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('DORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgl2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('DORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgl2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('DORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgl2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('DORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgl2(1, 1, 2, A, 1, X, W, INFO );
      chkxer('DORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorgl2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('DORGL2', INFOT, NOUT, LERR, OK );

      // DORMLQ

      SRNAMT = 'DORMLQ'
      INFOT = 1
      dormlq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dormlq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dormlq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dormlq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormlq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormlq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormlq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormlq('L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormlq('R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dormlq('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dormlq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dormlq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMLQ', INFOT, NOUT, LERR, OK );

      // DORML2

      SRNAMT = 'DORML2'
      INFOT = 1
      dorml2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorml2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorml2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dorml2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorml2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorml2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorml2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dorml2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dorml2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dorml2('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('DORML2', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of DERRLQ

      }
