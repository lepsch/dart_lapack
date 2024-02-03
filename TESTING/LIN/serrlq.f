      SUBROUTINE SERRLQ( PATH, NUNIT )

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
      // EXTERNAL ALAESM, CHKXER, SGELQ2, SGELQF, SORGL2, SORGLQ, SORML2, SORMLQ
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1. / REAL( I+J )
            AF( I, J ) = 1. / REAL( I+J )
   10    CONTINUE
         B( J ) = 0.
         W( J ) = 0.
         X( J ) = 0.
   20 CONTINUE
      OK = .TRUE.

      // Error exits for LQ factorization

      // SGELQF

      SRNAMT = 'SGELQF'
      INFOT = 1
      sgelqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('SGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgelqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('SGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sgelqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('SGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sgelqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('SGELQF', INFOT, NOUT, LERR, OK );

      // SGELQ2

      SRNAMT = 'SGELQ2'
      INFOT = 1
      sgelq2(-1, 0, A, 1, B, W, INFO );
      chkxer('SGELQ2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgelq2(0, -1, A, 1, B, W, INFO );
      chkxer('SGELQ2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sgelq2(2, 1, A, 1, B, W, INFO );
      chkxer('SGELQ2', INFOT, NOUT, LERR, OK );

      // SORGLQ

      SRNAMT = 'SORGLQ'
      INFOT = 1
      sorglq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorglq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorglq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('SORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorglq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('SORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorglq(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('SORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sorglq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('SORGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      sorglq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('SORGLQ', INFOT, NOUT, LERR, OK );

      // SORGL2

      SRNAMT = 'SORGL2'
      INFOT = 1
      sorgl2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('SORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorgl2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('SORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorgl2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('SORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorgl2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('SORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorgl2(1, 1, 2, A, 1, X, W, INFO );
      chkxer('SORGL2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sorgl2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('SORGL2', INFOT, NOUT, LERR, OK );

      // SORMLQ

      SRNAMT = 'SORMLQ'
      INFOT = 1
      sormlq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sormlq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sormlq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sormlq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormlq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormlq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sormlq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sormlq('L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sormlq('R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 10
      sormlq('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      sormlq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 12
      sormlq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('SORMLQ', INFOT, NOUT, LERR, OK );

      // SORML2

      SRNAMT = 'SORML2'
      INFOT = 1
      sorml2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sorml2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sorml2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sorml2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sorml2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sorml2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sorml2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sorml2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sorml2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );
      INFOT = 10
      sorml2('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('SORML2', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of SERRLQ

      }
