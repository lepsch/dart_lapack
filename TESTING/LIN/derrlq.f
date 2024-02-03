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
      PARAMETER          ( NMAX = 2 )
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
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
            AF( I, J ) = 1.D0 / DBLE( I+J )
   10    CONTINUE
         B( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
   20 CONTINUE
      OK = .TRUE.

      // Error exits for LQ factorization

      // DGELQF

      SRNAMT = 'DGELQF'
      INFOT = 1
      CALL DGELQF( -1, 0, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'DGELQF', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGELQF( 0, -1, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'DGELQF', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DGELQF( 2, 1, A, 1, B, W, 2, INFO )
      CALL CHKXER( 'DGELQF', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DGELQF( 2, 1, A, 2, B, W, 1, INFO )
      CALL CHKXER( 'DGELQF', INFOT, NOUT, LERR, OK )

      // DGELQ2

      SRNAMT = 'DGELQ2'
      INFOT = 1
      CALL DGELQ2( -1, 0, A, 1, B, W, INFO )
      CALL CHKXER( 'DGELQ2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGELQ2( 0, -1, A, 1, B, W, INFO )
      CALL CHKXER( 'DGELQ2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DGELQ2( 2, 1, A, 1, B, W, INFO )
      CALL CHKXER( 'DGELQ2', INFOT, NOUT, LERR, OK )

      // DORGLQ

      SRNAMT = 'DORGLQ'
      INFOT = 1
      CALL DORGLQ( -1, 0, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGLQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGLQ( 0, -1, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGLQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGLQ( 2, 1, 0, A, 2, X, W, 2, INFO )
      CALL CHKXER( 'DORGLQ', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGLQ( 0, 0, -1, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGLQ', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGLQ( 1, 1, 2, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGLQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORGLQ( 2, 2, 0, A, 1, X, W, 2, INFO )
      CALL CHKXER( 'DORGLQ', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL DORGLQ( 2, 2, 0, A, 2, X, W, 1, INFO )
      CALL CHKXER( 'DORGLQ', INFOT, NOUT, LERR, OK )

      // DORGL2

      SRNAMT = 'DORGL2'
      INFOT = 1
      CALL DORGL2( -1, 0, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGL2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGL2( 0, -1, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGL2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGL2( 2, 1, 0, A, 2, X, W, INFO )
      CALL CHKXER( 'DORGL2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGL2( 0, 0, -1, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGL2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGL2( 1, 1, 2, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGL2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORGL2( 2, 2, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORGL2', INFOT, NOUT, LERR, OK )

      // DORMLQ

      SRNAMT = 'DORMLQ'
      INFOT = 1
      CALL DORMLQ( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORMLQ( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORMLQ( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DORMLQ( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMLQ( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMLQ( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMLQ( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORMLQ( 'L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORMLQ( 'R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DORMLQ( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL DORMLQ( 'L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL DORMLQ( 'R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
      CALL CHKXER( 'DORMLQ', INFOT, NOUT, LERR, OK )

      // DORML2

      SRNAMT = 'DORML2'
      INFOT = 1
      CALL DORML2( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORML2( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORML2( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DORML2( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORML2( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORML2( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORML2( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORML2( 'L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORML2( 'R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DORML2( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORML2', INFOT, NOUT, LERR, OK )

      // Print a summary line.

      CALL ALAESM( PATH, OK, NOUT )

      RETURN

      // End of DERRLQ

      END
