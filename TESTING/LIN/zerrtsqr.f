      SUBROUTINE ZERRTSQR( PATH, NUNIT )
      IMPLICIT NONE

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
      int                I, INFO, J, MB, NB;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX ), TAU(NMAX)
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQR, ZGEMQR, ZGELQ, ZGEMLQ
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

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.D0 / DBLE( I+J )
            C( I, J ) = 1.D0 / DBLE( I+J )
            T( I, J ) = 1.D0 / DBLE( I+J )
         }
         W( J ) = 0.D0
      }
      OK = true;

      // Error exits for TS factorization

      // ZGEQR

      SRNAMT = 'ZGEQR'
      INFOT = 1
      zgeqr(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgeqr(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zgeqr(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('ZGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 6
      zgeqr(3, 2, A, 3, TAU, 1, W, 1, INFO );
      chkxer('ZGEQR', INFOT, NOUT, LERR, OK );
      INFOT = 8
      zgeqr(3, 2, A, 3, TAU, 8, W, 0, INFO );
      chkxer('ZGEQR', INFOT, NOUT, LERR, OK );

      // ZLATSQR

      MB = 1
      NB = 1
      SRNAMT = 'ZLATSQR'
      INFOT = 1
      zlatsqr(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zlatsqr(1, 2, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', INFOT, NOUT, LERR, OK );
      zlatsqr(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zlatsqr(2, 1, -1, NB, A, 2, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zlatsqr(2, 1, MB, 2, A, 2, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 6
      zlatsqr(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 8
      zlatsqr(2, 1, MB, NB, A, 2, TAU, 0, W, 1, INFO );
      chkxer('ZLATSQR', INFOT, NOUT, LERR, OK );
      INFOT = 10
      zlatsqr(2, 1, MB, NB, A, 2, TAU, 2, W, 0, INFO );
      chkxer('ZLATSQR', INFOT, NOUT, LERR, OK );

      // ZGEMQR

      TAU(1)=1
      TAU(2)=1
      SRNAMT = 'ZGEMQR'
      NB=1
      INFOT = 1
      zgemqr('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgemqr('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zgemqr('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zgemqr('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zgemqr('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zgemqr('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 7
      zgemqr('L', 'N', 2, 1, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 9
      zgemqr('R', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 9
      zgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 0, C, 1, W, 1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 11
      zgemqr('L', 'N', 2, 1, 1, A, 2, TAU, 6, C, 0, W, 1,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );
      INFOT = 13
      zgemqr('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('ZGEMQR', INFOT, NOUT, LERR, OK );

      // ZGELQ

      SRNAMT = 'ZGELQ'
      INFOT = 1
      zgelq(-1, 0, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgelq(0, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zgelq(1, 1, A, 0, TAU, 1, W, 1, INFO );
      chkxer('ZGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 6
      zgelq(2, 3, A, 3, TAU, 1, W, 1, INFO );
      chkxer('ZGELQ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      zgelq(2, 3, A, 3, TAU, 8, W, 0, INFO );
      chkxer('ZGELQ', INFOT, NOUT, LERR, OK );

      // ZLASWLQ

      MB = 1
      NB = 1
      SRNAMT = 'ZLASWLQ'
      INFOT = 1
      zlaswlq(-1, 0, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zlaswlq(2, 1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );
      zlaswlq(0, -1, MB, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zlaswlq(1, 2, -1, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );
      zlaswlq(1, 1, 2, NB, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zlaswlq(1, 2, MB, -1, A, 1, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 6
      zlaswlq(1, 2, MB, NB, A, 0, TAU, 1, W, 1, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      zlaswlq(1, 2, MB, NB, A, 1, TAU, 0, W, 1, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );
      INFOT = 10
      zlaswlq(1, 2, MB, NB, A, 1, TAU, 1, W, 0, INFO );
      chkxer('ZLASWLQ', INFOT, NOUT, LERR, OK );

      // ZGEMLQ

      TAU(1)=1
      TAU(2)=1
      SRNAMT = 'ZGEMLQ'
      NB=1
      INFOT = 1
      zgemlq('/', 'N', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zgemlq('L', '/', 0, 0, 0, A, 1, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zgemlq('L', 'N', -1, 0, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zgemlq('L', 'N', 0, -1, 0, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zgemlq('L', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zgemlq('R', 'N', 0, 0, -1, A, 1, TAU, 1, C, 1, W,1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      zgemlq('L', 'N', 1, 2, 0, A, 0, TAU, 1, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 9
      zgemlq('R', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 9
      zgemlq('L', 'N', 2, 2, 1, A, 1, TAU, 0, C, 1, W, 1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 11
      zgemlq('L', 'N', 1, 2, 1, A, 1, TAU, 6, C, 0, W, 1,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 13
      zgemlq('L', 'N', 2, 2, 1, A, 2, TAU, 6, C, 2, W, 0,INFO);
      chkxer('ZGEMLQ', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of ZERRTSQR

      }
