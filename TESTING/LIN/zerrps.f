      SUBROUTINE ZERRPS( PATH, NUNIT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NUNIT;
      String             PATH;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, RANK;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX )
      double             RWORK( 2*NMAX );
      int                PIV( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZPSTF2, ZPSTRF
      // ..
      // .. Scalars in Common ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String             SRNAMT;
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

      for (J = 1; J <= NMAX; J++) { // 110
         for (I = 1; I <= NMAX; I++) { // 100
            A( I, J ) = 1.D0 / DBLE( I+J )

         } // 100
         PIV( J ) = J
         RWORK( J ) = 0.D0
         RWORK( NMAX+J ) = 0.D0

      } // 110
      OK = .TRUE.


         // Test error exits of the routines that use the Cholesky
         // decomposition of an Hermitian positive semidefinite matrix.

         // ZPSTRF

      SRNAMT = 'ZPSTRF'
      INFOT = 1
      zpstrf('/', 0, A, 1, PIV, RANK, -1.D0, RWORK, INFO );
      chkxer('ZPSTRF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zpstrf('U', -1, A, 1, PIV, RANK, -1.D0, RWORK, INFO );
      chkxer('ZPSTRF', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zpstrf('U', 2, A, 1, PIV, RANK, -1.D0, RWORK, INFO );
      chkxer('ZPSTRF', INFOT, NOUT, LERR, OK );

         // ZPSTF2

      SRNAMT = 'ZPSTF2'
      INFOT = 1
      zpstf2('/', 0, A, 1, PIV, RANK, -1.D0, RWORK, INFO );
      chkxer('ZPSTF2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zpstf2('U', -1, A, 1, PIV, RANK, -1.D0, RWORK, INFO );
      chkxer('ZPSTF2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zpstf2('U', 2, A, 1, PIV, RANK, -1.D0, RWORK, INFO );
      chkxer('ZPSTF2', INFOT, NOUT, LERR, OK );


      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of ZERRPS

      }
