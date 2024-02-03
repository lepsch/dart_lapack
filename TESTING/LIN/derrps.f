      SUBROUTINE DERRPS( PATH, NUNIT );

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
      double             A( NMAX, NMAX ), WORK( 2*NMAX );
      int                PIV( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DPSTF2, DPSTRF
      // ..
      // .. Scalars in Common ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 110
         for (I = 1; I <= NMAX; I++) { // 100
            A( I, J ) = 1.0 / DBLE( I+J );

         } // 100
         PIV( J ) = J;
         WORK( J ) = 0.0;
         WORK( NMAX+J ) = 0.0;

      } // 110
      OK = true;


         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive semidefinite matrix.

         // DPSTRF

      SRNAMT = 'DPSTRF';
      INFOT = 1;
      dpstrf('/', 0, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTRF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dpstrf('U', -1, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTRF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dpstrf('U', 2, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTRF', INFOT, NOUT, LERR, OK );

         // DPSTF2

      SRNAMT = 'DPSTF2';
      INFOT = 1;
      dpstf2('/', 0, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTF2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dpstf2('U', -1, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTF2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dpstf2('U', 2, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTF2', INFOT, NOUT, LERR, OK );


      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of DERRPS

      }
