import 'common.dart';
      void derrps(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NUNIT;
      String             PATH;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      int                I, INFO, J, RANK;
      double             A( NMAX, NMAX ), WORK( 2*NMAX );
      int                PIV( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DPSTF2, DPSTRF
      // ..
      // .. Scalars in Common ..
      int                infoc.INFOT, infoc.NOUT;
      bool               infoc.LERR, infoc.OK;
      String             srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 110
         for (I = 1; I <= NMAX; I++) { // 100
            A[I][J] = 1.0 / (I+J).toDouble();

         } // 100
         PIV[J] = J;
         WORK[J] = 0.0;
         WORK[NMAX+J] = 0.0;

      } // 110
      infoc.OK = true;


         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive semidefinite matrix.

         // DPSTRF

      srnamc.SRNAMT = 'DPSTRF';
      infoc.INFOT = 1;
      dpstrf('/', 0, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dpstrf('U', -1, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dpstrf('U', 2, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DPSTF2

      srnamc.SRNAMT = 'DPSTF2';
      infoc.INFOT = 1;
      dpstf2('/', 0, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dpstf2('U', -1, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dpstf2('U', 2, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('DPSTF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );


      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      }
