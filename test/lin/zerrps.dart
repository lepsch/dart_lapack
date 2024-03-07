      void zerrps(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                infoc.NUNIT;
      String             PATH;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      int                I, INFO, J, RANK;
      Complex         A( NMAX, NMAX );
      double             RWORK( 2*NMAX );
      int                PIV( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZPSTF2, ZPSTRF
      // ..
      // .. Scalars in Common ..
      int                infoc.INFOT, NOUT;
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 110
         for (I = 1; I <= NMAX; I++) { // 100
            A[I][J] = 1.0 / (I+J).toDouble();

         } // 100
         PIV[J] = J;
         RWORK[J] = 0.0;
         RWORK[NMAX+J] = 0.0;

      } // 110
      infoc.OK.value = true;


         // Test error exits of the routines that use the Cholesky
         // decomposition of an Hermitian positive semidefinite matrix.

         // ZPSTRF

     srnamc.SRNAMT = 'ZPSTRF';
      infoc.INFOT = 1;
      zpstrf('/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('ZPSTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zpstrf('U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('ZPSTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zpstrf('U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('ZPSTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPSTF2

     srnamc.SRNAMT = 'ZPSTF2';
      infoc.INFOT = 1;
      zpstf2('/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('ZPSTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zpstf2('U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('ZPSTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zpstf2('U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO );
      chkxer('ZPSTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );


      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      }
