      void serrps(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NUNIT;
      String             PATH;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      int                I, INFO, J, RANK;
      double               A( NMAX, NMAX ), WORK( 2*NMAX );
      int                PIV( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SPSTF2, SPSTRF
      // ..
      // .. Scalars in Common ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 110
         for (I = 1; I <= NMAX; I++) { // 100
            A[I][J] = 1.0 / REAL( I+J );

         } // 100
         PIV[J] = J;
         WORK[J] = 0.;
         WORK[NMAX+J] = 0.;

      } // 110
      OK = true;


         // Test error exits of the routines that use the Cholesky
         // decomposition of a symmetric positive semidefinite matrix.

         // SPSTRF

     srnamc.SRNAMT = 'SPSTRF';
      INFOT = 1;
      spstrf('/', 0, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('SPSTRF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      spstrf('U', -1, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('SPSTRF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      spstrf('U', 2, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('SPSTRF', INFOT, NOUT, LERR, OK );

         // SPSTF2

     srnamc.SRNAMT = 'SPSTF2';
      INFOT = 1;
      spstf2('/', 0, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('SPSTF2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      spstf2('U', -1, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('SPSTF2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      spstf2('U', 2, A, 1, PIV, RANK, -1.0, WORK, INFO );
      chkxer('SPSTF2', INFOT, NOUT, LERR, OK );


      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
