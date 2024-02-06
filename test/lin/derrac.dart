import 'common.dart';

      void derrac(NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      int                I, INFO, ITER, J;
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), C( NMAX ), R( NMAX ), R1( NMAX ), R2( NMAX ), W( 2*NMAX ), X( NMAX );
      double             WORK(NMAX*NMAX);
      double               SWORK(NMAX*NMAX);
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DSPOSV
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
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

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = 1.0 / (I+J).toDouble();
            AF[I][J] = 1.0 / (I+J).toDouble();
         } // 10
         B[J] = 0.0;
         R1[J] = 0.0;
         R2[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
         C[J] = 0.0;
         R[J] = 0.0;
      } // 20
      infoc.OK = true;

      srnamc.SRNAMT = 'DSPOSV';
      infoc.INFOT = 1;
      dsposv('/',0,0,A,1,B,1,X,1,WORK,SWORK,ITER,INFO);
      chkxer('DSPOSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dsposv('U',-1,0,A,1,B,1,X,1,WORK,SWORK,ITER,INFO);
      chkxer('DSPOSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dsposv('U',0,-1,A,1,B,1,X,1,WORK,SWORK,ITER,INFO);
      chkxer('DSPOSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dsposv('U',2,1,A,1,B,2,X,2,WORK,SWORK,ITER,INFO);
      chkxer('DSPOSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dsposv('U',2,1,A,2,B,1,X,2,WORK,SWORK,ITER,INFO);
      chkxer('DSPOSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      dsposv('U',2,1,A,2,B,2,X,1,WORK,SWORK,ITER,INFO);
      chkxer('DSPOSV', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      if ( infoc.OK ) {
         WRITE( infoc.NOUT, FMT = 9999 )'DSPOSV';
      } else {
         WRITE( infoc.NOUT, FMT = 9998 )'DSPOSV';
      }

 9999 FORMAT( 1X, A6, ' drivers passed the tests of the error exits' );
 9998 FORMAT( ' *** ', A6, ' drivers failed the tests of the error ', 'exits ***' );

      return;
      }
