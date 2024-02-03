      void zchkbk(NIN, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NIN, NOUT;
      // ..

// ======================================================================

      // .. Parameters ..
      int                LDE;
      const              LDE = 20 ;
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IHI, ILO, INFO, J, KNT, N, NINFO;
      double             EPS, RMAX, SAFMIN, VMAX, X;
      Complex         CDUM;
      // ..
      // .. Local Arrays ..
      int                LMAX( 2 );
      double             SCALE( LDE );
      Complex         E( LDE, LDE ), EIN( LDE, LDE );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEBAK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) );
      // ..
      // .. Executable Statements ..

      LMAX( 1 ) = 0;
      LMAX( 2 ) = 0;
      NINFO = 0;
      KNT = 0;
      RMAX = ZERO;
      EPS = DLAMCH( 'E' );
      SAFMIN = DLAMCH( 'S' );

      } // 10

      READ( NIN, FMT = * )N, ILO, IHI;
      if (N == 0) GO TO 60;

      READ( NIN, FMT = * )( SCALE( I ), I = 1, N );
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( E( I, J ), J = 1, N );
      } // 20

      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( EIN( I, J ), J = 1, N );
      } // 30

      KNT = KNT + 1;
      zgebak('B', 'R', N, ILO, IHI, SCALE, N, E, LDE, INFO );

      if ( INFO != 0 ) {
         NINFO = NINFO + 1;
         LMAX( 1 ) = KNT;
      }

      VMAX = ZERO;
      for (I = 1; I <= N; I++) { // 50
         for (J = 1; J <= N; J++) { // 40
            X = CABS1( E( I, J )-EIN( I, J ) ) / EPS;
            if( CABS1( E( I, J ) ) > SAFMIN ) X = X / CABS1( E( I, J ) );
            VMAX = max( VMAX, X );
         } // 40
      } // 50

      if ( VMAX > RMAX ) {
         LMAX( 2 ) = KNT;
         RMAX = VMAX;
      }

      GO TO 10;

      } // 60

      WRITE( NOUT, FMT = 9999 );
 9999 FORMAT( 1X, '.. test output of ZGEBAK .. ' );

      WRITE( NOUT, FMT = 9998 )RMAX;
 9998 FORMAT( 1X, 'value of largest test error             = ', D12.3 );
      WRITE( NOUT, FMT = 9997 )LMAX( 1 );
 9997 FORMAT( 1X, 'example number where info is not zero   = ', I4 );
      WRITE( NOUT, FMT = 9996 )LMAX( 2 );
 9996 FORMAT( 1X, 'example number having largest error     = ', I4 );
      WRITE( NOUT, FMT = 9995 )NINFO;
 9995 FORMAT( 1X, 'number of examples where info is not 0  = ', I4 );
      WRITE( NOUT, FMT = 9994 )KNT;
 9994 FORMAT( 1X, 'total number of examples tested         = ', I4 );

      return;

      // End of ZCHKBK

      }
