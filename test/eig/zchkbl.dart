      void zchkbl(NIN, final int NOUT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NIN, NOUT;
      // ..

// ======================================================================

      // .. Parameters ..
      int                LDA;
      const              LDA = 20 ;
      double             ZERO;
      const              ZERO = 0.0 ;
      int                I, IHI, IHIIN, ILO, ILOIN, INFO, J, KNT, N, NINFO;
      double             ANORM, MEPS, RMAX, SFMIN, TEMP, VMAX;
      Complex         CDUM;
      int                LMAX( 3 );
      double             DUMMY( 1 ), SCALE( LDA ), SCALIN( LDA );
      Complex         A( LDA, LDA ), AIN( LDA, LDA );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEBAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[CDUM] = ( CDUM.toDouble() ).abs() + ( DIMAG( CDUM ) ).abs();

      LMAX[1] = 0;
      LMAX[2] = 0;
      LMAX[3] = 0;
      NINFO = 0;
      KNT = 0;
      RMAX = ZERO;
      VMAX = ZERO;
      SFMIN = dlamch( 'S' );
      MEPS = dlamch( 'E' );

      } // 10

      READ( NIN, FMT = * )N;
      if (N == 0) GO TO 70;
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( A( I, J ), J = 1, N );
      } // 20

      READ( NIN, FMT = * )ILOIN, IHIIN;
      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( AIN( I, J ), J = 1, N );
      } // 30
      READ( NIN, FMT = * )( SCALIN( I ), I = 1, N );

      ANORM = ZLANGE( 'M', N, N, A, LDA, DUMMY );
      KNT = KNT + 1;
      zgebal('B', N, A, LDA, ILO, IHI, SCALE, INFO );

      if ( INFO != 0 ) {
         NINFO = NINFO + 1;
         LMAX[1] = KNT;
      }

      if ( ILO != ILOIN || IHI != IHIIN ) {
         NINFO = NINFO + 1;
         LMAX[2] = KNT;
      }

      for (I = 1; I <= N; I++) { // 50
         for (J = 1; J <= N; J++) { // 40
            TEMP = max( CABS1( A( I, J ) ), CABS1( AIN( I, J ) ) );
            TEMP = max( TEMP, SFMIN );
            VMAX = max( VMAX, CABS1( A( I, J )-AIN( I, J ) ) / TEMP );
         } // 40
      } // 50

      for (I = 1; I <= N; I++) { // 60
         TEMP = max( SCALE( I ), SCALIN( I ) );
         TEMP = max( TEMP, SFMIN );
         VMAX = max( VMAX, ABS( SCALE( I )-SCALIN( I ) ) / TEMP );
      } // 60

      if ( VMAX > RMAX ) {
         LMAX[3] = KNT;
         RMAX = VMAX;
      }

      GO TO 10;

      } // 70

      WRITE( NOUT, FMT = 9999 );
 9999 FORMAT(' .. test output of ZGEBAL .. ' );

      WRITE( NOUT, FMT = 9998 )RMAX;
 9998 FORMAT(' value of largest test error            = ', D12.3 );
      WRITE( NOUT, FMT = 9997 )LMAX( 1 );
 9997 FORMAT(' example number where info is not zero  = ${.i4}');
      WRITE( NOUT, FMT = 9996 )LMAX( 2 );
 9996 FORMAT(' example number where ILO or IHI wrong  = ${.i4}');
      WRITE( NOUT, FMT = 9995 )LMAX( 3 );
 9995 FORMAT(' example number having largest error    = ${.i4}');
      WRITE( NOUT, FMT = 9994 )NINFO;
 9994 FORMAT(' number of examples where info is not 0 = ${.i4}');
      WRITE( NOUT, FMT = 9993 )KNT;
 9993 FORMAT(' total number of examples tested        = ${.i4}');

      }
