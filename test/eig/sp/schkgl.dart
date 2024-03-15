      void schkgl(final int NIN, final int NOUT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NIN, NOUT;
      // ..

      int                LDA, LDB, LWORK;
      const              LDA = 20, LDB = 20, LWORK = 6*LDA ;
      double               ZERO;
      const              ZERO = 0.0 ;
      int                I, IHI, IHIIN, ILO, ILOIN, INFO, J, KNT, N, NINFO;
      double               ANORM, BNORM, EPS, RMAX, VMAX;
      int                LMAX( 5 );
      double               A( LDA, LDA ), AIN( LDA, LDA ), B( LDB, LDB ), BIN( LDB, LDB ), LSCALE( LDA ), LSCLIN( LDA ), RSCALE( LDA ), RSCLIN( LDA ), WORK( LWORK );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGGBAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      LMAX[1] = 0;
      LMAX[2] = 0;
      LMAX[3] = 0;
      NINFO = 0;
      KNT = 0;
      RMAX = ZERO;

      EPS = SLAMCH( 'Precision' );

      } // 10

      READ( NIN, FMT = * )N;
      if (N == 0) GO TO 90;
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( A( I, J ), J = 1, N );
      } // 20

      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( B( I, J ), J = 1, N );
      } // 30

      READ( NIN, FMT = * )ILOIN, IHIIN;
      for (I = 1; I <= N; I++) { // 40
         READ( NIN, FMT = * )( AIN( I, J ), J = 1, N );
      } // 40
      for (I = 1; I <= N; I++) { // 50
         READ( NIN, FMT = * )( BIN( I, J ), J = 1, N );
      } // 50

      READ( NIN, FMT = * )( LSCLIN( I ), I = 1, N );
      READ( NIN, FMT = * )( RSCLIN( I ), I = 1, N );

      ANORM = SLANGE( 'M', N, N, A, LDA, WORK );
      BNORM = SLANGE( 'M', N, N, B, LDB, WORK );

      KNT = KNT + 1;

      sggbal('B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO );

      if ( INFO != 0 ) {
         NINFO = NINFO + 1;
         LMAX[1] = KNT;
      }

      if ( ILO != ILOIN || IHI != IHIIN ) {
         NINFO = NINFO + 1;
         LMAX[2] = KNT;
      }

      VMAX = ZERO;
      for (I = 1; I <= N; I++) { // 70
         for (J = 1; J <= N; J++) { // 60
            VMAX = max( VMAX, ABS( A( I, J )-AIN( I, J ) ) );
            VMAX = max( VMAX, ABS( B( I, J )-BIN( I, J ) ) );
         } // 60
      } // 70

      for (I = 1; I <= N; I++) { // 80
         VMAX = max( VMAX, ABS( LSCALE( I )-LSCLIN( I ) ) );
         VMAX = max( VMAX, ABS( RSCALE( I )-RSCLIN( I ) ) );
      } // 80

      VMAX = VMAX / ( EPS*max( ANORM, BNORM ) );

      if ( VMAX > RMAX ) {
         LMAX[3] = KNT;
         RMAX = VMAX;
      }

      GO TO 10;

      } // 90

      WRITE( NOUT, FMT = 9999 );
 9999 FORMAT(' .. test output of SGGBAL .. ' );

      WRITE( NOUT, FMT = 9998 )RMAX;
 9998 FORMAT(' value of largest test error            = ', E12.3 );
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