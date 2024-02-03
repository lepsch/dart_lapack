      SUBROUTINE CCHKBL( NIN, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NIN, NOUT;
      // ..

* ======================================================================

      // .. Parameters ..
      int                LDA;
      const              LDA = 20 ;
      REAL               ZERO
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IHI, IHIIN, ILO, ILOIN, INFO, J, KNT, N, NINFO;
      REAL               ANORM, MEPS, RMAX, SFMIN, TEMP, VMAX
      COMPLEX            CDUM
      // ..
      // .. Local Arrays ..
      int                LMAX( 3 );
      REAL               DUMMY( 1 ), SCALE( LDA ), SCALIN( LDA )
      COMPLEX            A( LDA, LDA ), AIN( LDA, LDA )
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEBAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      LMAX( 1 ) = 0
      LMAX( 2 ) = 0
      LMAX( 3 ) = 0
      NINFO = 0
      KNT = 0
      RMAX = ZERO
      VMAX = ZERO
      SFMIN = SLAMCH( 'S' )
      MEPS = SLAMCH( 'E' )

      } // 10

      READ( NIN, FMT = * )N
      if (N == 0) GO TO 70;
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( A( I, J ), J = 1, N )
      } // 20

      READ( NIN, FMT = * )ILOIN, IHIIN
      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( AIN( I, J ), J = 1, N )
      } // 30
      READ( NIN, FMT = * )( SCALIN( I ), I = 1, N )

      ANORM = CLANGE( 'M', N, N, A, LDA, DUMMY )
      KNT = KNT + 1
      cgebal('B', N, A, LDA, ILO, IHI, SCALE, INFO );

      if ( INFO != 0 ) {
         NINFO = NINFO + 1
         LMAX( 1 ) = KNT
      }

      if ( ILO != ILOIN || IHI != IHIIN ) {
         NINFO = NINFO + 1
         LMAX( 2 ) = KNT
      }

      for (I = 1; I <= N; I++) { // 50
         for (J = 1; J <= N; J++) { // 40
            TEMP = MAX( CABS1( A( I, J ) ), CABS1( AIN( I, J ) ) )
            TEMP = MAX( TEMP, SFMIN )
            VMAX = MAX( VMAX, CABS1( A( I, J )-AIN( I, J ) ) / TEMP )
         } // 40
      } // 50

      for (I = 1; I <= N; I++) { // 60
         TEMP = MAX( SCALE( I ), SCALIN( I ) )
         TEMP = MAX( TEMP, SFMIN )
         VMAX = MAX( VMAX, ABS( SCALE( I )-SCALIN( I ) ) / TEMP )
      } // 60

      if ( VMAX > RMAX ) {
         LMAX( 3 ) = KNT
         RMAX = VMAX
      }

      GO TO 10

      } // 70

      WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of CGEBAL .. ' )

      WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( 1X, 'value of largest test error            = ', E12.3 )
      WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( 1X, 'example number where info is not zero  = ', I4 )
      WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( 1X, 'example number where ILO or IHI wrong  = ', I4 )
      WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( 1X, 'example number having largest error    = ', I4 )
      WRITE( NOUT, FMT = 9994 )NINFO
 9994 FORMAT( 1X, 'number of examples where info is not 0 = ', I4 )
      WRITE( NOUT, FMT = 9993 )KNT
 9993 FORMAT( 1X, 'total number of examples tested        = ', I4 )

      RETURN

      // End of CCHKBL

      }
