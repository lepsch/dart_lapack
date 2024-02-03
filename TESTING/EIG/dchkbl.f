      SUBROUTINE DCHKBL( NIN, NOUT )

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
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IHI, IHIIN, ILO, ILOIN, INFO, J, KNT, N, NINFO;
      double             ANORM, MEPS, RMAX, SFMIN, TEMP, VMAX;
      // ..
      // .. Local Arrays ..
      int                LMAX( 3 );
      double             A( LDA, LDA ), AIN( LDA, LDA ), DUMMY( 1 ), SCALE( LDA ), SCALIN( LDA );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEBAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      LMAX( 1 ) = 0
      LMAX( 2 ) = 0
      LMAX( 3 ) = 0
      NINFO = 0
      KNT = 0
      RMAX = ZERO
      VMAX = ZERO
      SFMIN = DLAMCH( 'S' )
      MEPS = DLAMCH( 'E' )

      } // 10

      READ( NIN, FMT = * )N
      IF( N.EQ.0 ) GO TO 70
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( A( I, J ), J = 1, N )
      } // 20

      READ( NIN, FMT = * )ILOIN, IHIIN
      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( AIN( I, J ), J = 1, N )
      } // 30
      READ( NIN, FMT = * )( SCALIN( I ), I = 1, N )

      ANORM = DLANGE( 'M', N, N, A, LDA, DUMMY )
      KNT = KNT + 1

      dgebal('B', N, A, LDA, ILO, IHI, SCALE, INFO );

      if ( INFO.NE.0 ) {
         NINFO = NINFO + 1
         LMAX( 1 ) = KNT
      }

      if ( ILO.NE.ILOIN .OR. IHI.NE.IHIIN ) {
         NINFO = NINFO + 1
         LMAX( 2 ) = KNT
      }

      for (I = 1; I <= N; I++) { // 50
         for (J = 1; J <= N; J++) { // 40
            TEMP = MAX( A( I, J ), AIN( I, J ) )
            TEMP = MAX( TEMP, SFMIN )
            VMAX = MAX( VMAX, ABS( A( I, J )-AIN( I, J ) ) / TEMP )
         } // 40
      } // 50

      for (I = 1; I <= N; I++) { // 60
         TEMP = MAX( SCALE( I ), SCALIN( I ) )
         TEMP = MAX( TEMP, SFMIN )
         VMAX = MAX( VMAX, ABS( SCALE( I )-SCALIN( I ) ) / TEMP )
      } // 60


      if ( VMAX.GT.RMAX ) {
         LMAX( 3 ) = KNT
         RMAX = VMAX
      }

      GO TO 10

      } // 70

      WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of DGEBAL .. ' )

      WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( 1X, 'value of largest test error            = ', D12.3 )
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

      // End of DCHKBL

      }
