      SUBROUTINE DCHKBK( NIN, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NIN, NOUT;
      // ..

* ======================================================================

      // .. Parameters ..
      int                LDE;
      const              LDE = 20 ;
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, IHI, ILO, INFO, J, KNT, N, NINFO;
      double             EPS, RMAX, SAFMIN, VMAX, X;
      // ..
      // .. Local Arrays ..
      int                LMAX( 2 );
      double             E( LDE, LDE ), EIN( LDE, LDE ), SCALE( LDE );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEBAK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      LMAX( 1 ) = 0
      LMAX( 2 ) = 0
      NINFO = 0
      KNT = 0
      RMAX = ZERO
      EPS = DLAMCH( 'E' )
      SAFMIN = DLAMCH( 'S' )

      } // 10

      READ( NIN, FMT = * )N, ILO, IHI
      if (N.EQ.0) GO TO 60;

      READ( NIN, FMT = * )( SCALE( I ), I = 1, N )
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( E( I, J ), J = 1, N )
      } // 20

      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( EIN( I, J ), J = 1, N )
      } // 30

      KNT = KNT + 1
      dgebak('B', 'R', N, ILO, IHI, SCALE, N, E, LDE, INFO );

      if ( INFO.NE.0 ) {
         NINFO = NINFO + 1
         LMAX( 1 ) = KNT
      }

      VMAX = ZERO
      for (I = 1; I <= N; I++) { // 50
         for (J = 1; J <= N; J++) { // 40
            X = ABS( E( I, J )-EIN( I, J ) ) / EPS
            IF( ABS( E( I, J ) ).GT.SAFMIN ) X = X / ABS( E( I, J ) )
            VMAX = MAX( VMAX, X )
         } // 40
      } // 50

      if ( VMAX.GT.RMAX ) {
         LMAX( 2 ) = KNT
         RMAX = VMAX
      }

      GO TO 10

      } // 60

      WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of DGEBAK .. ' )

      WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( 1X, 'value of largest test error             = ', D12.3 )
      WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( 1X, 'example number where info is not zero   = ', I4 )
      WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( 1X, 'example number having largest error     = ', I4 )
      WRITE( NOUT, FMT = 9995 )NINFO
 9995 FORMAT( 1X, 'number of examples where info is not 0  = ', I4 )
      WRITE( NOUT, FMT = 9994 )KNT
 9994 FORMAT( 1X, 'total number of examples tested         = ', I4 )

      RETURN

      // End of DCHKBK

      }
