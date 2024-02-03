      SUBROUTINE DCHKGK( NIN, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NIN, NOUT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                LDA, LDB, LDVL, LDVR;
      const              LDA = 50, LDB = 50, LDVL = 50, LDVR = 50 ;
      int                LDE, LDF, LDWORK;
      const              LDE = 50, LDF = 50, LDWORK = 50 ;
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IHI, ILO, INFO, J, KNT, M, N, NINFO;
      double             ANORM, BNORM, EPS, RMAX, VMAX;
      // ..
      // .. Local Arrays ..
      int                LMAX( 4 );
      double             A( LDA, LDA ), AF( LDA, LDA ), B( LDB, LDB ), BF( LDB, LDB ), E( LDE, LDE ), F( LDF, LDF ), LSCALE( LDA ), RSCALE( LDA ), VL( LDVL, LDVL ), VLF( LDVL, LDVL ), VR( LDVR, LDVR ), VRF( LDVR, LDVR ), WORK( LDWORK, LDWORK );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DGGBAK, DGGBAL, DLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Initialization

      LMAX( 1 ) = 0
      LMAX( 2 ) = 0
      LMAX( 3 ) = 0
      LMAX( 4 ) = 0
      NINFO = 0
      KNT = 0
      RMAX = ZERO

      EPS = DLAMCH( 'Precision' )

      } // 10
      READ( NIN, FMT = * )N, M
      if (N == 0) GO TO 100;

      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( A( I, J ), J = 1, N )
      } // 20

      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( B( I, J ), J = 1, N )
      } // 30

      for (I = 1; I <= N; I++) { // 40
         READ( NIN, FMT = * )( VL( I, J ), J = 1, M )
      } // 40

      for (I = 1; I <= N; I++) { // 50
         READ( NIN, FMT = * )( VR( I, J ), J = 1, M )
      } // 50

      KNT = KNT + 1

      ANORM = DLANGE( 'M', N, N, A, LDA, WORK )
      BNORM = DLANGE( 'M', N, N, B, LDB, WORK )

      dlacpy('FULL', N, N, A, LDA, AF, LDA );
      dlacpy('FULL', N, N, B, LDB, BF, LDB );

      dggbal('B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO );
      if ( INFO != 0 ) {
         NINFO = NINFO + 1
         LMAX( 1 ) = KNT
      }

      dlacpy('FULL', N, M, VL, LDVL, VLF, LDVL );
      dlacpy('FULL', N, M, VR, LDVR, VRF, LDVR );

      dggbak('B', 'L', N, ILO, IHI, LSCALE, RSCALE, M, VL, LDVL, INFO );
      if ( INFO != 0 ) {
         NINFO = NINFO + 1
         LMAX( 2 ) = KNT
      }

      dggbak('B', 'R', N, ILO, IHI, LSCALE, RSCALE, M, VR, LDVR, INFO );
      if ( INFO != 0 ) {
         NINFO = NINFO + 1
         LMAX( 3 ) = KNT
      }

      // Test of DGGBAK

      // Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
      // where tilde(A) denotes the transformed matrix.

      dgemm('N', 'N', N, M, N, ONE, AF, LDA, VR, LDVR, ZERO, WORK, LDWORK )       CALL DGEMM( 'T', 'N', M, M, N, ONE, VL, LDVL, WORK, LDWORK, ZERO, E, LDE );

      dgemm('N', 'N', N, M, N, ONE, A, LDA, VRF, LDVR, ZERO, WORK, LDWORK )       CALL DGEMM( 'T', 'N', M, M, N, ONE, VLF, LDVL, WORK, LDWORK, ZERO, F, LDF );

      VMAX = ZERO
      for (J = 1; J <= M; J++) { // 70
         for (I = 1; I <= M; I++) { // 60
            VMAX = MAX( VMAX, ABS( E( I, J )-F( I, J ) ) )
         } // 60
      } // 70
      VMAX = VMAX / ( EPS*MAX( ANORM, BNORM ) )
      if ( VMAX > RMAX ) {
         LMAX( 4 ) = KNT
         RMAX = VMAX
      }

      // Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR

      dgemm('N', 'N', N, M, N, ONE, BF, LDB, VR, LDVR, ZERO, WORK, LDWORK )       CALL DGEMM( 'T', 'N', M, M, N, ONE, VL, LDVL, WORK, LDWORK, ZERO, E, LDE );

      dgemm('N', 'N', N, M, N, ONE, B, LDB, VRF, LDVR, ZERO, WORK, LDWORK )       CALL DGEMM( 'T', 'N', M, M, N, ONE, VLF, LDVL, WORK, LDWORK, ZERO, F, LDF );

      VMAX = ZERO
      for (J = 1; J <= M; J++) { // 90
         for (I = 1; I <= M; I++) { // 80
            VMAX = MAX( VMAX, ABS( E( I, J )-F( I, J ) ) )
         } // 80
      } // 90
      VMAX = VMAX / ( EPS*MAX( ANORM, BNORM ) )
      if ( VMAX > RMAX ) {
         LMAX( 4 ) = KNT
         RMAX = VMAX
      }

      GO TO 10

      } // 100

      WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of DGGBAK .. ' )

      WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( ' value of largest test error                  =', D12.3 )
      WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( ' example number where DGGBAL info is not 0    =', I4 )
      WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( ' example number where DGGBAK(L) info is not 0 =', I4 )
      WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( ' example number where DGGBAK(R) info is not 0 =', I4 )
      WRITE( NOUT, FMT = 9994 )LMAX( 4 )
 9994 FORMAT( ' example number having largest error          =', I4 )
      WRITE( NOUT, FMT = 9993 )NINFO
 9993 FORMAT( ' number of examples where info is not 0       =', I4 )
      WRITE( NOUT, FMT = 9992 )KNT
 9992 FORMAT( ' total number of examples tested              =', I4 )

      RETURN

      // End of DCHKGK

      }
