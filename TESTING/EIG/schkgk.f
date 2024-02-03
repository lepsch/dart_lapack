      SUBROUTINE SCHKGK( NIN, NOUT )

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
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IHI, ILO, INFO, J, KNT, M, N, NINFO;
      REAL               ANORM, BNORM, EPS, RMAX, VMAX
      // ..
      // .. Local Arrays ..
      int                LMAX( 4 );
      REAL               A( LDA, LDA ), AF( LDA, LDA ), B( LDB, LDB ), BF( LDB, LDB ), E( LDE, LDE ), F( LDF, LDF ), LSCALE( LDA ), RSCALE( LDA ), VL( LDVL, LDVL ), VLF( LDVL, LDVL ), VR( LDVR, LDVR ), VRF( LDVR, LDVR ), WORK( LDWORK, LDWORK )
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SGGBAK, SGGBAL, SLACPY
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

      EPS = SLAMCH( 'Precision' )

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

      ANORM = SLANGE( 'M', N, N, A, LDA, WORK )
      BNORM = SLANGE( 'M', N, N, B, LDB, WORK )

      slacpy('FULL', N, N, A, LDA, AF, LDA );
      slacpy('FULL', N, N, B, LDB, BF, LDB );

      sggbal('B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO );
      if ( INFO.NE.0 ) {
         NINFO = NINFO + 1
         LMAX( 1 ) = KNT
      }

      slacpy('FULL', N, M, VL, LDVL, VLF, LDVL );
      slacpy('FULL', N, M, VR, LDVR, VRF, LDVR );

      sggbak('B', 'L', N, ILO, IHI, LSCALE, RSCALE, M, VL, LDVL, INFO );
      if ( INFO.NE.0 ) {
         NINFO = NINFO + 1
         LMAX( 2 ) = KNT
      }

      sggbak('B', 'R', N, ILO, IHI, LSCALE, RSCALE, M, VR, LDVR, INFO );
      if ( INFO.NE.0 ) {
         NINFO = NINFO + 1
         LMAX( 3 ) = KNT
      }

      // Test of SGGBAK

      // Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
      // where tilde(A) denotes the transformed matrix.

      sgemm('N', 'N', N, M, N, ONE, AF, LDA, VR, LDVR, ZERO, WORK, LDWORK )       CALL SGEMM( 'T', 'N', M, M, N, ONE, VL, LDVL, WORK, LDWORK, ZERO, E, LDE );

      sgemm('N', 'N', N, M, N, ONE, A, LDA, VRF, LDVR, ZERO, WORK, LDWORK )       CALL SGEMM( 'T', 'N', M, M, N, ONE, VLF, LDVL, WORK, LDWORK, ZERO, F, LDF );

      VMAX = ZERO
      for (J = 1; J <= M; J++) { // 70
         for (I = 1; I <= M; I++) { // 60
            VMAX = MAX( VMAX, ABS( E( I, J )-F( I, J ) ) )
         } // 60
      } // 70
      VMAX = VMAX / ( EPS*MAX( ANORM, BNORM ) )
      if ( VMAX.GT.RMAX ) {
         LMAX( 4 ) = KNT
         RMAX = VMAX
      }

      // Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR

      sgemm('N', 'N', N, M, N, ONE, BF, LDB, VR, LDVR, ZERO, WORK, LDWORK )       CALL SGEMM( 'T', 'N', M, M, N, ONE, VL, LDVL, WORK, LDWORK, ZERO, E, LDE );

      sgemm('N', 'N', N, M, N, ONE, B, LDB, VRF, LDVR, ZERO, WORK, LDWORK )       CALL SGEMM( 'T', 'N', M, M, N, ONE, VLF, LDVL, WORK, LDWORK, ZERO, F, LDF );

      VMAX = ZERO
      for (J = 1; J <= M; J++) { // 90
         for (I = 1; I <= M; I++) { // 80
            VMAX = MAX( VMAX, ABS( E( I, J )-F( I, J ) ) )
         } // 80
      } // 90
      VMAX = VMAX / ( EPS*MAX( ANORM, BNORM ) )
      if ( VMAX.GT.RMAX ) {
         LMAX( 4 ) = KNT
         RMAX = VMAX
      }

      GO TO 10

      } // 100

      WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of SGGBAK .. ' )

      WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( ' value of largest test error                  =', E12.3 )
      WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( ' example number where SGGBAL info is not 0    =', I4 )
      WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( ' example number where SGGBAK(L) info is not 0 =', I4 )
      WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( ' example number where SGGBAK(R) info is not 0 =', I4 )
      WRITE( NOUT, FMT = 9994 )LMAX( 4 )
 9994 FORMAT( ' example number having largest error          =', I4 )
      WRITE( NOUT, FMT = 9992 )NINFO
 9992 FORMAT( ' number of examples where info is not 0       =', I4 )
      WRITE( NOUT, FMT = 9991 )KNT
 9991 FORMAT( ' total number of examples tested              =', I4 )

      RETURN

      // End of SCHKGK

      }
