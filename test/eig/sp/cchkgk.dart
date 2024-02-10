      void cchkgk(NIN, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                NIN, NOUT;
      // ..

      int                LDA, LDB, LDVL, LDVR;
      const              LDA = 50, LDB = 50, LDVL = 50, LDVR = 50 ;
      int                LDE, LDF, LDWORK, LRWORK;
      const              LDE = 50, LDF = 50, LDWORK = 50, LRWORK = 6*50 ;
      double               ZERO;
      const              ZERO = 0.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, IHI, ILO, INFO, J, KNT, M, N, NINFO;
      double               ANORM, BNORM, EPS, RMAX, VMAX;
      Complex            CDUM;
      int                LMAX( 4 );
      double               LSCALE( LDA ), RSCALE( LDA ), RWORK( LRWORK );
      Complex            A( LDA, LDA ), AF( LDA, LDA ), B( LDB, LDB ), BF( LDB, LDB ), E( LDE, LDE ), F( LDF, LDF ), VL( LDVL, LDVL ), VLF( LDVL, LDVL ), VR( LDVR, LDVR ), VRF( LDVR, LDVR ), WORK( LDWORK, LDWORK );
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGGBAK, CGGBAL, CLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[CDUM] = ( double( CDUM ) ).abs() + ( AIMAG( CDUM ) ).abs();

      LMAX[1] = 0;
      LMAX[2] = 0;
      LMAX[3] = 0;
      LMAX[4] = 0;
      NINFO = 0;
      KNT = 0;
      RMAX = ZERO;

      EPS = SLAMCH( 'Precision' );

      } // 10
      READ( NIN, FMT = * )N, M;
      if (N == 0) GO TO 100;

      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( A( I, J ), J = 1, N );
      } // 20

      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( B( I, J ), J = 1, N );
      } // 30

      for (I = 1; I <= N; I++) { // 40
         READ( NIN, FMT = * )( VL( I, J ), J = 1, M );
      } // 40

      for (I = 1; I <= N; I++) { // 50
         READ( NIN, FMT = * )( VR( I, J ), J = 1, M );
      } // 50

      KNT = KNT + 1;

      ANORM = CLANGE( 'M', N, N, A, LDA, RWORK );
      BNORM = CLANGE( 'M', N, N, B, LDB, RWORK );

      clacpy('FULL', N, N, A, LDA, AF, LDA );
      clacpy('FULL', N, N, B, LDB, BF, LDB );

      cggbal('B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, RWORK, INFO );
      if ( INFO != 0 ) {
         NINFO = NINFO + 1;
         LMAX[1] = KNT;
      }

      clacpy('FULL', N, M, VL, LDVL, VLF, LDVL );
      clacpy('FULL', N, M, VR, LDVR, VRF, LDVR );

      cggbak('B', 'L', N, ILO, IHI, LSCALE, RSCALE, M, VL, LDVL, INFO );
      if ( INFO != 0 ) {
         NINFO = NINFO + 1;
         LMAX[2] = KNT;
      }

      cggbak('B', 'R', N, ILO, IHI, LSCALE, RSCALE, M, VR, LDVR, INFO );
      if ( INFO != 0 ) {
         NINFO = NINFO + 1;
         LMAX[3] = KNT;
      }

      // Test of CGGBAK

      // Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
      // where tilde(A) denotes the transformed matrix.

      cgemm('N', 'N', N, M, N, CONE, AF, LDA, VR, LDVR, CZERO, WORK, LDWORK )       CALL CGEMM( 'C', 'N', M, M, N, CONE, VL, LDVL, WORK, LDWORK, CZERO, E, LDE );

      cgemm('N', 'N', N, M, N, CONE, A, LDA, VRF, LDVR, CZERO, WORK, LDWORK )       CALL CGEMM( 'C', 'N', M, M, N, CONE, VLF, LDVL, WORK, LDWORK, CZERO, F, LDF );

      VMAX = ZERO;
      for (J = 1; J <= M; J++) { // 70
         for (I = 1; I <= M; I++) { // 60
            VMAX = max( VMAX, CABS1( E( I, J )-F( I, J ) ) );
         } // 60
      } // 70
      VMAX = VMAX / ( EPS*max( ANORM, BNORM ) );
      if ( VMAX > RMAX ) {
         LMAX[4] = KNT;
         RMAX = VMAX;
      }

      // Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR

      cgemm('N', 'N', N, M, N, CONE, BF, LDB, VR, LDVR, CZERO, WORK, LDWORK )       CALL CGEMM( 'C', 'N', M, M, N, CONE, VL, LDVL, WORK, LDWORK, CZERO, E, LDE );

      cgemm('n', 'n', N, M, N, CONE, B, LDB, VRF, LDVR, CZERO, WORK, LDWORK )       CALL CGEMM( 'C', 'N', M, M, N, CONE, VLF, LDVL, WORK, LDWORK, CZERO, F, LDF );

      VMAX = ZERO;
      for (J = 1; J <= M; J++) { // 90
         for (I = 1; I <= M; I++) { // 80
            VMAX = max( VMAX, CABS1( E( I, J )-F( I, J ) ) );
         } // 80
      } // 90
      VMAX = VMAX / ( EPS*max( ANORM, BNORM ) );
      if ( VMAX > RMAX ) {
         LMAX[4] = KNT;
         RMAX = VMAX;
      }

      GO TO 10;

      } // 100

      WRITE( NOUT, FMT = 9999 );
 9999 FORMAT(' .. test output of CGGBAK .. ' );

      WRITE( NOUT, FMT = 9998 )RMAX;
 9998 FORMAT( ' value of largest test error                  =', E12.3 );
      WRITE( NOUT, FMT = 9997 )LMAX( 1 );
 9997 FORMAT( ' example number where CGGBAL info is not 0    =${.i4}');
      WRITE( NOUT, FMT = 9996 )LMAX( 2 );
 9996 FORMAT( ' example number where CGGBAK(L) info is not 0 =${.i4}');
      WRITE( NOUT, FMT = 9995 )LMAX( 3 );
 9995 FORMAT( ' example number where CGGBAK(R) info is not 0 =${.i4}');
      WRITE( NOUT, FMT = 9994 )LMAX( 4 );
 9994 FORMAT( ' example number having largest error          =${.i4}');
      WRITE( NOUT, FMT = 9992 )NINFO;
 9992 FORMAT( ' number of examples where info is not 0       =${.i4}');
      WRITE( NOUT, FMT = 9991 )KNT;
 9991 FORMAT( ' total number of examples tested              =${.i4}');

      }
