      SUBROUTINE ZCHKGK( NIN, NOUT )

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
      int                LDE, LDF, LDWORK, LRWORK;
      const              LDE = 50, LDF = 50, LDWORK = 50, LRWORK = 6*50 ;
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, IHI, ILO, INFO, J, KNT, M, N, NINFO;
      double             ANORM, BNORM, EPS, RMAX, VMAX;
      COMPLEX*16         CDUM
      // ..
      // .. Local Arrays ..
      int                LMAX( 4 );
      double             LSCALE( LDA ), RSCALE( LDA ), RWORK( LRWORK );
      COMPLEX*16         A( LDA, LDA ), AF( LDA, LDA ), B( LDB, LDB ), BF( LDB, LDB ), E( LDE, LDE ), F( LDF, LDF ), VL( LDVL, LDVL ), VLF( LDVL, LDVL ), VR( LDVR, LDVR ), VRF( LDVR, LDVR ), WORK( LDWORK, LDWORK )
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZGGBAK, ZGGBAL, ZLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      LMAX( 1 ) = 0
      LMAX( 2 ) = 0
      LMAX( 3 ) = 0
      LMAX( 4 ) = 0
      NINFO = 0
      KNT = 0
      RMAX = ZERO

      EPS = DLAMCH( 'Precision' )

   10 CONTINUE
      READ( NIN, FMT = * )N, M
      IF( N.EQ.0 ) GO TO 100

      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( A( I, J ), J = 1, N )
   20 CONTINUE

      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( B( I, J ), J = 1, N )
   30 CONTINUE

      for (I = 1; I <= N; I++) { // 40
         READ( NIN, FMT = * )( VL( I, J ), J = 1, M )
   40 CONTINUE

      for (I = 1; I <= N; I++) { // 50
         READ( NIN, FMT = * )( VR( I, J ), J = 1, M )
   50 CONTINUE

      KNT = KNT + 1

      ANORM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      BNORM = ZLANGE( 'M', N, N, B, LDB, RWORK )

      zlacpy('FULL', N, N, A, LDA, AF, LDA );
      zlacpy('FULL', N, N, B, LDB, BF, LDB );

      zggbal('B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, RWORK, INFO );
      if ( INFO.NE.0 ) {
         NINFO = NINFO + 1
         LMAX( 1 ) = KNT
      }

      zlacpy('FULL', N, M, VL, LDVL, VLF, LDVL );
      zlacpy('FULL', N, M, VR, LDVR, VRF, LDVR );

      zggbak('B', 'L', N, ILO, IHI, LSCALE, RSCALE, M, VL, LDVL, INFO );
      if ( INFO.NE.0 ) {
         NINFO = NINFO + 1
         LMAX( 2 ) = KNT
      }

      zggbak('B', 'R', N, ILO, IHI, LSCALE, RSCALE, M, VR, LDVR, INFO );
      if ( INFO.NE.0 ) {
         NINFO = NINFO + 1
         LMAX( 3 ) = KNT
      }

      // Test of ZGGBAK

      // Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
      // where tilde(A) denotes the transformed matrix.

      zgemm('N', 'N', N, M, N, CONE, AF, LDA, VR, LDVR, CZERO, WORK, LDWORK )       CALL ZGEMM( 'C', 'N', M, M, N, CONE, VL, LDVL, WORK, LDWORK, CZERO, E, LDE );

      zgemm('N', 'N', N, M, N, CONE, A, LDA, VRF, LDVR, CZERO, WORK, LDWORK )       CALL ZGEMM( 'C', 'N', M, M, N, CONE, VLF, LDVL, WORK, LDWORK, CZERO, F, LDF );

      VMAX = ZERO
      for (J = 1; J <= M; J++) { // 70
         for (I = 1; I <= M; I++) { // 60
            VMAX = MAX( VMAX, CABS1( E( I, J )-F( I, J ) ) )
   60    CONTINUE
   70 CONTINUE
      VMAX = VMAX / ( EPS*MAX( ANORM, BNORM ) )
      if ( VMAX.GT.RMAX ) {
         LMAX( 4 ) = KNT
         RMAX = VMAX
      }

      // Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR

      zgemm('N', 'N', N, M, N, CONE, BF, LDB, VR, LDVR, CZERO, WORK, LDWORK )       CALL ZGEMM( 'C', 'N', M, M, N, CONE, VL, LDVL, WORK, LDWORK, CZERO, E, LDE );

      zgemm('n', 'n', N, M, N, CONE, B, LDB, VRF, LDVR, CZERO, WORK, LDWORK )       CALL ZGEMM( 'C', 'N', M, M, N, CONE, VLF, LDVL, WORK, LDWORK, CZERO, F, LDF );

      VMAX = ZERO
      for (J = 1; J <= M; J++) { // 90
         for (I = 1; I <= M; I++) { // 80
            VMAX = MAX( VMAX, CABS1( E( I, J )-F( I, J ) ) )
   80    CONTINUE
   90 CONTINUE
      VMAX = VMAX / ( EPS*MAX( ANORM, BNORM ) )
      if ( VMAX.GT.RMAX ) {
         LMAX( 4 ) = KNT
         RMAX = VMAX
      }

      GO TO 10

  100 CONTINUE

      WRITE( NOUT, FMT = 9999 )
 9999 FORMAT( 1X, '.. test output of ZGGBAK .. ' )

      WRITE( NOUT, FMT = 9998 )RMAX
 9998 FORMAT( ' value of largest test error                  =', D12.3 )
      WRITE( NOUT, FMT = 9997 )LMAX( 1 )
 9997 FORMAT( ' example number where ZGGBAL info is not 0    =', I4 )
      WRITE( NOUT, FMT = 9996 )LMAX( 2 )
 9996 FORMAT( ' example number where ZGGBAK(L) info is not 0 =', I4 )
      WRITE( NOUT, FMT = 9995 )LMAX( 3 )
 9995 FORMAT( ' example number where ZGGBAK(R) info is not 0 =', I4 )
      WRITE( NOUT, FMT = 9994 )LMAX( 4 )
 9994 FORMAT( ' example number having largest error          =', I4 )
      WRITE( NOUT, FMT = 9992 )NINFO
 9992 FORMAT( ' number of examples where info is not 0       =', I4 )
      WRITE( NOUT, FMT = 9991 )KNT
 9991 FORMAT( ' total number of examples tested              =', I4 )

      RETURN

      // End of ZCHKGK

      }
