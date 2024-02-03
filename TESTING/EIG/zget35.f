      SUBROUTINE ZGET35( RMAX, LMAX, NINFO, KNT, NIN )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NIN, NINFO;
      double             RMAX;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                LDT;
      const              LDT = 10 ;
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 ;
      double             LARGE;
      const              LARGE = 1.0D6 ;
      COMPLEX*16         CONE
      const              CONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      String             TRANA, TRANB;
      int                I, IMLA, IMLAD, IMLB, IMLC, INFO, ISGN, ITRANA, ITRANB, J, M, N;
      double             BIGNUM, EPS, RES, RES1, SCALE, SMLNUM, TNRM, XNRM;
      COMPLEX*16         RMUL
      // ..
      // .. Local Arrays ..
      double             DUM( 1 ), VM1( 3 ), VM2( 3 );
      COMPLEX*16         A( LDT, LDT ), ATMP( LDT, LDT ), B( LDT, LDT ), BTMP( LDT, LDT ), C( LDT, LDT ), CSAV( LDT, LDT ), CTMP( LDT, LDT )
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZTRSYL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Get machine parameters

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // Set up test case parameters

      VM1( 1 ) = SQRT( SMLNUM )
      VM1( 2 ) = ONE
      VM1( 3 ) = LARGE
      VM2( 1 ) = ONE
      VM2( 2 ) = ONE + TWO*EPS
      VM2( 3 ) = TWO

      KNT = 0
      NINFO = 0
      LMAX = 0
      RMAX = ZERO

      // Begin test loop

      } // 10
      READ( NIN, FMT = * )M, N
      if (N == 0) RETURN;
      for (I = 1; I <= M; I++) { // 20
         READ( NIN, FMT = * )( ATMP( I, J ), J = 1, M )
      } // 20
      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( BTMP( I, J ), J = 1, N )
      } // 30
      for (I = 1; I <= M; I++) { // 40
         READ( NIN, FMT = * )( CTMP( I, J ), J = 1, N )
      } // 40
      for (IMLA = 1; IMLA <= 3; IMLA++) { // 170
         for (IMLAD = 1; IMLAD <= 3; IMLAD++) { // 160
            for (IMLB = 1; IMLB <= 3; IMLB++) { // 150
               for (IMLC = 1; IMLC <= 3; IMLC++) { // 140
                  for (ITRANA = 1; ITRANA <= 2; ITRANA++) { // 130
                     for (ITRANB = 1; ITRANB <= 2; ITRANB++) { // 120
                        DO 110 ISGN = -1, 1, 2
                           if (ITRANA == 1) TRANA = 'N'                            IF( ITRANA == 2 ) TRANA = 'C'                            IF( ITRANB == 1 ) TRANB = 'N'                            IF( ITRANB == 2 ) TRANB = 'C';
                           TNRM = ZERO
                           for (I = 1; I <= M; I++) { // 60
                              for (J = 1; J <= M; J++) { // 50
                                 A( I, J ) = ATMP( I, J )*VM1( IMLA )
                                 TNRM = MAX( TNRM, ABS( A( I, J ) ) )
                              } // 50
                              A( I, I ) = A( I, I )*VM2( IMLAD )
                              TNRM = MAX( TNRM, ABS( A( I, I ) ) )
                           } // 60
                           for (I = 1; I <= N; I++) { // 80
                              for (J = 1; J <= N; J++) { // 70
                                 B( I, J ) = BTMP( I, J )*VM1( IMLB )
                                 TNRM = MAX( TNRM, ABS( B( I, J ) ) )
                              } // 70
                           } // 80
                           if (TNRM == ZERO) TNRM = ONE;
                           for (I = 1; I <= M; I++) { // 100
                              for (J = 1; J <= N; J++) { // 90
                                 C( I, J ) = CTMP( I, J )*VM1( IMLC )
                                 CSAV( I, J ) = C( I, J )
                              } // 90
                           } // 100
                           KNT = KNT + 1
                           ztrsyl(TRANA, TRANB, ISGN, M, N, A, LDT, B, LDT, C, LDT, SCALE, INFO );
                           if (INFO.NE.0) NINFO = NINFO + 1;
                           XNRM = ZLANGE( 'M', M, N, C, LDT, DUM )
                           RMUL = CONE
                           if ( XNRM.GT.ONE .AND. TNRM.GT.ONE ) {
                              if ( XNRM.GT.BIGNUM / TNRM ) {
                                 RMUL = MAX( XNRM, TNRM )
                                 RMUL = CONE / RMUL
                              }
                           }
                           zgemm(TRANA, 'N', M, N, M, RMUL, A, LDT, C, LDT, -SCALE*RMUL, CSAV, LDT );
                           zgemm('N', TRANB, M, N, N, DBLE( ISGN )*RMUL, C, LDT, B, LDT, CONE, CSAV, LDT );
                           RES1 = ZLANGE( 'M', M, N, CSAV, LDT, DUM )
                           RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )
                           if ( RES.GT.RMAX ) {
                              LMAX = KNT
                              RMAX = RES
                           }
                        } // 110
                     } // 120
                  } // 130
               } // 140
            } // 150
         } // 160
      } // 170
      GO TO 10

      // End of ZGET35

      }
