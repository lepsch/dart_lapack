      SUBROUTINE DGET35( RMAX, LMAX, NINFO, KNT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, LMAX, NINFO;
      double             RMAX;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      double             TWO, FOUR;
      const              TWO = 2.0D0, FOUR = 4.0D0 ;
      // ..
      // .. Local Scalars ..
      String             TRANA, TRANB;
      int                I, IMA, IMB, IMLDA1, IMLDA2, IMLDB1, IMLOFF, INFO, ISGN, ITRANA, ITRANB, J, M, N;
      double             BIGNUM, CNRM, EPS, RES, RES1, RMUL, SCALE, SMLNUM, TNRM, XNRM;
      // ..
      // .. Local Arrays ..
      int                IDIM( 8 ), IVAL( 6, 6, 8 );
      double             A( 6, 6 ), B( 6, 6 ), C( 6, 6 ), CC( 6, 6 ), DUM( 1 ), VM1( 3 ), VM2( 3 );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DTRSYL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, SIN, SQRT
      // ..
      // .. Data statements ..
      DATA               IDIM / 1, 2, 3, 4, 3, 3, 6, 4 /
      DATA               IVAL / 1, 35*0, 1, 2, 4*0, -2, 0, 28*0, 1, 5*0, 5, 1, 2, 3*0, -8, -2, 1, 21*0, 3, 4, 4*0, -5, 3, 4*0, 1, 2, 1, 4, 2*0, -3, -9, -1, 1, 14*0, 1, 5*0, 2, 3, 4*0, 5, 6, 7, 21*0, 1, 5*0, 1, 3, -4, 3*0, 2, 5, 2, 21*0, 1, 2, 4*0, -2, 0, 4*0, 5, 6, 3, 4, 2*0, -1, -9, -5, 2, 2*0, 4*8, 5, 6, 4*9, -7, 5, 1, 5*0, 1, 5, 2, 3*0, 2, -21, 5, 3*0, 1, 2, 3, 4, 14*0 /
      // ..
      // .. Executable Statements ..

      // Get machine parameters

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )*FOUR / EPS
      BIGNUM = ONE / SMLNUM

      // Set up test case parameters

      VM1( 1 ) = SQRT( SMLNUM )
      VM1( 2 ) = ONE
      VM1( 3 ) = SQRT( BIGNUM )
      VM2( 1 ) = ONE
      VM2( 2 ) = ONE + TWO*EPS
      VM2( 3 ) = TWO

      KNT = 0
      NINFO = 0
      LMAX = 0
      RMAX = ZERO

      // Begin test loop

      for (ITRANA = 1; ITRANA <= 2; ITRANA++) { // 150
         for (ITRANB = 1; ITRANB <= 2; ITRANB++) { // 140
            DO 130 ISGN = -1, 1, 2
               for (IMA = 1; IMA <= 8; IMA++) { // 120
                  for (IMLDA1 = 1; IMLDA1 <= 3; IMLDA1++) { // 110
                     for (IMLDA2 = 1; IMLDA2 <= 3; IMLDA2++) { // 100
                        for (IMLOFF = 1; IMLOFF <= 2; IMLOFF++) { // 90
                           for (IMB = 1; IMB <= 8; IMB++) { // 80
                              for (IMLDB1 = 1; IMLDB1 <= 3; IMLDB1++) { // 70
                                 if (ITRANA == 1) TRANA = 'N'                                  IF( ITRANA == 2 ) TRANA = 'T'                                  IF( ITRANB == 1 ) TRANB = 'N'                                  IF( ITRANB == 2 ) TRANB = 'T';
                                 M = IDIM( IMA )
                                 N = IDIM( IMB )
                                 TNRM = ZERO
                                 for (I = 1; I <= M; I++) { // 20
                                    for (J = 1; J <= M; J++) { // 10
                                       A( I, J ) = IVAL( I, J, IMA )
                                       if ( ABS( I-J ).LE.1 ) {
                                          A( I, J ) = A( I, J )* VM1( IMLDA1 )                                           A( I, J ) = A( I, J )* VM2( IMLDA2 )
                                       } else {
                                          A( I, J ) = A( I, J )* VM1( IMLOFF )
                                       }
                                       TNRM = MAX( TNRM, ABS( A( I, J ) ) )
                                    } // 10
                                 } // 20
                                 for (I = 1; I <= N; I++) { // 40
                                    for (J = 1; J <= N; J++) { // 30
                                       B( I, J ) = IVAL( I, J, IMB )
                                       if ( ABS( I-J ).LE.1 ) {
                                          B( I, J ) = B( I, J )* VM1( IMLDB1 )
                                       } else {
                                          B( I, J ) = B( I, J )* VM1( IMLOFF )
                                       }
                                       TNRM = MAX( TNRM, ABS( B( I, J ) ) )
                                    } // 30
                                 } // 40
                                 CNRM = ZERO
                                 for (I = 1; I <= M; I++) { // 60
                                    for (J = 1; J <= N; J++) { // 50
                                       C( I, J ) = SIN( DBLE( I*J ) )
                                       CNRM = MAX( CNRM, C( I, J ) )
                                       CC( I, J ) = C( I, J )
                                    } // 50
                                 } // 60
                                 KNT = KNT + 1
                                 dtrsyl(TRANA, TRANB, ISGN, M, N, A, 6, B, 6, C, 6, SCALE, INFO );
                                 if (INFO != 0) NINFO = NINFO + 1;
                                 XNRM = DLANGE( 'M', M, N, C, 6, DUM )
                                 RMUL = ONE
                                 if ( XNRM > ONE && TNRM > ONE ) {
                                    if ( XNRM > BIGNUM / TNRM ) {
                                       RMUL = ONE / MAX( XNRM, TNRM )
                                    }
                                 }
                                 dgemm(TRANA, 'N', M, N, M, RMUL, A, 6, C, 6, -SCALE*RMUL, CC, 6 );
                                 dgemm('N', TRANB, M, N, N, DBLE( ISGN )*RMUL, C, 6, B, 6, ONE, CC, 6 );
                                 RES1 = DLANGE( 'M', M, N, CC, 6, DUM )
                                 RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM, ( ( RMUL*TNRM )*EPS )*XNRM )
                                 if ( RES > RMAX ) {
                                    LMAX = KNT
                                    RMAX = RES
                                 }
                              } // 70
                           } // 80
                        } // 90
                     } // 100
                  } // 110
               } // 120
            } // 130
         } // 140
      } // 150

      RETURN

      // End of DGET35

      }
