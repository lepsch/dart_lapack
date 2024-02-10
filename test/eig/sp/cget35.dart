      void cget35(final int RMAX, final int LMAX, final int NINFO, final int KNT, final int NIN) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KNT, LMAX, NIN, NINFO;
      double               RMAX;
      // ..

      int                LDT;
      const              LDT = 10 ;
      double               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double               LARGE;
      const              LARGE = 1.0e6 ;
      Complex            CONE;
      const              CONE = 1.0 ;
      String             TRANA, TRANB;
      int                I, IMLA, IMLAD, IMLB, IMLC, INFO, ISGN, ITRANA, ITRANB, J, M, N;
      double               BIGNUM, EPS, RES, RES1, SCALE, SMLNUM, TNRM, XNRM;
      Complex            RMUL;
      double               DUM( 1 ), VM1( 3 ), VM2( 3 );
      Complex            A( LDT, LDT ), ATMP( LDT, LDT ), B( LDT, LDT ), BTMP( LDT, LDT ), C( LDT, LDT ), CSAV( LDT, LDT ), CTMP( LDT, LDT );
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CTRSYL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, SQRT

      // Get machine parameters

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Set up test case parameters

      VM1[1] = sqrt( SMLNUM );
      VM1[2] = ONE;
      VM1[3] = LARGE;
      VM2[1] = ONE;
      VM2[2] = ONE + TWO*EPS;
      VM2[3] = TWO;

      KNT = 0;
      NINFO = 0;
      LMAX = 0;
      RMAX = ZERO;

      // Begin test loop

      } // 10
      READ( NIN, FMT = * )M, N;
      if (N == 0) return;
      for (I = 1; I <= M; I++) { // 20
         READ( NIN, FMT = * )( ATMP( I, J ), J = 1, M );
      } // 20
      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )( BTMP( I, J ), J = 1, N );
      } // 30
      for (I = 1; I <= M; I++) { // 40
         READ( NIN, FMT = * )( CTMP( I, J ), J = 1, N );
      } // 40
      for (IMLA = 1; IMLA <= 3; IMLA++) { // 170
         for (IMLAD = 1; IMLAD <= 3; IMLAD++) { // 160
            for (IMLB = 1; IMLB <= 3; IMLB++) { // 150
               for (IMLC = 1; IMLC <= 3; IMLC++) { // 140
                  for (ITRANA = 1; ITRANA <= 2; ITRANA++) { // 130
                     for (ITRANB = 1; ITRANB <= 2; ITRANB++) { // 120
                        for (ISGN = -1; 2 < 0 ? ISGN >= 1 : ISGN <= 1; ISGN += 2) { // 110
                           if (ITRANA == 1) TRANA = 'N';
                           if( ITRANA == 2 ) TRANA = 'C';
                           if( ITRANB == 1 ) TRANB = 'N';
                           IF( ITRANB == 2 ) TRANB = 'C';
                           TNRM = ZERO;
                           for (I = 1; I <= M; I++) { // 60
                              for (J = 1; J <= M; J++) { // 50
                                 A[I][J] = ATMP( I, J )*VM1( IMLA );
                                 TNRM = max( TNRM, ( A( I, J ) ).abs() );
                              } // 50
                              A[I][I] = A( I, I )*VM2( IMLAD );
                              TNRM = max( TNRM, ( A( I, I ) ).abs() );
                           } // 60
                           for (I = 1; I <= N; I++) { // 80
                              for (J = 1; J <= N; J++) { // 70
                                 B[I][J] = BTMP( I, J )*VM1( IMLB );
                                 TNRM = max( TNRM, ( B( I, J ) ).abs() );
                              } // 70
                           } // 80
                           if (TNRM == ZERO) TNRM = ONE;
                           for (I = 1; I <= M; I++) { // 100
                              for (J = 1; J <= N; J++) { // 90
                                 C[I][J] = CTMP( I, J )*VM1( IMLC );
                                 CSAV[I][J] = C( I, J );
                              } // 90
                           } // 100
                           KNT = KNT + 1;
                           ctrsyl(TRANA, TRANB, ISGN, M, N, A, LDT, B, LDT, C, LDT, SCALE, INFO );
                           if (INFO != 0) NINFO = NINFO + 1;
                           XNRM = CLANGE( 'M', M, N, C, LDT, DUM );
                           RMUL = CONE;
                           if ( XNRM > ONE && TNRM > ONE ) {
                              if ( XNRM > BIGNUM / TNRM ) {
                                 RMUL = max( XNRM, TNRM );
                                 RMUL = CONE / RMUL;
                              }
                           }
                           cgemm(TRANA, 'N', M, N, M, RMUL, A, LDT, C, LDT, -SCALE*RMUL, CSAV, LDT );
                           cgemm('N', TRANB, M, N, N, double( ISGN )*RMUL, C, LDT, B, LDT, CONE, CSAV, LDT );
                           RES1 = CLANGE( 'M', M, N, CSAV, LDT, DUM );
                           RES = RES1 / max( SMLNUM, SMLNUM*XNRM, ( ( ( RMUL ).abs()*TNRM )*EPS )*XNRM );
                           if ( RES > RMAX ) {
                              LMAX = KNT;
                              RMAX = RES;
                           }
                        } // 110
                     } // 120
                  } // 130
               } // 140
            } // 150
         } // 160
      } // 170
      GO TO 10;
      }
