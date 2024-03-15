      void sget23(final int COMP, final int BALANC, final int JTYPE, final int THRESH, final Array<int> ISEED_, final int NOUNIT, final int N, final Matrix<double> A_, final int LDA, final int H, final int WR, final int WI, final int WR1, final int WI1, final Matrix<double> VL_, final int LDVL, final Matrix<double> VR_, final int LDVR, final Matrix<double> LRE_, final int LDLRE, final int RCONDV, final int RCNDV1, final int RCDVIN, final int RCONDE, final int RCNDE1, final int RCDEIN, final int SCALE, final int SCALE1, final int RESULT, final Array<double> WORK_, final int LWORK, final Array<int> IWORK_, final Box<int> INFO,) {
  final ISEED = ISEED_.dim();
  final A = A_.dim();
  final VL = VL_.dim();
  final VR = VR_.dim();
  final LRE = LRE_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               COMP;
      String             BALANC;
      int                INFO, JTYPE, LDA, LDLRE, LDVL, LDVR, LWORK, N, NOUNIT;
      double               THRESH;
      int                ISEED( 4 ), IWORK( * );
      double               A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), RCDEIN( * ), RCDVIN( * ), RCNDE1( * ), RCNDV1( * ), RCONDE( * ), RCONDV( * ), RESULT( 11 ), SCALE( * ), SCALE1( * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WI1( * ), WORK( * ), WR( * ), WR1( * );
      // ..

// =====================================================================


      // .. Parameters ..
      double               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double               EPSIN;
      const              EPSIN = 5.9605e-8 ;
      bool               BALOK, NOBAL;
      String             SENSE;
      int                I, IHI, IHI1, IINFO, ILO, ILO1, ISENS, ISENSM, J, JJ, KMIN;
      double               ABNRM, ABNRM1, EPS, SMLNUM, TNRM, TOL, TOLIN, ULP, ULPINV, V, VIMIN, VMAX, VMX, VRMIN, VRMX, VTST;
      String             SENS( 2 );
      double               DUM( 1 ), RES( 2 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLAPY2, SNRM2;
      // EXTERNAL lsame, SLAMCH, SLAPY2, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEEVX, SGET22, SLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL
      // ..
      // .. Data statements ..
      const SENS = [ 'N', 'V' ];

      // Check for errors

      NOBAL = lsame( BALANC, 'N' );
      BALOK = NOBAL || lsame( BALANC, 'P' ) || lsame( BALANC, 'S' ) || lsame( BALANC, 'B' );
      INFO = 0;
      if ( !BALOK ) {
         INFO = -2;
      } else if ( THRESH < ZERO ) {
         INFO = -4;
      } else if ( NOUNIT <= 0 ) {
         INFO = -6;
      } else if ( N < 0 ) {
         INFO = -7;
      } else if ( LDA < 1 || LDA < N ) {
         INFO = -9;
      } else if ( LDVL < 1 || LDVL < N ) {
         INFO = -16;
      } else if ( LDVR < 1 || LDVR < N ) {
         INFO = -18;
      } else if ( LDLRE < 1 || LDLRE < N ) {
         INFO = -20;
      } else if ( LWORK < 3*N || ( COMP && LWORK < 6*N+N*N ) ) {
         INFO = -31;
      }

      if ( INFO != 0 ) {
         xerbla('SGET23', -INFO );
         return;
      }

      // Quick return if nothing to do

      for (I = 1; I <= 11; I++) { // 10
         RESULT[I] = -ONE;
      } // 10

      if (N == 0) return;

      // More Important constants

      ULP = SLAMCH( 'Precision' );
      SMLNUM = SLAMCH( 'S' );
      ULPINV = ONE / ULP;

      // Compute eigenvalues and eigenvectors, and test them

      if ( LWORK >= 6*N+N*N ) {
         SENSE = 'B';
         ISENSM = 2;
      } else {
         SENSE = 'E';
         ISENSM = 1;
      }
      slacpy('F', N, N, A, LDA, H, LDA );
      sgeevx(BALANC, 'V', 'V', SENSE, N, H, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, IINFO );
      if ( IINFO != 0 ) {
         RESULT[1] = ULPINV;
         if ( JTYPE != 22 ) {
            WRITE( NOUNIT, FMT = 9998 )'SGEEVX1', IINFO, N, JTYPE, BALANC, ISEED;
         } else {
            WRITE( NOUNIT, FMT = 9999 )'SGEEVX1', IINFO, N, ISEED( 1 );
         }
         INFO = ( IINFO ).abs();
         return;
      }

      // Do Test (1)

      sget22('N', 'N', 'N', N, A, LDA, VR, LDVR, WR, WI, WORK, RES );
      RESULT[1] = RES( 1 );

      // Do Test (2)

      sget22('T', 'N', 'T', N, A, LDA, VL, LDVL, WR, WI, WORK, RES );
      RESULT[2] = RES( 1 );

      // Do Test (3)

      for (J = 1; J <= N; J++) { // 30
         TNRM = ONE;
         if ( WI( J ) == ZERO ) {
            TNRM = SNRM2( N, VR( 1, J ), 1 );
         } else if ( WI( J ) > ZERO ) {
            TNRM = SLAPY2( SNRM2( N, VR( 1, J ), 1 ), SNRM2( N, VR( 1, J+1 ), 1 ) );
         }
         RESULT[3] = max( RESULT( 3 ), min( ULPINV, ( TNRM-ONE ).abs() / ULP ) );
         if ( WI( J ) > ZERO ) {
            VMX = ZERO;
            VRMX = ZERO;
            for (JJ = 1; JJ <= N; JJ++) { // 20
               VTST = SLAPY2( VR( JJ, J ), VR( JJ, J+1 ) );
               if (VTST > VMX) VMX = VTST;
               IF( VR( JJ, J+1 ) == ZERO && ( VR( JJ, J ) ).abs() > VRMX )VRMX = ( VR( JJ, J ) ).abs();
            } // 20
            if (VRMX / VMX < ONE-TWO*ULP) RESULT( 3 ) = ULPINV;
         }
      } // 30

      // Do Test (4)

      for (J = 1; J <= N; J++) { // 50
         TNRM = ONE;
         if ( WI( J ) == ZERO ) {
            TNRM = SNRM2( N, VL( 1, J ), 1 );
         } else if ( WI( J ) > ZERO ) {
            TNRM = SLAPY2( SNRM2( N, VL( 1, J ), 1 ), SNRM2( N, VL( 1, J+1 ), 1 ) );
         }
         RESULT[4] = max( RESULT( 4 ), min( ULPINV, ( TNRM-ONE ).abs() / ULP ) );
         if ( WI( J ) > ZERO ) {
            VMX = ZERO;
            VRMX = ZERO;
            for (JJ = 1; JJ <= N; JJ++) { // 40
               VTST = SLAPY2( VL( JJ, J ), VL( JJ, J+1 ) );
               if (VTST > VMX) VMX = VTST;
               IF( VL( JJ, J+1 ) == ZERO && ( VL( JJ, J ) ).abs() > VRMX )VRMX = ( VL( JJ, J ) ).abs();
            } // 40
            if (VRMX / VMX < ONE-TWO*ULP) RESULT( 4 ) = ULPINV;
         }
      } // 50

      // Test for all options of computing condition numbers

      for (ISENS = 1; ISENS <= ISENSM; ISENS++) { // 200

         SENSE = SENS( ISENS );

         // Compute eigenvalues only, and test them

         slacpy('F', N, N, A, LDA, H, LDA );
         sgeevx(BALANC, 'N', 'N', SENSE, N, H, LDA, WR1, WI1, DUM, 1, DUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, IWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEEVX2', IINFO, N, JTYPE, BALANC, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEEVX2', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 190;
         }

         // Do Test (5)

         for (J = 1; J <= N; J++) { // 60
            if( WR( J ) != WR1( J ) || WI( J ) != WI1( J ) ) RESULT( 5 ) = ULPINV;
         } // 60

         // Do Test (8)

         if ( !NOBAL ) {
            for (J = 1; J <= N; J++) { // 70
               if( SCALE( J ) != SCALE1( J ) ) RESULT( 8 ) = ULPINV;
            } // 70
            if (ILO != ILO1) RESULT( 8 ) = ULPINV;
            if( IHI != IHI1 ) RESULT( 8 ) = ULPINV;
            IF( ABNRM != ABNRM1 ) RESULT( 8 ) = ULPINV;
         }

         // Do Test (9)

         if ( ISENS == 2 && N > 1 ) {
            for (J = 1; J <= N; J++) { // 80
               if( RCONDV( J ) != RCNDV1( J ) ) RESULT( 9 ) = ULPINV;
            } // 80
         }

         // Compute eigenvalues and right eigenvectors, and test them

         slacpy('F', N, N, A, LDA, H, LDA );
         sgeevx(BALANC, 'N', 'V', SENSE, N, H, LDA, WR1, WI1, DUM, 1, LRE, LDLRE, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, IWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEEVX3', IINFO, N, JTYPE, BALANC, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEEVX3', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 190;
         }

         // Do Test (5) again

         for (J = 1; J <= N; J++) { // 90
            if( WR( J ) != WR1( J ) || WI( J ) != WI1( J ) ) RESULT( 5 ) = ULPINV;
         } // 90

         // Do Test (6)

         for (J = 1; J <= N; J++) { // 110
            for (JJ = 1; JJ <= N; JJ++) { // 100
               if( VR( J, JJ ) != LRE( J, JJ ) ) RESULT( 6 ) = ULPINV;
            } // 100
         } // 110

         // Do Test (8) again

         if ( !NOBAL ) {
            for (J = 1; J <= N; J++) { // 120
               if( SCALE( J ) != SCALE1( J ) ) RESULT( 8 ) = ULPINV;
            } // 120
            if (ILO != ILO1) RESULT( 8 ) = ULPINV;
            if( IHI != IHI1 ) RESULT( 8 ) = ULPINV;
            IF( ABNRM != ABNRM1 ) RESULT( 8 ) = ULPINV;
         }

         // Do Test (9) again

         if ( ISENS == 2 && N > 1 ) {
            for (J = 1; J <= N; J++) { // 130
               if( RCONDV( J ) != RCNDV1( J ) ) RESULT( 9 ) = ULPINV;
            } // 130
         }

         // Compute eigenvalues and left eigenvectors, and test them

         slacpy('F', N, N, A, LDA, H, LDA );
         sgeevx(BALANC, 'V', 'N', SENSE, N, H, LDA, WR1, WI1, LRE, LDLRE, DUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, IWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'SGEEVX4', IINFO, N, JTYPE, BALANC, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'SGEEVX4', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 190;
         }

         // Do Test (5) again

         for (J = 1; J <= N; J++) { // 140
            if( WR( J ) != WR1( J ) || WI( J ) != WI1( J ) ) RESULT( 5 ) = ULPINV;
         } // 140

         // Do Test (7)

         for (J = 1; J <= N; J++) { // 160
            for (JJ = 1; JJ <= N; JJ++) { // 150
               if( VL( J, JJ ) != LRE( J, JJ ) ) RESULT( 7 ) = ULPINV;
            } // 150
         } // 160

         // Do Test (8) again

         if ( !NOBAL ) {
            for (J = 1; J <= N; J++) { // 170
               if( SCALE( J ) != SCALE1( J ) ) RESULT( 8 ) = ULPINV;
            } // 170
            if (ILO != ILO1) RESULT( 8 ) = ULPINV;
            if( IHI != IHI1 ) RESULT( 8 ) = ULPINV;
            IF( ABNRM != ABNRM1 ) RESULT( 8 ) = ULPINV;
         }

         // Do Test (9) again

         if ( ISENS == 2 && N > 1 ) {
            for (J = 1; J <= N; J++) { // 180
               if( RCONDV( J ) != RCNDV1( J ) ) RESULT( 9 ) = ULPINV;
            } // 180
         }

         } // 190

      } // 200

      // If COMP, compare condition numbers to precomputed ones

      if ( COMP ) {
         slacpy('F', N, N, A, LDA, H, LDA );
         sgeevx('N', 'V', 'V', 'B', N, H, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = ULPINV;
            WRITE( NOUNIT, FMT = 9999 )'SGEEVX5', IINFO, N, ISEED( 1 );
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         // Sort eigenvalues and condition numbers lexicographically
         // to compare with inputs

         for (I = 1; I <= N - 1; I++) { // 220
            KMIN = I;
            VRMIN = WR( I );
            VIMIN = WI( I );
            for (J = I + 1; J <= N; J++) { // 210
               if ( WR( J ) < VRMIN ) {
                  KMIN = J;
                  VRMIN = WR( J );
                  VIMIN = WI( J );
               }
            } // 210
            WR[KMIN] = WR( I );
            WI[KMIN] = WI( I );
            WR[I] = VRMIN;
            WI[I] = VIMIN;
            VRMIN = RCONDE( KMIN );
            RCONDE[KMIN] = RCONDE( I );
            RCONDE[I] = VRMIN;
            VRMIN = RCONDV( KMIN );
            RCONDV[KMIN] = RCONDV( I );
            RCONDV[I] = VRMIN;
         } // 220

         // Compare condition numbers for eigenvectors
         // taking their condition numbers into account

         RESULT[10] = ZERO;
         EPS = max( EPSIN, ULP );
         V = max( double( N )*EPS*ABNRM, SMLNUM );
         if (ABNRM == ZERO) V = ONE;
         for (I = 1; I <= N; I++) { // 230
            if ( V > RCONDV( I )*RCONDE( I ) ) {
               TOL = RCONDV( I );
            } else {
               TOL = V / RCONDE( I );
            }
            if ( V > RCDVIN( I )*RCDEIN( I ) ) {
               TOLIN = RCDVIN( I );
            } else {
               TOLIN = V / RCDEIN( I );
            }
            TOL = max( TOL, SMLNUM / EPS );
            TOLIN = max( TOLIN, SMLNUM / EPS );
            if ( EPS*( RCDVIN( I )-TOLIN ) > RCONDV( I )+TOL ) {
               VMAX = ONE / EPS;
            } else if ( RCDVIN( I )-TOLIN > RCONDV( I )+TOL ) {
               VMAX = ( RCDVIN( I )-TOLIN ) / ( RCONDV( I )+TOL );
            } else if ( RCDVIN( I )+TOLIN < EPS*( RCONDV( I )-TOL ) ) {
               VMAX = ONE / EPS;
            } else if ( RCDVIN( I )+TOLIN < RCONDV( I )-TOL ) {
               VMAX = ( RCONDV( I )-TOL ) / ( RCDVIN( I )+TOLIN );
            } else {
               VMAX = ONE;
            }
            RESULT[10] = max( RESULT( 10 ), VMAX );
         } // 230

         // Compare condition numbers for eigenvalues
         // taking their condition numbers into account

         RESULT[11] = ZERO;
         for (I = 1; I <= N; I++) { // 240
            if ( V > RCONDV( I ) ) {
               TOL = ONE;
            } else {
               TOL = V / RCONDV( I );
            }
            if ( V > RCDVIN( I ) ) {
               TOLIN = ONE;
            } else {
               TOLIN = V / RCDVIN( I );
            }
            TOL = max( TOL, SMLNUM / EPS );
            TOLIN = max( TOLIN, SMLNUM / EPS );
            if ( EPS*( RCDEIN( I )-TOLIN ) > RCONDE( I )+TOL ) {
               VMAX = ONE / EPS;
            } else if ( RCDEIN( I )-TOLIN > RCONDE( I )+TOL ) {
               VMAX = ( RCDEIN( I )-TOLIN ) / ( RCONDE( I )+TOL );
            } else if ( RCDEIN( I )+TOLIN < EPS*( RCONDE( I )-TOL ) ) {
               VMAX = ONE / EPS;
            } else if ( RCDEIN( I )+TOLIN < RCONDE( I )-TOL ) {
               VMAX = ( RCONDE( I )-TOL ) / ( RCDEIN( I )+TOLIN );
            } else {
               VMAX = ONE;
            }
            RESULT[11] = max( RESULT( 11 ), VMAX );
         } // 240
         } // 250

      }

 9999 FORMAT( ' SGET23: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, INPUT EXAMPLE NUMBER = ${.i4}');
 9998 FORMAT( ' SGET23: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, BALANC = ${}, ISEED=(${.i5(4, ',')})' );

      }