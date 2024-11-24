// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cget23(final int COMP, final int ISRT, final int BALANC, final int JTYPE, final int THRESH, final Array<int> ISEED_, final int NOUNIT, final int N, final Matrix<double> A_, final int LDA, final int H, final int W, final int W1, final Matrix<double> VL_, final int LDVL, final Matrix<double> VR_, final int LDVR, final Matrix<double> LRE_, final int LDLRE, final int RCONDV, final int RCNDV1, final int RCDVIN, final int RCONDE, final int RCNDE1, final int RCDEIN, final int SCALE, final int SCALE1, final int RESULT, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final Box<int> INFO,) {
  final ISEED = ISEED_.dim();
  final A = A_.dim();
  final VL = VL_.dim();
  final VR = VR_.dim();
  final LRE = LRE_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               COMP;
      String             BALANC;
      int                INFO, ISRT, JTYPE, LDA, LDLRE, LDVL, LDVR, LWORK, N, NOUNIT;
      double               THRESH;
      int                ISEED( 4 );
      double               RCDEIN( * ), RCDVIN( * ), RCNDE1( * ), RCNDV1( * ), RCONDE( * ), RCONDV( * ), RESULT( 11 ), RWORK( * ), SCALE( * ), SCALE1( * );
      Complex            A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), W1( * ), WORK( * );
      // ..

      double               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double               EPSIN;
      const              EPSIN = 5.9605e-8 ;
      bool               BALOK, NOBAL;
      String             SENSE;
      int                I, IHI, IHI1, IINFO, ILO, ILO1, ISENS, ISENSM, J, JJ, KMIN;
      double               ABNRM, ABNRM1, EPS, SMLNUM, TNRM, TOL, TOLIN, ULP, ULPINV, V, VMAX, VMX, VRICMP, VRIMIN, VRMX, VTST;
      Complex            CTMP;
      String             SENS( 2 );
      double               RES( 2 );
      Complex            CDUM( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SCNRM2, SLAMCH;
      // EXTERNAL lsame, SCNRM2, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEEVX, CGET22, CLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Data statements ..
      const SENS = [ 'N', 'V' ];

      // Check for errors

      NOBAL = lsame( BALANC, 'N' );
      BALOK = NOBAL || lsame( BALANC, 'P' ) || lsame( BALANC, 'S' ) || lsame( BALANC, 'B' );
      INFO = 0;
      if ( ISRT != 0 && ISRT != 1 ) {
         INFO = -2;
      } else if ( !BALOK ) {
         INFO = -3;
      } else if ( THRESH < ZERO ) {
         INFO = -5;
      } else if ( NOUNIT <= 0 ) {
         INFO = -7;
      } else if ( N < 0 ) {
         INFO = -8;
      } else if ( LDA < 1 || LDA < N ) {
         INFO = -10;
      } else if ( LDVL < 1 || LDVL < N ) {
         INFO = -15;
      } else if ( LDVR < 1 || LDVR < N ) {
         INFO = -17;
      } else if ( LDLRE < 1 || LDLRE < N ) {
         INFO = -19;
      } else if ( LWORK < 2*N || ( COMP && LWORK < 2*N+N*N ) ) {
         INFO = -30;
      }

      if ( INFO != 0 ) {
         xerbla('CGET23', -INFO );
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

      if ( LWORK >= 2*N+N*N ) {
         SENSE = 'B';
         ISENSM = 2;
      } else {
         SENSE = 'E';
         ISENSM = 1;
      }
      clacpy('F', N, N, A, LDA, H, LDA );
      cgeevx(BALANC, 'V', 'V', SENSE, N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO );
      if ( IINFO != 0 ) {
         RESULT[1] = ULPINV;
         if ( JTYPE != 22 ) {
            WRITE( NOUNIT, FMT = 9998 )'CGEEVX1', IINFO, N, JTYPE, BALANC, ISEED;
         } else {
            WRITE( NOUNIT, FMT = 9999 )'CGEEVX1', IINFO, N, ISEED( 1 );
         }
         INFO = ( IINFO ).abs();
         return;
      }

      // Do Test (1)

      cget22('N', 'N', 'N', N, A, LDA, VR, LDVR, W, WORK, RWORK, RES );
      RESULT[1] = RES( 1 );

      // Do Test (2)

      cget22('C', 'N', 'C', N, A, LDA, VL, LDVL, W, WORK, RWORK, RES );
      RESULT[2] = RES( 1 );

      // Do Test (3)

      for (J = 1; J <= N; J++) { // 30
         TNRM = SCNRM2( N, VR( 1, J ), 1 );
         RESULT[3] = max( RESULT( 3 ), min( ULPINV, ( TNRM-ONE ).abs() / ULP ) );
         VMX = ZERO;
         VRMX = ZERO;
         for (JJ = 1; JJ <= N; JJ++) { // 20
            VTST = ( VR( JJ, J ) ).abs();
            if (VTST > VMX) VMX = VTST;
            IF( AIMAG( VR( JJ, J ) ) == ZERO && ABS( REAL( VR( JJ, J ) ) ) > VRMX ) VRMX = ABS( double( VR( JJ, J ) ) );
         } // 20
         if (VRMX / VMX < ONE-TWO*ULP) RESULT( 3 ) = ULPINV;
      } // 30

      // Do Test (4)

      for (J = 1; J <= N; J++) { // 50
         TNRM = SCNRM2( N, VL( 1, J ), 1 );
         RESULT[4] = max( RESULT( 4 ), min( ULPINV, ( TNRM-ONE ).abs() / ULP ) );
         VMX = ZERO;
         VRMX = ZERO;
         for (JJ = 1; JJ <= N; JJ++) { // 40
            VTST = ( VL( JJ, J ) ).abs();
            if (VTST > VMX) VMX = VTST;
            IF( AIMAG( VL( JJ, J ) ) == ZERO && ABS( REAL( VL( JJ, J ) ) ) > VRMX ) VRMX = ABS( double( VL( JJ, J ) ) );
         } // 40
         if (VRMX / VMX < ONE-TWO*ULP) RESULT( 4 ) = ULPINV;
      } // 50

      // Test for all options of computing condition numbers

      for (ISENS = 1; ISENS <= ISENSM; ISENS++) { // 200

         SENSE = SENS( ISENS );

         // Compute eigenvalues only, and test them

         clacpy('F', N, N, A, LDA, H, LDA );
         cgeevx(BALANC, 'N', 'N', SENSE, N, H, LDA, W1, CDUM, 1, CDUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX2', IINFO, N, JTYPE, BALANC, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX2', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 190;
         }

         // Do Test (5)

         for (J = 1; J <= N; J++) { // 60
            if( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV;
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

         clacpy('F', N, N, A, LDA, H, LDA );
         cgeevx(BALANC, 'N', 'V', SENSE, N, H, LDA, W1, CDUM, 1, LRE, LDLRE, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX3', IINFO, N, JTYPE, BALANC, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX3', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 190;
         }

         // Do Test (5) again

         for (J = 1; J <= N; J++) { // 90
            if( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV;
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

         clacpy('F', N, N, A, LDA, H, LDA );
         cgeevx(BALANC, 'V', 'N', SENSE, N, H, LDA, W1, LRE, LDLRE, CDUM, 1, ILO1, IHI1, SCALE1, ABNRM1, RCNDE1, RCNDV1, WORK, LWORK, RWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = ULPINV;
            if ( JTYPE != 22 ) {
               WRITE( NOUNIT, FMT = 9998 )'CGEEVX4', IINFO, N, JTYPE, BALANC, ISEED;
            } else {
               WRITE( NOUNIT, FMT = 9999 )'CGEEVX4', IINFO, N, ISEED( 1 );
            }
            INFO = ( IINFO ).abs();
            GO TO 190;
         }

         // Do Test (5) again

         for (J = 1; J <= N; J++) { // 140
            if( W( J ) != W1( J ) ) RESULT( 5 ) = ULPINV;
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
         clacpy('F', N, N, A, LDA, H, LDA );
         cgeevx('N', 'V', 'V', 'B', N, H, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IINFO );
         if ( IINFO != 0 ) {
            RESULT[1] = ULPINV;
            WRITE( NOUNIT, FMT = 9999 )'CGEEVX5', IINFO, N, ISEED( 1 );
            INFO = ( IINFO ).abs();
            GO TO 250;
         }

         // Sort eigenvalues and condition numbers lexicographically
         // to compare with inputs

         for (I = 1; I <= N - 1; I++) { // 220
            KMIN = I;
            if ( ISRT == 0 ) {
               VRIMIN = double( W( I ) );
            } else {
               VRIMIN = AIMAG( W( I ) );
            }
            for (J = I + 1; J <= N; J++) { // 210
               if ( ISRT == 0 ) {
                  VRICMP = double( W( J ) );
               } else {
                  VRICMP = AIMAG( W( J ) );
               }
               if ( VRICMP < VRIMIN ) {
                  KMIN = J;
                  VRIMIN = VRICMP;
               }
            } // 210
            CTMP = W( KMIN );
            W[KMIN] = W( I );
            W[I] = CTMP;
            VRIMIN = RCONDE( KMIN );
            RCONDE[KMIN] = RCONDE( I );
            RCONDE[I] = VRIMIN;
            VRIMIN = RCONDV( KMIN );
            RCONDV[KMIN] = RCONDV( I );
            RCONDV[I] = VRIMIN;
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

 9999 FORMAT( ' CGET23: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, INPUT EXAMPLE NUMBER = ${.i4}');
 9998 FORMAT( ' CGET23: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, BALANC = ${}, ISEED=(${.i5(4, ',')})' );

      }
