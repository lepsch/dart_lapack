import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlaed6(KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               ORGATI;
      int                INFO, KNITER;
      double             FINIT, RHO, TAU;
      double             D( 3 ), Z( 3 );
      // ..

      int                MAXIT;
      const              MAXIT = 40 ;
      double             ZERO, ONE, TWO, THREE, FOUR, EIGHT;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0, FOUR = 4.0, EIGHT = 8.0 ;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      double             DSCALE( 3 ), ZSCALE( 3 );
      bool               SCALE;
      int                I, ITER, NITER;
      double             A, B, BASE, C, DDF, DF, EPS, ERRETM, ETA, F, FC, SCLFAC, SCLINV, SMALL1, SMALL2, SMINV1, SMINV2, TEMP, TEMP1, TEMP2, TEMP3, TEMP4, LBD, UBD;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, MIN, SQRT

      INFO = 0;

      if ( ORGATI ) {
         LBD = D(2);
         UBD = D(3);
      } else {
         LBD = D(1);
         UBD = D(2);
      }
      if ( FINIT < ZERO ) {
         LBD = ZERO;
      } else {
         UBD = ZERO;
      }

      NITER = 1;
      TAU = ZERO;
      if ( KNITER == 2 ) {
         if ( ORGATI ) {
            TEMP = ( D( 3 )-D( 2 ) ) / TWO;
            C = RHO + Z( 1 ) / ( ( D( 1 )-D( 2 ) )-TEMP );
            A = C*( D( 2 )+D( 3 ) ) + Z( 2 ) + Z( 3 );
            B = C*D( 2 )*D( 3 ) + Z( 2 )*D( 3 ) + Z( 3 )*D( 2 );
         } else {
            TEMP = ( D( 1 )-D( 2 ) ) / TWO;
            C = RHO + Z( 3 ) / ( ( D( 3 )-D( 2 ) )-TEMP );
            A = C*( D( 1 )+D( 2 ) ) + Z( 1 ) + Z( 2 );
            B = C*D( 1 )*D( 2 ) + Z( 1 )*D( 2 ) + Z( 2 )*D( 1 );
         }
         TEMP = max( ( A ).abs(), ( B ).abs(), ( C ).abs() );
         A = A / TEMP;
         B = B / TEMP;
         C = C / TEMP;
         if ( C == ZERO ) {
            TAU = B / A;
         } else if ( A <= ZERO ) {
            TAU = ( A-sqrt( ( A*A-FOUR*B*C ).abs() ) ) / ( TWO*C );
         } else {
            TAU = TWO*B / ( A+sqrt( ( A*A-FOUR*B*C ).abs() ) );
         }
         if (TAU < LBD || TAU > UBD) TAU = ( LBD+UBD )/TWO;
         if ( D(1) == TAU || D(2) == TAU || D(3) == TAU ) {
            TAU = ZERO;
         } else {
            TEMP = FINIT + TAU*Z(1)/( D(1)*( D( 1 )-TAU ) ) + TAU*Z(2)/( D(2)*( D( 2 )-TAU ) ) + TAU*Z(3)/( D(3)*( D( 3 )-TAU ) );
            if ( TEMP <= ZERO ) {
               LBD = TAU;
            } else {
               UBD = TAU;
            }
            if( ( FINIT ).abs() <= ( TEMP ).abs() ) TAU = ZERO;
         }
      }

      // get machine parameters for possible scaling to avoid overflow

      // modified by Sven: parameters SMALL1, SMINV1, SMALL2,
      // SMINV2, EPS are not SAVEd anymore between one call to the
      // others but recomputed at each call

      EPS = dlamch( 'Epsilon' );
      BASE = dlamch( 'Base' );
      SMALL1 = BASE**( INT( LOG( dlamch( 'SafMin' ) ) / LOG( BASE ) / THREE ) );
      SMINV1 = ONE / SMALL1;
      SMALL2 = SMALL1*SMALL1;
      SMINV2 = SMINV1*SMINV1;

      // Determine if scaling of inputs necessary to avoid overflow
      // when computing 1/TEMP**3

      if ( ORGATI ) {
         TEMP = min( ABS( D( 2 )-TAU ), ( D( 3 )-TAU ).abs() );
      } else {
         TEMP = min( ABS( D( 1 )-TAU ), ( D( 2 )-TAU ).abs() );
      }
      SCALE = false;
      if ( TEMP <= SMALL1 ) {
         SCALE = true;
         if ( TEMP <= SMALL2 ) {

         // Scale up by power of radix nearest 1/SAFMIN**(2/3)

            SCLFAC = SMINV2;
            SCLINV = SMALL2;
         } else {

         // Scale up by power of radix nearest 1/SAFMIN**(1/3)

            SCLFAC = SMINV1;
            SCLINV = SMALL1;
         }

         // Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)

         for (I = 1; I <= 3; I++) { // 10
            DSCALE[I] = D( I )*SCLFAC;
            ZSCALE[I] = Z( I )*SCLFAC;
         } // 10
         TAU = TAU*SCLFAC;
         LBD = LBD*SCLFAC;
         UBD = UBD*SCLFAC;
      } else {

         // Copy D and Z to DSCALE and ZSCALE

         for (I = 1; I <= 3; I++) { // 20
            DSCALE[I] = D( I );
            ZSCALE[I] = Z( I );
         } // 20
      }

      FC = ZERO;
      DF = ZERO;
      DDF = ZERO;
      for (I = 1; I <= 3; I++) { // 30
         TEMP = ONE / ( DSCALE( I )-TAU );
         TEMP1 = ZSCALE( I )*TEMP;
         TEMP2 = TEMP1*TEMP;
         TEMP3 = TEMP2*TEMP;
         FC = FC + TEMP1 / DSCALE( I );
         DF = DF + TEMP2;
         DDF = DDF + TEMP3;
      } // 30
      F = FINIT + TAU*FC;

      if( ( F ).abs() <= ZERO ) GO TO 60;
      if ( F <= ZERO ) {
         LBD = TAU;
      } else {
         UBD = TAU;
      }

         // Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
                             // scheme

      // It is not hard to see that

            // 1) Iterations will go up monotonically
               // if FINIT < 0;

            // 2) Iterations will go down monotonically
               // if FINIT > 0.

      ITER = NITER + 1;

      for (NITER = ITER; NITER <= MAXIT; NITER++) { // 50

         if ( ORGATI ) {
            TEMP1 = DSCALE( 2 ) - TAU;
            TEMP2 = DSCALE( 3 ) - TAU;
         } else {
            TEMP1 = DSCALE( 1 ) - TAU;
            TEMP2 = DSCALE( 2 ) - TAU;
         }
         A = ( TEMP1+TEMP2 )*F - TEMP1*TEMP2*DF;
         B = TEMP1*TEMP2*F;
         C = F - ( TEMP1+TEMP2 )*DF + TEMP1*TEMP2*DDF;
         TEMP = max( ( A ).abs(), ( B ).abs(), ( C ).abs() );
         A = A / TEMP;
         B = B / TEMP;
         C = C / TEMP;
         if ( C == ZERO ) {
            ETA = B / A;
         } else if ( A <= ZERO ) {
            ETA = ( A-sqrt( ( A*A-FOUR*B*C ).abs() ) ) / ( TWO*C );
         } else {
            ETA = TWO*B / ( A+sqrt( ( A*A-FOUR*B*C ).abs() ) );
         }
         if ( F*ETA >= ZERO ) {
            ETA = -F / DF;
         }

         TAU = TAU + ETA;
         if (TAU < LBD || TAU > UBD) TAU = ( LBD + UBD )/TWO;

         FC = ZERO;
         ERRETM = ZERO;
         DF = ZERO;
         DDF = ZERO;
         for (I = 1; I <= 3; I++) { // 40
            if ( ( DSCALE( I )-TAU ) != ZERO ) {
               TEMP = ONE / ( DSCALE( I )-TAU );
               TEMP1 = ZSCALE( I )*TEMP;
               TEMP2 = TEMP1*TEMP;
               TEMP3 = TEMP2*TEMP;
               TEMP4 = TEMP1 / DSCALE( I );
               FC = FC + TEMP4;
               ERRETM = ERRETM + ( TEMP4 ).abs();
               DF = DF + TEMP2;
               DDF = DDF + TEMP3;
            } else {
               GO TO 60;
            }
         } // 40
         F = FINIT + TAU*FC;
         ERRETM = EIGHT*( ( FINIT ).abs()+( TAU ).abs()*ERRETM ) + ( TAU ).abs()*DF          IF( ( ( F ).abs() <= FOUR*EPS*ERRETM ) || ( (UBD-LBD) <= FOUR*EPS*(TAU).abs() )  ) GO TO 60;
         if ( F <= ZERO ) {
            LBD = TAU;
         } else {
            UBD = TAU;
         }
      } // 50
      INFO = 1;
      } // 60

      // Undo scaling

      if (SCALE) TAU = TAU*SCLINV;
      return;
      }
