import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgtsv(N, NRHS, DL, D, DU, B, LDB, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDB, N, NRHS;
      double             B( LDB, * ), D( * ), DL( * ), DU( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      int                I, J;
      double             FACT, TEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( NRHS < 0 ) {
         INFO = -2;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('DGTSV ', -INFO );
         return;
      }

      if (N == 0) return;

      if ( NRHS == 1 ) {
         for (I = 1; I <= N - 2; I++) { // 10
            if ( ( D( I ) ).abs() >= ( DL( I ) ).abs() ) {

               // No row interchange required

               if ( D( I ) != ZERO ) {
                  FACT = DL( I ) / D( I );
                  D[I+1] = D( I+1 ) - FACT*DU( I );
                  B[I+1][1] = B( I+1, 1 ) - FACT*B( I, 1 );
               } else {
                  INFO = I;
                  return;
               }
               DL[I] = ZERO;
            } else {

               // Interchange rows I and I+1

               FACT = D( I ) / DL( I );
               D[I] = DL( I );
               TEMP = D( I+1 );
               D[I+1] = DU( I ) - FACT*TEMP;
               DL[I] = DU( I+1 );
               DU[I+1] = -FACT*DL( I );
               DU[I] = TEMP;
               TEMP = B( I, 1 );
               B[I][1] = B( I+1, 1 );
               B[I+1][1] = TEMP - FACT*B( I+1, 1 );
            }
         } // 10
         if ( N > 1 ) {
            I = N - 1;
            if ( ( D( I ) ).abs() >= ( DL( I ) ).abs() ) {
               if ( D( I ) != ZERO ) {
                  FACT = DL( I ) / D( I );
                  D[I+1] = D( I+1 ) - FACT*DU( I );
                  B[I+1][1] = B( I+1, 1 ) - FACT*B( I, 1 );
               } else {
                  INFO = I;
                  return;
               }
            } else {
               FACT = D( I ) / DL( I );
               D[I] = DL( I );
               TEMP = D( I+1 );
               D[I+1] = DU( I ) - FACT*TEMP;
               DU[I] = TEMP;
               TEMP = B( I, 1 );
               B[I][1] = B( I+1, 1 );
               B[I+1][1] = TEMP - FACT*B( I+1, 1 );
            }
         }
         if ( D( N ) == ZERO ) {
            INFO = N;
            return;
         }
      } else {
         for (I = 1; I <= N - 2; I++) { // 40
            if ( ( D( I ) ).abs() >= ( DL( I ) ).abs() ) {

               // No row interchange required

               if ( D( I ) != ZERO ) {
                  FACT = DL( I ) / D( I );
                  D[I+1] = D( I+1 ) - FACT*DU( I );
                  for (J = 1; J <= NRHS; J++) { // 20
                     B[I+1][J] = B( I+1, J ) - FACT*B( I, J );
                  } // 20
               } else {
                  INFO = I;
                  return;
               }
               DL[I] = ZERO;
            } else {

               // Interchange rows I and I+1

               FACT = D( I ) / DL( I );
               D[I] = DL( I );
               TEMP = D( I+1 );
               D[I+1] = DU( I ) - FACT*TEMP;
               DL[I] = DU( I+1 );
               DU[I+1] = -FACT*DL( I );
               DU[I] = TEMP;
               for (J = 1; J <= NRHS; J++) { // 30
                  TEMP = B( I, J );
                  B[I][J] = B( I+1, J );
                  B[I+1][J] = TEMP - FACT*B( I+1, J );
               } // 30
            }
         } // 40
         if ( N > 1 ) {
            I = N - 1;
            if ( ( D( I ) ).abs() >= ( DL( I ) ).abs() ) {
               if ( D( I ) != ZERO ) {
                  FACT = DL( I ) / D( I );
                  D[I+1] = D( I+1 ) - FACT*DU( I );
                  for (J = 1; J <= NRHS; J++) { // 50
                     B[I+1][J] = B( I+1, J ) - FACT*B( I, J );
                  } // 50
               } else {
                  INFO = I;
                  return;
               }
            } else {
               FACT = D( I ) / DL( I );
               D[I] = DL( I );
               TEMP = D( I+1 );
               D[I+1] = DU( I ) - FACT*TEMP;
               DU[I] = TEMP;
               for (J = 1; J <= NRHS; J++) { // 60
                  TEMP = B( I, J );
                  B[I][J] = B( I+1, J );
                  B[I+1][J] = TEMP - FACT*B( I+1, J );
               } // 60
            }
         }
         if ( D( N ) == ZERO ) {
            INFO = N;
            return;
         }
      }

      // Back solve with the matrix U from the factorization.

      if ( NRHS <= 2 ) {
         J = 1;
         } // 70
         B[N][J] = B( N, J ) / D( N );
         if (N > 1) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 );
         for (I = N - 2; I >= 1; I--) { // 80
            B[I][J] = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* B( I+2, J ) ) / D( I );
         } // 80
         if ( J < NRHS ) {
            J = J + 1;
            GO TO 70;
         }
      } else {
         for (J = 1; J <= NRHS; J++) { // 100
            B[N][J] = B( N, J ) / D( N );
            if (N > 1) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 );
            for (I = N - 2; I >= 1; I--) { // 90
               B[I][J] = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )* B( I+2, J ) ) / D( I );
            } // 90
         } // 100
      }

      }
