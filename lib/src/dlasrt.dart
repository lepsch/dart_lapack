import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlasrt(ID, N, D, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             ID;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             D( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                SELECT;
      const              SELECT = 20 ;
      // ..
      // .. Local Scalars ..
      int                DIR, ENDD, I, J, START, STKPNT;
      double             D1, D2, D3, DMNMX, TMP;
      // ..
      // .. Local Arrays ..
      int                STACK( 2, 32 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      DIR = -1;
      if ( lsame( ID, 'D' ) ) {
         DIR = 0;
      } else if ( lsame( ID, 'I' ) ) {
         DIR = 1;
      }
      if ( DIR == -1 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('DLASRT', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 1) return;

      STKPNT = 1;
      STACK[1, 1] = 1;
      STACK[2, 1] = N;
      } // 10
      START = STACK( 1, STKPNT );
      ENDD = STACK( 2, STKPNT );
      STKPNT = STKPNT - 1;
      if ( ENDD-START <= SELECT && ENDD-START > 0 ) {

         // Do Insertion sort on D( START:ENDD )

         if ( DIR == 0 ) {

            // Sort into decreasing order

            for (I = START + 1; I <= ENDD; I++) { // 30
               for (J = I; J >= START + 1; J--) { // 20
                  if ( D( J ) > D( J-1 ) ) {
                     DMNMX = D( J );
                     D[J] = D( J-1 );
                     D[J-1] = DMNMX;
                  } else {
                     GO TO 30;
                  }
               } // 20
            } // 30

         } else {

            // Sort into increasing order

            for (I = START + 1; I <= ENDD; I++) { // 50
               for (J = I; J >= START + 1; J--) { // 40
                  if ( D( J ) < D( J-1 ) ) {
                     DMNMX = D( J );
                     D[J] = D( J-1 );
                     D[J-1] = DMNMX;
                  } else {
                     GO TO 50;
                  }
               } // 40
            } // 50

         }

      } else if ( ENDD-START > SELECT ) {

         // Partition D( START:ENDD ) and stack parts, largest one first

         // Choose partition entry as median of 3

         D1 = D( START );
         D2 = D( ENDD );
         I = ( START+ENDD ) / 2;
         D3 = D( I );
         if ( D1 < D2 ) {
            if ( D3 < D1 ) {
               DMNMX = D1;
            } else if ( D3 < D2 ) {
               DMNMX = D3;
            } else {
               DMNMX = D2;
            }
         } else {
            if ( D3 < D2 ) {
               DMNMX = D2;
            } else if ( D3 < D1 ) {
               DMNMX = D3;
            } else {
               DMNMX = D1;
            }
         }

         if ( DIR == 0 ) {

            // Sort into decreasing order

            I = START - 1;
            J = ENDD + 1;
            } // 60
            } // 70
            J = J - 1;
            if( D( J ) < DMNMX ) GO TO 70;
            } // 80
            I = I + 1;
            if( D( I ) > DMNMX ) GO TO 80;
            if ( I < J ) {
               TMP = D( I );
               D[I] = D( J );
               D[J] = TMP;
               GO TO 60;
            }
            if ( J-START > ENDD-J-1 ) {
               STKPNT = STKPNT + 1;
               STACK[1, STKPNT] = START;
               STACK[2, STKPNT] = J;
               STKPNT = STKPNT + 1;
               STACK[1, STKPNT] = J + 1;
               STACK[2, STKPNT] = ENDD;
            } else {
               STKPNT = STKPNT + 1;
               STACK[1, STKPNT] = J + 1;
               STACK[2, STKPNT] = ENDD;
               STKPNT = STKPNT + 1;
               STACK[1, STKPNT] = START;
               STACK[2, STKPNT] = J;
            }
         } else {

            // Sort into increasing order

            I = START - 1;
            J = ENDD + 1;
            } // 90
            } // 100
            J = J - 1;
            if( D( J ) > DMNMX ) GO TO 100;
            } // 110
            I = I + 1;
            if( D( I ) < DMNMX ) GO TO 110;
            if ( I < J ) {
               TMP = D( I );
               D[I] = D( J );
               D[J] = TMP;
               GO TO 90;
            }
            if ( J-START > ENDD-J-1 ) {
               STKPNT = STKPNT + 1;
               STACK[1, STKPNT] = START;
               STACK[2, STKPNT] = J;
               STKPNT = STKPNT + 1;
               STACK[1, STKPNT] = J + 1;
               STACK[2, STKPNT] = ENDD;
            } else {
               STKPNT = STKPNT + 1;
               STACK[1, STKPNT] = J + 1;
               STACK[2, STKPNT] = ENDD;
               STKPNT = STKPNT + 1;
               STACK[1, STKPNT] = START;
               STACK[2, STKPNT] = J;
            }
         }
      }
      if (STKPNT > 0) GO TO 10;
      return;
      }
