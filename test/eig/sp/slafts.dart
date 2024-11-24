// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slafts(final int TYPE, final int M, final int N, final int IMAT, final int NTESTS, final int RESULT, final Array<int> ISEED_, final int THRESH, final int IOUNIT, final int IE,) {
  final ISEED = ISEED_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TYPE;
      int                IE, IMAT, IOUNIT, M, N, NTESTS;
      double               THRESH;
      int                ISEED( 4 );
      double               RESULT( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                K;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAHD2

      if ( M == N ) {

      // Output for square matrices:

         for (K = 1; K <= NTESTS; K++) { // 10
            if ( RESULT( K ) >= THRESH ) {

            // If this is the first test to fail, call SLAHD2
            // to print a header to the data file.

               if (IE == 0) slahd2( IOUNIT, TYPE );
               IE = IE + 1;
               if ( RESULT( K ) < 10000.0 ) {
                  WRITE( IOUNIT, FMT = 9999 )N, IMAT, ISEED, K, RESULT( K );
 9999             FORMAT( ' Matrix order=${.i5}, type=${.i2}, seed=${i4(4, ',')}', ' result ${.i3} is${.f8_2}');
               } else {
                  WRITE( IOUNIT, FMT = 9998 )N, IMAT, ISEED, K, RESULT( K );
 9998             FORMAT( ' Matrix order=${.i5}, type=${.i2}, seed=${i4(4, ',')}', ' result ${.i3} is', 1P, E10.3 );
               }
            }
         } // 10
      } else {

      // Output for rectangular matrices

         for (K = 1; K <= NTESTS; K++) { // 20
            if ( RESULT( K ) >= THRESH ) {

               // If this is the first test to fail, call SLAHD2
               // to print a header to the data file.

               if (IE == 0) slahd2( IOUNIT, TYPE );
               IE = IE + 1;
               if ( RESULT( K ) < 10000.0 ) {
                  WRITE( IOUNIT, FMT = 9997 )M, N, IMAT, ISEED, K, RESULT( K );
 9997             FORMAT(' ${.i5} x${.i5} matrix, type=${.i2}, s', 'eed=${i4(3, ',')}', I4, ': result ${.i3} is${.f8_2}');
               } else {
                  WRITE( IOUNIT, FMT = 9996 )M, N, IMAT, ISEED, K, RESULT( K );
 9996             FORMAT(' ${.i5} x${.i5} matrix, type=${.i2}, s', 'eed=${i4(3, ',')}', I4, ': result ${.i3} is', 1P, E10.3 );
               }
            }
         } // 20

      }
      }
