// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void stpsv(final int UPLO, final int TRANS, final int DIAG, final int N, final int AP, final int X, final int INCX,) {

// -- Reference BLAS level2 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     INCX,N;
      String    DIAG,TRANS,UPLO;
      double AP(*),X(*);
      // ..

      double ZERO;
      const     ZERO=0.0;
      double TEMP;
      int     I,INFO,IX,J,JX,K,KK,KX;
      bool    NOUNIT;
      // ..
      // .. External Functions ..
      //- bool    lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..

      // Test the input parameters.

      INFO = 0;
      if ( !lsame(UPLO,'U') && !lsame(UPLO,'L')) {
          INFO = 1;
      } else if ( !lsame(TRANS,'N') && !lsame(TRANS,'T') && !lsame(TRANS,'C')) {
          INFO = 2;
      } else if ( !lsame(DIAG,'U') && !lsame(DIAG,'N')) {
          INFO = 3;
      } else if (N < 0) {
          INFO = 4;
      } else if (INCX == 0) {
          INFO = 7;
      }
      if (INFO != 0) {
          xerbla('STPSV ',INFO);
          return;
      }

      // Quick return if possible.

      if (N == 0) return;

      NOUNIT = lsame(DIAG,'N');

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX <= 0) {
          KX = 1 - (N-1)*INCX;
      } else if (INCX != 1) {
          KX = 1;
      }

      // Start the operations. In this version the elements of AP are
      // accessed sequentially with one pass through AP.

      if (lsame(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (lsame(UPLO,'U')) {
              KK = (N* (N+1))/2;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 20
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/AP(KK);
                          TEMP = X(J);
                          K = KK - 1;
                          for (I = J - 1; I >= 1; I--) { // 10
                              X[I] = X(I) - TEMP*AP(K);
                              K = K - 1;
                          } // 10
                      }
                      KK = KK - J;
                  } // 20
              } else {
                  JX = KX + (N-1)*INCX;
                  for (J = N; J >= 1; J--) { // 40
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/AP(KK);
                          TEMP = X(JX);
                          IX = JX;
                          for (K = KK - 1; K >= KK - J + 1; K--) { // 30
                              IX = IX - INCX;
                              X[IX] = X(IX) - TEMP*AP(K);
                          } // 30
                      }
                      JX = JX - INCX;
                      KK = KK - J;
                  } // 40
              }
          } else {
              KK = 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 60
                      if (X(J) != ZERO) {
                          if (NOUNIT) X(J) = X(J)/AP(KK);
                          TEMP = X(J);
                          K = KK + 1;
                          for (I = J + 1; I <= N; I++) { // 50
                              X[I] = X(I) - TEMP*AP(K);
                              K = K + 1;
                          } // 50
                      }
                      KK = KK + (N-J+1);
                  } // 60
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 80
                      if (X(JX) != ZERO) {
                          if (NOUNIT) X(JX) = X(JX)/AP(KK);
                          TEMP = X(JX);
                          IX = JX;
                          for (K = KK + 1; K <= KK + N - J; K++) { // 70
                              IX = IX + INCX;
                              X[IX] = X(IX) - TEMP*AP(K);
                          } // 70
                      }
                      JX = JX + INCX;
                      KK = KK + (N-J+1);
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T )*x.

          if (lsame(UPLO,'U')) {
              KK = 1;
              if (INCX == 1) {
                  for (J = 1; J <= N; J++) { // 100
                      TEMP = X(J);
                      K = KK;
                      for (I = 1; I <= J - 1; I++) { // 90
                          TEMP = TEMP - AP(K)*X(I);
                          K = K + 1;
                      } // 90
                      if (NOUNIT) TEMP = TEMP/AP(KK+J-1);
                      X[J] = TEMP;
                      KK = KK + J;
                  } // 100
              } else {
                  JX = KX;
                  for (J = 1; J <= N; J++) { // 120
                      TEMP = X(JX);
                      IX = KX;
                      for (K = KK; K <= KK + J - 2; K++) { // 110
                          TEMP = TEMP - AP(K)*X(IX);
                          IX = IX + INCX;
                      } // 110
                      if (NOUNIT) TEMP = TEMP/AP(KK+J-1);
                      X[JX] = TEMP;
                      JX = JX + INCX;
                      KK = KK + J;
                  } // 120
              }
          } else {
              KK = (N* (N+1))/2;
              if (INCX == 1) {
                  for (J = N; J >= 1; J--) { // 140
                      TEMP = X(J);
                      K = KK;
                      for (I = N; I >= J + 1; I--) { // 130
                          TEMP = TEMP - AP(K)*X(I);
                          K = K - 1;
                      } // 130
                      if (NOUNIT) TEMP = TEMP/AP(KK-N+J);
                      X[J] = TEMP;
                      KK = KK - (N-J+1);
                  } // 140
              } else {
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for (J = N; J >= 1; J--) { // 160
                      TEMP = X(JX);
                      IX = KX;
                      for (K = KK; K >= KK - (N- (J+1)); K--) { // 150
                          TEMP = TEMP - AP(K)*X(IX);
                          IX = IX - INCX;
                      } // 150
                      if (NOUNIT) TEMP = TEMP/AP(KK-N+J);
                      X[JX] = TEMP;
                      JX = JX - INCX;
                      KK = KK - (N-J+1);
                  } // 160
              }
          }
      }

      }
