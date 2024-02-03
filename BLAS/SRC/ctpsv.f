      SUBROUTINE CTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX AP(*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP
      int     I,INFO,IX,J,JX,K,KK,KX;
      bool    NOCONJ,NOUNIT;
      // ..
      // .. External Functions ..
      bool    LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..

      // Test the input parameters.

      INFO = 0
      if (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) {
          INFO = 1
      } else if (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.LSAME(TRANS,'C')) {
          INFO = 2
      } else if (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) {
          INFO = 3
      } else if (N.LT.0) {
          INFO = 4
      } else if (INCX.EQ.0) {
          INFO = 7
      }
      if (INFO.NE.0) {
          xerbla('CTPSV ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF (N.EQ.0) RETURN

      NOCONJ = LSAME(TRANS,'T')
      NOUNIT = LSAME(DIAG,'N')

      // Set up the start point in X if the increment is not unity. This
      // will be  ( N - 1 )*INCX  too small for descending loops.

      if (INCX.LE.0) {
          KX = 1 - (N-1)*INCX
      } else if (INCX.NE.1) {
          KX = 1
      }

      // Start the operations. In this version the elements of AP are
      // accessed sequentially with one pass through AP.

      if (LSAME(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (LSAME(UPLO,'U')) {
              KK = (N* (N+1))/2
              if (INCX.EQ.1) {
                  DO 20 J = N,1,-1
                      if (X(J).NE.ZERO) {
                          IF (NOUNIT) X(J) = X(J)/AP(KK)
                          TEMP = X(J)
                          K = KK - 1
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*AP(K)
                              K = K - 1
                          } // 10
                      }
                      KK = KK - J
                  } // 20
              } else {
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      if (X(JX).NE.ZERO) {
                          IF (NOUNIT) X(JX) = X(JX)/AP(KK)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 K = KK - 1,KK - J + 1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*AP(K)
                          } // 30
                      }
                      JX = JX - INCX
                      KK = KK - J
                  } // 40
              }
          } else {
              KK = 1
              if (INCX.EQ.1) {
                  for (J = 1; J <= N; J++) { // 60
                      if (X(J).NE.ZERO) {
                          IF (NOUNIT) X(J) = X(J)/AP(KK)
                          TEMP = X(J)
                          K = KK + 1
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*AP(K)
                              K = K + 1
                          } // 50
                      }
                      KK = KK + (N-J+1)
                  } // 60
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 80
                      if (X(JX).NE.ZERO) {
                          IF (NOUNIT) X(JX) = X(JX)/AP(KK)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 K = KK + 1,KK + N - J
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*AP(K)
                          } // 70
                      }
                      JX = JX + INCX
                      KK = KK + (N-J+1)
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.

          if (LSAME(UPLO,'U')) {
              KK = 1
              if (INCX.EQ.1) {
                  for (J = 1; J <= N; J++) { // 110
                      TEMP = X(J)
                      K = KK
                      if (NOCONJ) {
                          DO 90 I = 1,J - 1
                              TEMP = TEMP - AP(K)*X(I)
                              K = K + 1
                          } // 90
                          IF (NOUNIT) TEMP = TEMP/AP(KK+J-1)
                      } else {
                          DO 100 I = 1,J - 1
                              TEMP = TEMP - CONJG(AP(K))*X(I)
                              K = K + 1
                          } // 100
                          IF (NOUNIT) TEMP = TEMP/CONJG(AP(KK+J-1))
                      }
                      X(J) = TEMP
                      KK = KK + J
                  } // 110
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 140
                      TEMP = X(JX)
                      IX = KX
                      if (NOCONJ) {
                          DO 120 K = KK,KK + J - 2
                              TEMP = TEMP - AP(K)*X(IX)
                              IX = IX + INCX
                          } // 120
                          IF (NOUNIT) TEMP = TEMP/AP(KK+J-1)
                      } else {
                          DO 130 K = KK,KK + J - 2
                              TEMP = TEMP - CONJG(AP(K))*X(IX)
                              IX = IX + INCX
                          } // 130
                          IF (NOUNIT) TEMP = TEMP/CONJG(AP(KK+J-1))
                      }
                      X(JX) = TEMP
                      JX = JX + INCX
                      KK = KK + J
                  } // 140
              }
          } else {
              KK = (N* (N+1))/2
              if (INCX.EQ.1) {
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      K = KK
                      if (NOCONJ) {
                          DO 150 I = N,J + 1,-1
                              TEMP = TEMP - AP(K)*X(I)
                              K = K - 1
                          } // 150
                          IF (NOUNIT) TEMP = TEMP/AP(KK-N+J)
                      } else {
                          DO 160 I = N,J + 1,-1
                              TEMP = TEMP - CONJG(AP(K))*X(I)
                              K = K - 1
                          } // 160
                          IF (NOUNIT) TEMP = TEMP/CONJG(AP(KK-N+J))
                      }
                      X(J) = TEMP
                      KK = KK - (N-J+1)
                  } // 170
              } else {
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      if (NOCONJ) {
                          DO 180 K = KK,KK - (N- (J+1)),-1
                              TEMP = TEMP - AP(K)*X(IX)
                              IX = IX - INCX
                          } // 180
                          IF (NOUNIT) TEMP = TEMP/AP(KK-N+J)
                      } else {
                          DO 190 K = KK,KK - (N- (J+1)),-1
                              TEMP = TEMP - CONJG(AP(K))*X(IX)
                              IX = IX - INCX
                          } // 190
                          IF (NOUNIT) TEMP = TEMP/CONJG(AP(KK-N+J))
                      }
                      X(JX) = TEMP
                      JX = JX - INCX
                      KK = KK - (N-J+1)
                  } // 200
              }
          }
      }

      RETURN

      // End of CTPSV

      }
