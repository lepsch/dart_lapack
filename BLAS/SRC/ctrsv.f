      SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,LDA,N;
      String    DIAG,TRANS,UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),X(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP
      int     I,INFO,IX,J,JX,KX;
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
      // INTRINSIC CONJG,MAX
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
      } else if (LDA.LT.MAX(1,N)) {
          INFO = 6
      } else if (INCX.EQ.0) {
          INFO = 8
      }
      if (INFO.NE.0) {
          xerbla('CTRSV ',INFO);
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

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      if (LSAME(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (LSAME(UPLO,'U')) {
              if (INCX.EQ.1) {
                  DO 20 J = N,1,-1
                      if (X(J).NE.ZERO) {
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,1,-1
                              X(I) = X(I) - TEMP*A(I,J)
                          } // 10
                      }
                  } // 20
              } else {
                  JX = KX + (N-1)*INCX
                  DO 40 J = N,1,-1
                      if (X(JX).NE.ZERO) {
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 30 I = J - 1,1,-1
                              IX = IX - INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
                          } // 30
                      }
                      JX = JX - INCX
                  } // 40
              }
          } else {
              if (INCX.EQ.1) {
                  for (J = 1; J <= N; J++) { // 60
                      if (X(J).NE.ZERO) {
                          IF (NOUNIT) X(J) = X(J)/A(J,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,N
                              X(I) = X(I) - TEMP*A(I,J)
                          } // 50
                      }
                  } // 60
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 80
                      if (X(JX).NE.ZERO) {
                          IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                          TEMP = X(JX)
                          IX = JX
                          DO 70 I = J + 1,N
                              IX = IX + INCX
                              X(IX) = X(IX) - TEMP*A(I,J)
                          } // 70
                      }
                      JX = JX + INCX
                  } // 80
              }
          }
      } else {

         // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.

          if (LSAME(UPLO,'U')) {
              if (INCX.EQ.1) {
                  for (J = 1; J <= N; J++) { // 110
                      TEMP = X(J)
                      if (NOCONJ) {
                          DO 90 I = 1,J - 1
                              TEMP = TEMP - A(I,J)*X(I)
                          } // 90
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      } else {
                          DO 100 I = 1,J - 1
                              TEMP = TEMP - CONJG(A(I,J))*X(I)
                          } // 100
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
                      }
                      X(J) = TEMP
                  } // 110
              } else {
                  JX = KX
                  for (J = 1; J <= N; J++) { // 140
                      IX = KX
                      TEMP = X(JX)
                      if (NOCONJ) {
                          DO 120 I = 1,J - 1
                              TEMP = TEMP - A(I,J)*X(IX)
                              IX = IX + INCX
                          } // 120
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      } else {
                          DO 130 I = 1,J - 1
                              TEMP = TEMP - CONJG(A(I,J))*X(IX)
                              IX = IX + INCX
                          } // 130
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
                      }
                      X(JX) = TEMP
                      JX = JX + INCX
                  } // 140
              }
          } else {
              if (INCX.EQ.1) {
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      if (NOCONJ) {
                          DO 150 I = N,J + 1,-1
                              TEMP = TEMP - A(I,J)*X(I)
                          } // 150
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      } else {
                          DO 160 I = N,J + 1,-1
                              TEMP = TEMP - CONJG(A(I,J))*X(I)
                          } // 160
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
                      }
                      X(J) = TEMP
                  } // 170
              } else {
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      IX = KX
                      TEMP = X(JX)
                      if (NOCONJ) {
                          DO 180 I = N,J + 1,-1
                              TEMP = TEMP - A(I,J)*X(IX)
                              IX = IX - INCX
                          } // 180
                          IF (NOUNIT) TEMP = TEMP/A(J,J)
                      } else {
                          DO 190 I = N,J + 1,-1
                              TEMP = TEMP - CONJG(A(I,J))*X(IX)
                              IX = IX - INCX
                          } // 190
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(J,J))
                      }
                      X(JX) = TEMP
                      JX = JX - INCX
                  } // 200
              }
          }
      }

      RETURN

      // End of CTRSV

      }
