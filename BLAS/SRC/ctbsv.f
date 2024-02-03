      SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INCX,K,LDA,N;
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
      int     I,INFO,IX,J,JX,KPLUS1,KX,L;
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
      // INTRINSIC CONJG,MAX,MIN
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
      } else if (K.LT.0) {
          INFO = 5
      } else if (LDA.LT. (K+1)) {
          INFO = 7
      } else if (INCX.EQ.0) {
          INFO = 9
      }
      if (INFO.NE.0) {
          CALL XERBLA('CTBSV ',INFO)
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
      // accessed by sequentially with one pass through A.

      if (LSAME(TRANS,'N')) {

         // Form  x := inv( A )*x.

          if (LSAME(UPLO,'U')) {
              KPLUS1 = K + 1
              if (INCX.EQ.1) {
                  DO 20 J = N,1,-1
                      if (X(J).NE.ZERO) {
                          L = KPLUS1 - J
                          IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,MAX(1,J-K),-1
                              X(I) = X(I) - TEMP*A(L+I,J)
   10                     CONTINUE
                      }
   20             CONTINUE
              } else {
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                      KX = KX - INCX
                      if (X(JX).NE.ZERO) {
                          IX = KX
                          L = KPLUS1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                          TEMP = X(JX)
                          DO 30 I = J - 1,MAX(1,J-K),-1
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX - INCX
   30                     CONTINUE
                      }
                      JX = JX - INCX
   40             CONTINUE
              }
          } else {
              if (INCX.EQ.1) {
                  DO 60 J = 1,N
                      if (X(J).NE.ZERO) {
                          L = 1 - J
                          IF (NOUNIT) X(J) = X(J)/A(1,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,MIN(N,J+K)
                              X(I) = X(I) - TEMP*A(L+I,J)
   50                     CONTINUE
                      }
   60             CONTINUE
              } else {
                  JX = KX
                  DO 80 J = 1,N
                      KX = KX + INCX
                      if (X(JX).NE.ZERO) {
                          IX = KX
                          L = 1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                          TEMP = X(JX)
                          DO 70 I = J + 1,MIN(N,J+K)
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX + INCX
   70                     CONTINUE
                      }
                      JX = JX + INCX
   80             CONTINUE
              }
          }
      } else {

         // Form  x := inv( A**T )*x  or  x := inv( A**H )*x.

          if (LSAME(UPLO,'U')) {
              KPLUS1 = K + 1
              if (INCX.EQ.1) {
                  DO 110 J = 1,N
                      TEMP = X(J)
                      L = KPLUS1 - J
                      if (NOCONJ) {
                          DO 90 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(I)
   90                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      } else {
                          DO 100 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - CONJG(A(L+I,J))*X(I)
  100                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(KPLUS1,J))
                      }
                      X(J) = TEMP
  110             CONTINUE
              } else {
                  JX = KX
                  DO 140 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      L = KPLUS1 - J
                      if (NOCONJ) {
                          DO 120 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX + INCX
  120                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      } else {
                          DO 130 I = MAX(1,J-K),J - 1
                              TEMP = TEMP - CONJG(A(L+I,J))*X(IX)
                              IX = IX + INCX
  130                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(KPLUS1,J))
                      }
                      X(JX) = TEMP
                      JX = JX + INCX
                      IF (J.GT.K) KX = KX + INCX
  140             CONTINUE
              }
          } else {
              if (INCX.EQ.1) {
                  DO 170 J = N,1,-1
                      TEMP = X(J)
                      L = 1 - J
                      if (NOCONJ) {
                          DO 150 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(I)
  150                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(1,J)
                      } else {
                          DO 160 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - CONJG(A(L+I,J))*X(I)
  160                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(1,J))
                      }
                      X(J) = TEMP
  170             CONTINUE
              } else {
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 200 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      L = 1 - J
                      if (NOCONJ) {
                          DO 180 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - A(L+I,J)*X(IX)
                              IX = IX - INCX
  180                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(1,J)
                      } else {
                          DO 190 I = MIN(N,J+K),J + 1,-1
                              TEMP = TEMP - CONJG(A(L+I,J))*X(IX)
                              IX = IX - INCX
  190                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/CONJG(A(1,J))
                      }
                      X(JX) = TEMP
                      JX = JX - INCX
                      IF ((N-J).GE.K) KX = KX - INCX
  200             CONTINUE
              }
          }
      }

      RETURN

      // End of CTBSV

      }
