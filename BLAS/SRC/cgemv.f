      SUBROUTINE CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

*  -- Reference BLAS level2 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX ALPHA,BETA
      int     INCX,INCY,LDA,M,N;
      String    TRANS;
      // ..
      // .. Array Arguments ..
      COMPLEX A(LDA,*),X(*),Y(*)
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX ONE
      const     ONE= (1.0E+0,0.0E+0);
      COMPLEX ZERO
      const     ZERO= (0.0E+0,0.0E+0);
      // ..
      // .. Local Scalars ..
      COMPLEX TEMP
      int     I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY;
      bool    NOCONJ;
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
      if (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.LSAME(TRANS,'C')) {
          INFO = 1
      } else if (M.LT.0) {
          INFO = 2
      } else if (N.LT.0) {
          INFO = 3
      } else if (LDA.LT.MAX(1,M)) {
          INFO = 6
      } else if (INCX.EQ.0) {
          INFO = 8
      } else if (INCY.EQ.0) {
          INFO = 11
      }
      if (INFO.NE.0) {
          xerbla('CGEMV ',INFO);
          RETURN
      }

      // Quick return if possible.

      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN

      NOCONJ = LSAME(TRANS,'T')

      // Set  LENX  and  LENY, the lengths of the vectors x and y, and set
      // up the start points in  X  and  Y.

      if (LSAME(TRANS,'N')) {
          LENX = N
          LENY = M
      } else {
          LENX = M
          LENY = N
      }
      if (INCX.GT.0) {
          KX = 1
      } else {
          KX = 1 - (LENX-1)*INCX
      }
      if (INCY.GT.0) {
          KY = 1
      } else {
          KY = 1 - (LENY-1)*INCY
      }

      // Start the operations. In this version the elements of A are
      // accessed sequentially with one pass through A.

      // First form  y := beta*y.

      if (BETA.NE.ONE) {
          if (INCY.EQ.1) {
              if (BETA.EQ.ZERO) {
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              } else {
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              }
          } else {
              IY = KY
              if (BETA.EQ.ZERO) {
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              } else {
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              }
          }
      }
      IF (ALPHA.EQ.ZERO) RETURN
      if (LSAME(TRANS,'N')) {

         // Form  y := alpha*A*x + y.

          JX = KX
          if (INCY.EQ.1) {
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          } else {
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          }
      } else {

         // Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.

          JY = KY
          if (INCX.EQ.1) {
              DO 110 J = 1,N
                  TEMP = ZERO
                  if (NOCONJ) {
                      DO 90 I = 1,M
                          TEMP = TEMP + A(I,J)*X(I)
   90                 CONTINUE
                  } else {
                      DO 100 I = 1,M
                          TEMP = TEMP + CONJG(A(I,J))*X(I)
  100                 CONTINUE
                  }
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  110         CONTINUE
          } else {
              DO 140 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  if (NOCONJ) {
                      DO 120 I = 1,M
                          TEMP = TEMP + A(I,J)*X(IX)
                          IX = IX + INCX
  120                 CONTINUE
                  } else {
                      DO 130 I = 1,M
                          TEMP = TEMP + CONJG(A(I,J))*X(IX)
                          IX = IX + INCX
  130                 CONTINUE
                  }
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  140         CONTINUE
          }
      }

      RETURN

      // End of CGEMV

      }
