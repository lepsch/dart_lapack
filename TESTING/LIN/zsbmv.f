      SUBROUTINE ZSBMV( UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCX, INCY, K, LDA, N;
      COMPLEX*16         ALPHA, BETA
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                I, INFO, IX, IY, J, JX, JY, KPLUS1, KX, KY, L;
      COMPLEX*16         TEMP1, TEMP2
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = 1
      ELSE IF( N.LT.0 ) THEN
         INFO = 2
      ELSE IF( K.LT.0 ) THEN
         INFO = 3
      ELSE IF( LDA.LT.( K+1 ) ) THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 ) THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 ) THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSBMV ', INFO )
         RETURN
      END IF

      // Quick return if possible.

      IF( ( N.EQ.0 ) .OR. ( ( ALPHA.EQ.ZERO ) .AND. ( BETA.EQ.ONE ) ) ) RETURN

      // Set up the start points in  X  and  Y.

      IF( INCX.GT.0 ) THEN
         KX = 1
      ELSE
         KX = 1 - ( N-1 )*INCX
      END IF
      IF( INCY.GT.0 ) THEN
         KY = 1
      ELSE
         KY = 1 - ( N-1 )*INCY
      END IF

      // Start the operations. In this version the elements of the array A
      // are accessed sequentially with one pass through A.

      // First form  y := beta*y.

      IF( BETA.NE.ONE ) THEN
         IF( INCY.EQ.1 ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 10 I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20 I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO ) THEN
               DO 30 I = 1, N
                  Y( IY ) = ZERO
                  IY = IY + INCY
   30          CONTINUE
            ELSE
               DO 40 I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY = IY + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO ) RETURN
      IF( LSAME( UPLO, 'U' ) ) THEN

         // Form  y  when upper triangle of A is stored.

         KPLUS1 = K + 1
         IF( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) THEN
            DO 60 J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               L = KPLUS1 - J
               DO 50 I = MAX( 1, J-K ), J - 1
                  Y( I ) = Y( I ) + TEMP1*A( L+I, J )
                  TEMP2 = TEMP2 + A( L+I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( KPLUS1, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80 J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX = KX
               IY = KY
               L = KPLUS1 - J
               DO 70 I = MAX( 1, J-K ), J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( L+I, J )
                  TEMP2 = TEMP2 + A( L+I, J )*X( IX )
                  IX = IX + INCX
                  IY = IY + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( KPLUS1, J ) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
               IF( J.GT.K ) THEN
                  KX = KX + INCX
                  KY = KY + INCY
               END IF
   80       CONTINUE
         END IF
      ELSE

         // Form  y  when lower triangle of A is stored.

         IF( ( INCX.EQ.1 ) .AND. ( INCY.EQ.1 ) ) THEN
            DO 100 J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               Y( J ) = Y( J ) + TEMP1*A( 1, J )
               L = 1 - J
               DO 90 I = J + 1, MIN( N, J+K )
                  Y( I ) = Y( I ) + TEMP1*A( L+I, J )
                  TEMP2 = TEMP2 + A( L+I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120 J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               Y( JY ) = Y( JY ) + TEMP1*A( 1, J )
               L = 1 - J
               IX = JX
               IY = JY
               DO 110 I = J + 1, MIN( N, J+K )
                  IX = IX + INCX
                  IY = IY + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( L+I, J )
                  TEMP2 = TEMP2 + A( L+I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX = JX + INCX
               JY = JY + INCY
  120       CONTINUE
         END IF
      END IF

      RETURN

      // End of ZSBMV

      END
