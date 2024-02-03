      SUBROUTINE ZLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, OFFSET;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double             VN1( * ), VN2( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      COMPLEX*16         CONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, CONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                I, ITEMP, J, MN, OFFPI, PVT;
      double             TEMP, TEMP2, TOL3Z;
      COMPLEX*16         AII
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARF, ZLARFG, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH, DZNRM2;
      // EXTERNAL IDAMAX, DLAMCH, DZNRM2
      // ..
      // .. Executable Statements ..

      MN = MIN( M-OFFSET, N )
      TOL3Z = SQRT(DLAMCH('Epsilon'))

      // Compute factorization.

      DO 20 I = 1, MN

         OFFPI = OFFSET + I

         // Determine ith pivot column and swap if necessary.

         PVT = ( I-1 ) + IDAMAX( N-I+1, VN1( I ), 1 )

         IF( PVT.NE.I ) THEN
            CALL ZSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( I )
            JPVT( I ) = ITEMP
            VN1( PVT ) = VN1( I )
            VN2( PVT ) = VN2( I )
         END IF

         // Generate elementary reflector H(i).

         IF( OFFPI.LT.M ) THEN
            CALL ZLARFG( M-OFFPI+1, A( OFFPI, I ), A( OFFPI+1, I ), 1, TAU( I ) )
         ELSE
            CALL ZLARFG( 1, A( M, I ), A( M, I ), 1, TAU( I ) )
         END IF

         IF( I.LT.N ) THEN

            // Apply H(i)**H to A(offset+i:m,i+1:n) from the left.

            AII = A( OFFPI, I )
            A( OFFPI, I ) = CONE
            CALL ZLARF( 'Left', M-OFFPI+1, N-I, A( OFFPI, I ), 1, DCONJG( TAU( I ) ), A( OFFPI, I+1 ), LDA, WORK( 1 ) )
            A( OFFPI, I ) = AII
         END IF

         // Update partial column norms.

         DO 10 J = I + 1, N
            IF( VN1( J ).NE.ZERO ) THEN

               // NOTE: The following 4 lines follow from the analysis in
               // Lapack Working Note 176.

               TEMP = ONE - ( ABS( A( OFFPI, J ) ) / VN1( J ) )**2
               TEMP = MAX( TEMP, ZERO )
               TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
               IF( TEMP2 .LE. TOL3Z ) THEN
                  IF( OFFPI.LT.M ) THEN
                     VN1( J ) = DZNRM2( M-OFFPI, A( OFFPI+1, J ), 1 )
                     VN2( J ) = VN1( J )
                  ELSE
                     VN1( J ) = ZERO
                     VN2( J ) = ZERO
                  END IF
               ELSE
                  VN1( J ) = VN1( J )*SQRT( TEMP )
               END IF
            END IF
   10    CONTINUE

   20 CONTINUE

      RETURN

      // End of ZLAQP2

      END
