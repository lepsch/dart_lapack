      SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMBDA, W, S, LDS, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, KSTART, KSTOP, LDQ, LDS, N;
      double             RHO;
      // ..
      // .. Array Arguments ..
      double             D( * ), DLAMBDA( * ), Q( LDQ, * ), S( LDS, * ), W( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             TEMP;
      // ..
      // .. External Functions ..
      double             DNRM2;
      // EXTERNAL DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAED4, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( KSTART.LT.1 .OR. KSTART.GT.MAX( 1, K ) ) THEN
         INFO = -2
      ELSE IF( MAX( 1, KSTOP ).LT.KSTART .OR. KSTOP.GT.MAX( 1, K ) ) THEN
         INFO = -3
      ELSE IF( N.LT.K ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDS.LT.MAX( 1, K ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED9', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( K.EQ.0 ) RETURN

      DO 20 J = KSTART, KSTOP
         CALL DLAED4( K, J, DLAMBDA, W, Q( 1, J ), RHO, D( J ), INFO )

         // If the zero finder fails, the computation is terminated.

         IF( INFO.NE.0 ) GO TO 120
   20 CONTINUE

      IF( K.EQ.1 .OR. K.EQ.2 ) THEN
         DO 40 I = 1, K
            DO 30 J = 1, K
               S( J, I ) = Q( J, I )
   30       CONTINUE
   40    CONTINUE
         GO TO 120
      END IF

      // Compute updated W.

      CALL DCOPY( K, W, 1, S, 1 )

      // Initialize W(I) = Q(I,I)

      CALL DCOPY( K, Q, LDQ+1, W, 1 )
      DO 70 J = 1, K
         DO 50 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) )
   50    CONTINUE
         DO 60 I = J + 1, K
            W( I ) = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) )
   60    CONTINUE
   70 CONTINUE
      DO 80 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I, 1 ) )
   80 CONTINUE

      // Compute eigenvectors of the modified rank-1 modification.

      DO 110 J = 1, K
         DO 90 I = 1, K
            Q( I, J ) = W( I ) / Q( I, J )
   90    CONTINUE
         TEMP = DNRM2( K, Q( 1, J ), 1 )
         DO 100 I = 1, K
            S( I, J ) = Q( I, J ) / TEMP
  100    CONTINUE
  110 CONTINUE

  120 CONTINUE
      RETURN

      // End of DLAED9

      }
