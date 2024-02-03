      SUBROUTINE SLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMBDA, Q2, INDX, CTOT, W, S, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDQ, N, N1;
      REAL               RHO
      // ..
      // .. Array Arguments ..
      int                CTOT( * ), INDX( * );
      REAL               D( * ), DLAMBDA( * ), Q( LDQ, * ), Q2( * ), S( * ), W( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E0, ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I, II, IQ2, J, N12, N2, N23;
      REAL               TEMP
      // ..
      // .. External Functions ..
      REAL               SNRM2
      // EXTERNAL SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMM, SLACPY, SLAED4, SLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.K ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAED3', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( K.EQ.0 ) RETURN

      DO 20 J = 1, K
         CALL SLAED4( K, J, DLAMBDA, W, Q( 1, J ), RHO, D( J ), INFO )

         // If the zero finder fails, the computation is terminated.

         IF( INFO.NE.0 ) GO TO 120
   20 CONTINUE

      IF( K.EQ.1 ) GO TO 110
      IF( K.EQ.2 ) THEN
         DO 30 J = 1, K
            W( 1 ) = Q( 1, J )
            W( 2 ) = Q( 2, J )
            II = INDX( 1 )
            Q( 1, J ) = W( II )
            II = INDX( 2 )
            Q( 2, J ) = W( II )
   30    CONTINUE
         GO TO 110
      END IF

      // Compute updated W.

      CALL SCOPY( K, W, 1, S, 1 )

      // Initialize W(I) = Q(I,I)

      CALL SCOPY( K, Q, LDQ+1, W, 1 )
      DO 60 J = 1, K
         DO 40 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) )
   40    CONTINUE
         DO 50 I = J + 1, K
            W( I ) = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) )
   50    CONTINUE
   60 CONTINUE
      DO 70 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I ) )
   70 CONTINUE

      // Compute eigenvectors of the modified rank-1 modification.

      DO 100 J = 1, K
         DO 80 I = 1, K
            S( I ) = W( I ) / Q( I, J )
   80    CONTINUE
         TEMP = SNRM2( K, S, 1 )
         DO 90 I = 1, K
            II = INDX( I )
            Q( I, J ) = S( II ) / TEMP
   90    CONTINUE
  100 CONTINUE

      // Compute the updated eigenvectors.

  110 CONTINUE

      N2 = N - N1
      N12 = CTOT( 1 ) + CTOT( 2 )
      N23 = CTOT( 2 ) + CTOT( 3 )

      CALL SLACPY( 'A', N23, K, Q( CTOT( 1 )+1, 1 ), LDQ, S, N23 )
      IQ2 = N1*N12 + 1
      IF( N23.NE.0 ) THEN
         CALL SGEMM( 'N', 'N', N2, K, N23, ONE, Q2( IQ2 ), N2, S, N23, ZERO, Q( N1+1, 1 ), LDQ )
      ELSE
         CALL SLASET( 'A', N2, K, ZERO, ZERO, Q( N1+1, 1 ), LDQ )
      END IF

      CALL SLACPY( 'A', N12, K, Q, LDQ, S, N12 )
      IF( N12.NE.0 ) THEN
         CALL SGEMM( 'N', 'N', N1, K, N12, ONE, Q2, N1, S, N12, ZERO, Q, LDQ )
      ELSE
         CALL SLASET( 'A', N1, K, ZERO, ZERO, Q( 1, 1 ), LDQ )
      END IF


  120 CONTINUE
      RETURN

      // End of SLAED3

      }
