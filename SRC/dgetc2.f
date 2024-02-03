      SUBROUTINE DGETC2( N, A, LDA, IPIV, JPIV, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      int                I, IP, IPV, J, JP, JPV;
      double             BIGNUM, EPS, SMIN, SMLNUM, XMAX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGER, DSWAP
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Set constants to control overflow

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // Handle the case N=1 by itself

      IF( N.EQ.1 ) THEN
         IPIV( 1 ) = 1
         JPIV( 1 ) = 1
         IF( ABS( A( 1, 1 ) ).LT.SMLNUM ) THEN
            INFO = 1
            A( 1, 1 ) = SMLNUM
         END IF
         RETURN
      END IF

      // Factorize A using complete pivoting.
      // Set pivots less than SMIN to SMIN.

      DO 40 I = 1, N - 1

         // Find max element in matrix A

         XMAX = ZERO
         DO 20 IP = I, N
            DO 10 JP = I, N
               IF( ABS( A( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( A( IP, JP ) )
                  IPV = IP
                  JPV = JP
               END IF
   10       CONTINUE
   20    CONTINUE
         IF( I.EQ.1 ) SMIN = MAX( EPS*XMAX, SMLNUM )

         // Swap rows

         IF( IPV.NE.I ) CALL DSWAP( N, A( IPV, 1 ), LDA, A( I, 1 ), LDA )
         IPIV( I ) = IPV

         // Swap columns

         IF( JPV.NE.I ) CALL DSWAP( N, A( 1, JPV ), 1, A( 1, I ), 1 )
         JPIV( I ) = JPV

         // Check for singularity

         IF( ABS( A( I, I ) ).LT.SMIN ) THEN
            INFO = I
            A( I, I ) = SMIN
         END IF
         DO 30 J = I + 1, N
            A( J, I ) = A( J, I ) / A( I, I )
   30    CONTINUE
         CALL DGER( N-I, N-I, -ONE, A( I+1, I ), 1, A( I, I+1 ), LDA, A( I+1, I+1 ), LDA )
   40 CONTINUE

      IF( ABS( A( N, N ) ).LT.SMIN ) THEN
         INFO = N
         A( N, N ) = SMIN
      END IF

      // Set last pivots to N

      IPIV( N ) = N
      JPIV( N ) = N

      RETURN

      // End of DGETC2

      END
