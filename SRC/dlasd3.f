      SUBROUTINE DLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2, LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDQ, LDU, LDU2, LDVT, LDVT2, NL, NR, SQRE;
      // ..
      // .. Array Arguments ..
      int                CTOT( * ), IDXC( * );
      double             D( * ), DSIGMA( * ), Q( LDQ, * ), U( LDU, * ), U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ), Z( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO, NEGONE;
      const              ONE = 1.0D+0, ZERO = 0.0D+0, NEGONE = -1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                CTEMP, I, J, JC, KTEMP, M, N, NLP1, NLP2, NRP1;
      double             RHO, TEMP;
      // ..
      // .. External Functions ..
      double             DNRM2;
      // EXTERNAL DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DLACPY, DLASCL, DLASD4, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( NL.LT.1 ) {
         INFO = -1
      } else if ( NR.LT.1 ) {
         INFO = -2
      } else if ( ( SQRE.NE.1 ) .AND. ( SQRE.NE.0 ) ) {
         INFO = -3
      }

      N = NL + NR + 1
      M = N + SQRE
      NLP1 = NL + 1
      NLP2 = NL + 2

      if ( ( K.LT.1 ) .OR. ( K.GT.N ) ) {
         INFO = -4
      } else if ( LDQ.LT.K ) {
         INFO = -7
      } else if ( LDU.LT.N ) {
         INFO = -10
      } else if ( LDU2.LT.N ) {
         INFO = -12
      } else if ( LDVT.LT.M ) {
         INFO = -14
      } else if ( LDVT2.LT.M ) {
         INFO = -16
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DLASD3', -INFO )
         RETURN
      }

      // Quick return if possible

      if ( K.EQ.1 ) {
         D( 1 ) = ABS( Z( 1 ) )
         CALL DCOPY( M, VT2( 1, 1 ), LDVT2, VT( 1, 1 ), LDVT )
         if ( Z( 1 ).GT.ZERO ) {
            CALL DCOPY( N, U2( 1, 1 ), 1, U( 1, 1 ), 1 )
         } else {
            DO 10 I = 1, N
               U( I, 1 ) = -U2( I, 1 )
   10       CONTINUE
         }
         RETURN
      }

      // Keep a copy of Z.

      CALL DCOPY( K, Z, 1, Q, 1 )

      // Normalize Z.

      RHO = DNRM2( K, Z, 1 )
      CALL DLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO

      // Find the new singular values.

      DO 30 J = 1, K
         CALL DLASD4( K, J, DSIGMA, Z, U( 1, J ), RHO, D( J ), VT( 1, J ), INFO )

         // If the zero finder fails, report the convergence failure.

         if ( INFO.NE.0 ) {
            RETURN
         }
   30 CONTINUE

      // Compute updated Z.

      DO 60 I = 1, K
         Z( I ) = U( I, K )*VT( I, K )
         DO 40 J = 1, I - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / ( DSIGMA( I )-DSIGMA( J ) ) / ( DSIGMA( I )+DSIGMA( J ) ) )
   40    CONTINUE
         DO 50 J = I, K - 1
            Z( I ) = Z( I )*( U( I, J )*VT( I, J ) / ( DSIGMA( I )-DSIGMA( J+1 ) ) / ( DSIGMA( I )+DSIGMA( J+1 ) ) )
   50    CONTINUE
         Z( I ) = SIGN( SQRT( ABS( Z( I ) ) ), Q( I, 1 ) )
   60 CONTINUE

      // Compute left singular vectors of the modified diagonal matrix,
      // and store related information for the right singular vectors.

      DO 90 I = 1, K
         VT( 1, I ) = Z( 1 ) / U( 1, I ) / VT( 1, I )
         U( 1, I ) = NEGONE
         DO 70 J = 2, K
            VT( J, I ) = Z( J ) / U( J, I ) / VT( J, I )
            U( J, I ) = DSIGMA( J )*VT( J, I )
   70    CONTINUE
         TEMP = DNRM2( K, U( 1, I ), 1 )
         Q( 1, I ) = U( 1, I ) / TEMP
         DO 80 J = 2, K
            JC = IDXC( J )
            Q( J, I ) = U( JC, I ) / TEMP
   80    CONTINUE
   90 CONTINUE

      // Update the left singular vector matrix.

      if ( K.EQ.2 ) {
         CALL DGEMM( 'N', 'N', N, K, K, ONE, U2, LDU2, Q, LDQ, ZERO, U, LDU )
         GO TO 100
      }
      if ( CTOT( 1 ).GT.0 ) {
         CALL DGEMM( 'N', 'N', NL, K, CTOT( 1 ), ONE, U2( 1, 2 ), LDU2, Q( 2, 1 ), LDQ, ZERO, U( 1, 1 ), LDU )
         if ( CTOT( 3 ).GT.0 ) {
            KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
            CALL DGEMM( 'N', 'N', NL, K, CTOT( 3 ), ONE, U2( 1, KTEMP ), LDU2, Q( KTEMP, 1 ), LDQ, ONE, U( 1, 1 ), LDU )
         }
      } else if ( CTOT( 3 ).GT.0 ) {
         KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
         CALL DGEMM( 'N', 'N', NL, K, CTOT( 3 ), ONE, U2( 1, KTEMP ), LDU2, Q( KTEMP, 1 ), LDQ, ZERO, U( 1, 1 ), LDU )
      } else {
         CALL DLACPY( 'F', NL, K, U2, LDU2, U, LDU )
      }
      CALL DCOPY( K, Q( 1, 1 ), LDQ, U( NLP1, 1 ), LDU )
      KTEMP = 2 + CTOT( 1 )
      CTEMP = CTOT( 2 ) + CTOT( 3 )
      CALL DGEMM( 'N', 'N', NR, K, CTEMP, ONE, U2( NLP2, KTEMP ), LDU2, Q( KTEMP, 1 ), LDQ, ZERO, U( NLP2, 1 ), LDU )

      // Generate the right singular vectors.

  100 CONTINUE
      DO 120 I = 1, K
         TEMP = DNRM2( K, VT( 1, I ), 1 )
         Q( I, 1 ) = VT( 1, I ) / TEMP
         DO 110 J = 2, K
            JC = IDXC( J )
            Q( I, J ) = VT( JC, I ) / TEMP
  110    CONTINUE
  120 CONTINUE

      // Update the right singular vector matrix.

      if ( K.EQ.2 ) {
         CALL DGEMM( 'N', 'N', K, M, K, ONE, Q, LDQ, VT2, LDVT2, ZERO, VT, LDVT )
         RETURN
      }
      KTEMP = 1 + CTOT( 1 )
      CALL DGEMM( 'N', 'N', K, NLP1, KTEMP, ONE, Q( 1, 1 ), LDQ, VT2( 1, 1 ), LDVT2, ZERO, VT( 1, 1 ), LDVT )
      KTEMP = 2 + CTOT( 1 ) + CTOT( 2 )
      IF( KTEMP.LE.LDVT2 ) CALL DGEMM( 'N', 'N', K, NLP1, CTOT( 3 ), ONE, Q( 1, KTEMP ), LDQ, VT2( KTEMP, 1 ), LDVT2, ONE, VT( 1, 1 ), LDVT )

      KTEMP = CTOT( 1 ) + 1
      NRP1 = NR + SQRE
      if ( KTEMP.GT.1 ) {
         DO 130 I = 1, K
            Q( I, KTEMP ) = Q( I, 1 )
  130    CONTINUE
         DO 140 I = NLP2, M
            VT2( KTEMP, I ) = VT2( 1, I )
  140    CONTINUE
      }
      CTEMP = 1 + CTOT( 2 ) + CTOT( 3 )
      CALL DGEMM( 'N', 'N', K, NRP1, CTEMP, ONE, Q( 1, KTEMP ), LDQ, VT2( KTEMP, NLP2 ), LDVT2, ZERO, VT( 1, NLP2 ), LDVT )

      RETURN

      // End of DLASD3

      }
