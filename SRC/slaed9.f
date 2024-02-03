      SUBROUTINE SLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMBDA, W, S, LDS, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, KSTART, KSTOP, LDQ, LDS, N;
      REAL               RHO
      // ..
      // .. Array Arguments ..
      REAL               D( * ), DLAMBDA( * ), Q( LDQ, * ), S( LDS, * ), W( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      REAL               TEMP
      // ..
      // .. External Functions ..
      REAL               SNRM2
      // EXTERNAL SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLAED4, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( K.LT.0 ) {
         INFO = -1
      } else if ( KSTART.LT.1 .OR. KSTART.GT.MAX( 1, K ) ) {
         INFO = -2
      } else if ( MAX( 1, KSTOP ).LT.KSTART .OR. KSTOP.GT.MAX( 1, K ) ) {
         INFO = -3
      } else if ( N.LT.K ) {
         INFO = -4
      } else if ( LDQ.LT.MAX( 1, K ) ) {
         INFO = -7
      } else if ( LDS.LT.MAX( 1, K ) ) {
         INFO = -12
      }
      if ( INFO.NE.0 ) {
         xerbla('SLAED9', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( K.EQ.0 ) RETURN

      DO 20 J = KSTART, KSTOP
         slaed4(K, J, DLAMBDA, W, Q( 1, J ), RHO, D( J ), INFO );

         // If the zero finder fails, the computation is terminated.

         IF( INFO.NE.0 ) GO TO 120
   20 CONTINUE

      if ( K.EQ.1 .OR. K.EQ.2 ) {
         DO 40 I = 1, K
            DO 30 J = 1, K
               S( J, I ) = Q( J, I )
   30       CONTINUE
   40    CONTINUE
         GO TO 120
      }

      // Compute updated W.

      scopy(K, W, 1, S, 1 );

      // Initialize W(I) = Q(I,I)

      scopy(K, Q, LDQ+1, W, 1 );
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
         TEMP = SNRM2( K, Q( 1, J ), 1 )
         DO 100 I = 1, K
            S( I, J ) = Q( I, J ) / TEMP
  100    CONTINUE
  110 CONTINUE

  120 CONTINUE
      RETURN

      // End of SLAED9

      }
