      SUBROUTINE DGTT01( N, DL, D, DU, DLF, DF, DUF, DU2, IPIV, WORK, LDWORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDWORK, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), RWORK( * ), WORK( LDWORK, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IP, J, LASTJ;
      double             ANORM, EPS, LI;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGT, DLANHS;
      // EXTERNAL DLAMCH, DLANGT, DLANHS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DSWAP
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      EPS = DLAMCH( 'Epsilon' )

      // Copy the matrix U to WORK.

      DO 20 J = 1, N
         DO 10 I = 1, N
            WORK( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, N
         if ( I.EQ.1 ) {
            WORK( I, I ) = DF( I )
            IF( N.GE.2 ) WORK( I, I+1 ) = DUF( I )             IF( N.GE.3 ) WORK( I, I+2 ) = DU2( I )
         } else if ( I.EQ.N ) {
            WORK( I, I ) = DF( I )
         } else {
            WORK( I, I ) = DF( I )
            WORK( I, I+1 ) = DUF( I )
            IF( I.LT.N-1 ) WORK( I, I+2 ) = DU2( I )
         }
   30 CONTINUE

      // Multiply on the left by L.

      LASTJ = N
      DO 40 I = N - 1, 1, -1
         LI = DLF( I )
         CALL DAXPY( LASTJ-I+1, LI, WORK( I, I ), LDWORK, WORK( I+1, I ), LDWORK )
         IP = IPIV( I )
         if ( IP.EQ.I ) {
            LASTJ = MIN( I+2, N )
         } else {
            CALL DSWAP( LASTJ-I+1, WORK( I, I ), LDWORK, WORK( I+1, I ), LDWORK )
         }
   40 CONTINUE

      // Subtract the matrix A.

      WORK( 1, 1 ) = WORK( 1, 1 ) - D( 1 )
      if ( N.GT.1 ) {
         WORK( 1, 2 ) = WORK( 1, 2 ) - DU( 1 )
         WORK( N, N-1 ) = WORK( N, N-1 ) - DL( N-1 )
         WORK( N, N ) = WORK( N, N ) - D( N )
         DO 50 I = 2, N - 1
            WORK( I, I-1 ) = WORK( I, I-1 ) - DL( I-1 )
            WORK( I, I ) = WORK( I, I ) - D( I )
            WORK( I, I+1 ) = WORK( I, I+1 ) - DU( I )
   50    CONTINUE
      }

      // Compute the 1-norm of the tridiagonal matrix A.

      ANORM = DLANGT( '1', N, DL, D, DU )

      // Compute the 1-norm of WORK, which is only guaranteed to be
      // upper Hessenberg.

      RESID = DLANHS( '1', N, WORK, LDWORK, RWORK )

      // Compute norm(L*U - A) / (norm(A) * EPS)

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( RESID / ANORM ) / EPS
      }

      RETURN

      // End of DGTT01

      }
