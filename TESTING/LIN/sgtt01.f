      SUBROUTINE SGTT01( N, DL, D, DU, DLF, DF, DUF, DU2, IPIV, WORK, LDWORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDWORK, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), RWORK( * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IP, J, LASTJ;
      REAL               ANORM, EPS, LI
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGT, SLANHS
      // EXTERNAL SLAMCH, SLANGT, SLANHS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSWAP
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      EPS = SLAMCH( 'Epsilon' )

      // Copy the matrix U to WORK.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= N; I++) { // 10
            WORK( I, J ) = ZERO
         } // 10
      } // 20
      for (I = 1; I <= N; I++) { // 30
         if ( I == 1 ) {
            WORK( I, I ) = DF( I )
            if (N.GE.2) WORK( I, I+1 ) = DUF( I )             IF( N.GE.3 ) WORK( I, I+2 ) = DU2( I );
         } else if ( I == N ) {
            WORK( I, I ) = DF( I )
         } else {
            WORK( I, I ) = DF( I )
            WORK( I, I+1 ) = DUF( I )
            if (I < N-1) WORK( I, I+2 ) = DU2( I );
         }
      } // 30

      // Multiply on the left by L.

      LASTJ = N
      DO 40 I = N - 1, 1, -1
         LI = DLF( I )
         saxpy(LASTJ-I+1, LI, WORK( I, I ), LDWORK, WORK( I+1, I ), LDWORK );
         IP = IPIV( I )
         if ( IP == I ) {
            LASTJ = MIN( I+2, N )
         } else {
            sswap(LASTJ-I+1, WORK( I, I ), LDWORK, WORK( I+1, I ), LDWORK );
         }
      } // 40

      // Subtract the matrix A.

      WORK( 1, 1 ) = WORK( 1, 1 ) - D( 1 )
      if ( N > 1 ) {
         WORK( 1, 2 ) = WORK( 1, 2 ) - DU( 1 )
         WORK( N, N-1 ) = WORK( N, N-1 ) - DL( N-1 )
         WORK( N, N ) = WORK( N, N ) - D( N )
         for (I = 2; I <= N - 1; I++) { // 50
            WORK( I, I-1 ) = WORK( I, I-1 ) - DL( I-1 )
            WORK( I, I ) = WORK( I, I ) - D( I )
            WORK( I, I+1 ) = WORK( I, I+1 ) - DU( I )
         } // 50
      }

      // Compute the 1-norm of the tridiagonal matrix A.

      ANORM = SLANGT( '1', N, DL, D, DU )

      // Compute the 1-norm of WORK, which is only guaranteed to be
      // upper Hessenberg.

      RESID = SLANHS( '1', N, WORK, LDWORK, RWORK )

      // Compute norm(L*U - A) / (norm(A) * EPS)

      if ( ANORM.LE.ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( RESID / ANORM ) / EPS
      }

      RETURN

      // End of SGTT01

      }
