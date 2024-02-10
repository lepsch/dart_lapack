      void sgtt01(N, DL, D, DU, DLF, DF, DUF, DU2, final Array<int> IPIV, final Matrix<double> WORK, final int LDWORK, final Array<double> RWORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDWORK, N;
      double               RESID;
      int                IPIV( * );
      double               D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), RWORK( * ), WORK( LDWORK, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, IP, J, LASTJ;
      double               ANORM, EPS, LI;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGT, SLANHS;
      // EXTERNAL SLAMCH, SLANGT, SLANHS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSWAP

      // Quick return if possible

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      EPS = SLAMCH( 'Epsilon' );

      // Copy the matrix U to WORK.

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= N; I++) { // 10
            WORK[I][J] = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= N; I++) { // 30
         if ( I == 1 ) {
            WORK[I][I] = DF( I );
            if (N >= 2) WORK( I, I+1 ) = DUF( I );
            IF( N >= 3 ) WORK( I, I+2 ) = DU2( I );
         } else if ( I == N ) {
            WORK[I][I] = DF( I );
         } else {
            WORK[I][I] = DF( I );
            WORK[I][I+1] = DUF( I );
            if (I < N-1) WORK( I, I+2 ) = DU2( I );
         }
      } // 30

      // Multiply on the left by L.

      LASTJ = N;
      for (I = N - 1; I >= 1; I--) { // 40
         LI = DLF( I );
         saxpy(LASTJ-I+1, LI, WORK( I, I ), LDWORK, WORK( I+1, I ), LDWORK );
         IP = IPIV( I );
         if ( IP == I ) {
            LASTJ = min( I+2, N );
         } else {
            sswap(LASTJ-I+1, WORK( I, I ), LDWORK, WORK( I+1, I ), LDWORK );
         }
      } // 40

      // Subtract the matrix A.

      WORK[1][1] = WORK( 1, 1 ) - D( 1 );
      if ( N > 1 ) {
         WORK[1][2] = WORK( 1, 2 ) - DU( 1 );
         WORK[N][N-1] = WORK( N, N-1 ) - DL( N-1 );
         WORK[N][N] = WORK( N, N ) - D( N );
         for (I = 2; I <= N - 1; I++) { // 50
            WORK[I][I-1] = WORK( I, I-1 ) - DL( I-1 );
            WORK[I][I] = WORK( I, I ) - D( I );
            WORK[I][I+1] = WORK( I, I+1 ) - DU( I );
         } // 50
      }

      // Compute the 1-norm of the tridiagonal matrix A.

      ANORM = SLANGT( '1', N, DL, D, DU );

      // Compute the 1-norm of WORK, which is only guaranteed to be
      // upper Hessenberg.

      RESID = SLANHS( '1', N, WORK, LDWORK, RWORK );

      // Compute norm(L*U - A) / (norm(A) * EPS)

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( RESID / ANORM ) / EPS;
      }

      }
