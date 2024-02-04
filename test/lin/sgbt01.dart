      void sgbt01(M, N, KL, KU, A, LDA, AFAC, LDAFAC, IPIV, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KL, KU, LDA, LDAFAC, M, N;
      double               RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double               A( LDA, * ), AFAC( LDAFAC, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, I1, I2, IL, IP, IW, J, JL, JU, JUA, KD, LENJ;
      double               ANORM, EPS, T;
      // ..
      // .. External Functions ..
      //- REAL               SASUM, SLAMCH;
      // EXTERNAL SASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if M = 0 or N = 0.

      RESID = ZERO;
      if (M <= 0 || N <= 0) return;

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      KD = KU + 1;
      ANORM = ZERO;
      for (J = 1; J <= N; J++) { // 10
         I1 = max( KD+1-J, 1 );
         I2 = min( KD+M-J, KL+KD );
         if (I2 >= I1) ANORM = max( ANORM, SASUM( I2-I1+1, A( I1, J ), 1 ) );
      } // 10

      // Compute one column at a time of L*U - A.

      KD = KL + KU + 1;
      for (J = 1; J <= N; J++) { // 40

         // Copy the J-th column of U to WORK.

         JU = min( KL+KU, J-1 );
         JL = min( KL, M-J );
         LENJ = min( M, J ) - J + JU + 1;
         if ( LENJ > 0 ) {
            scopy(LENJ, AFAC( KD-JU, J ), 1, WORK, 1 );
            for (I = LENJ + 1; I <= JU + JL + 1; I++) { // 20
               WORK[I] = ZERO;
            } // 20

            // Multiply by the unit lower triangular matrix L.  Note that L
            // is stored as a product of transformations and permutations.

            for (I = min( M-1, J ); I >= J - JU; I--) { // 30
               IL = min( KL, M-I );
               if ( IL > 0 ) {
                  IW = I - J + JU + 1;
                  T = WORK( IW );
                  saxpy(IL, T, AFAC( KD+1, I ), 1, WORK( IW+1 ), 1 );
                  IP = IPIV( I );
                  if ( I != IP ) {
                     IP = IP - J + JU + 1;
                     WORK[IW] = WORK( IP );
                     WORK[IP] = T;
                  }
               }
            } // 30

            // Subtract the corresponding column of A.

            JUA = min( JU, KU );
            if (JUA+JL+1 > 0) saxpy( JUA+JL+1, -ONE, A( KU+1-JUA, J ), 1, WORK( JU+1-JUA ), 1 );

            // Compute the 1-norm of the column.

            RESID = max( RESID, SASUM( JU+JL+1, WORK, 1 ) );
         }
      } // 40

      // Compute norm(L*U - A) / ( N * norm(A) * EPS )

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;
      }

      return;
      }