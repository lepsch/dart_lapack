      SUBROUTINE ZGBT01( M, N, KL, KU, A, LDA, AFAC, LDAFAC, IPIV, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KL, KU, LDA, LDAFAC, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, I1, I2, IL, IP, IW, J, JL, JU, JUA, KD, LENJ;
      double             ANORM, EPS;
      COMPLEX*16         T
      // ..
      // .. External Functions ..
      double             DLAMCH, DZASUM;
      // EXTERNAL DLAMCH, DZASUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick exit if M = 0 or N = 0.

      RESID = ZERO
      if (M <= 0 || N <= 0) RETURN;

      // Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' )
      KD = KU + 1
      ANORM = ZERO
      for (J = 1; J <= N; J++) { // 10
         I1 = MAX( KD+1-J, 1 )
         I2 = MIN( KD+M-J, KL+KD )
         if (I2 >= I1) ANORM = MAX( ANORM, DZASUM( I2-I1+1, A( I1, J ), 1 ) );
      } // 10

      // Compute one column at a time of L*U - A.

      KD = KL + KU + 1
      for (J = 1; J <= N; J++) { // 40

         // Copy the J-th column of U to WORK.

         JU = MIN( KL+KU, J-1 )
         JL = MIN( KL, M-J )
         LENJ = MIN( M, J ) - J + JU + 1
         if ( LENJ > 0 ) {
            zcopy(LENJ, AFAC( KD-JU, J ), 1, WORK, 1 );
            for (I = LENJ + 1; I <= JU + JL + 1; I++) { // 20
               WORK( I ) = ZERO
            } // 20

            // Multiply by the unit lower triangular matrix L.  Note that L
            // is stored as a product of transformations and permutations.

            DO 30 I = MIN( M-1, J ), J - JU, -1
               IL = MIN( KL, M-I )
               if ( IL > 0 ) {
                  IW = I - J + JU + 1
                  T = WORK( IW )
                  zaxpy(IL, T, AFAC( KD+1, I ), 1, WORK( IW+1 ), 1 );
                  IP = IPIV( I )
                  if ( I != IP ) {
                     IP = IP - J + JU + 1
                     WORK( IW ) = WORK( IP )
                     WORK( IP ) = T
                  }
               }
            } // 30

            // Subtract the corresponding column of A.

            JUA = MIN( JU, KU )
            if (JUA+JL+1 > 0) CALL ZAXPY( JUA+JL+1, -DCMPLX( ONE ), A( KU+1-JUA, J ), 1, WORK( JU+1-JUA ), 1 );

            // Compute the 1-norm of the column.

            RESID = MAX( RESID, DZASUM( JU+JL+1, WORK, 1 ) )
         }
      } // 40

      // Compute norm(L*U - A) / ( N * norm(A) * EPS )

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      }

      RETURN

      // End of ZGBT01

      }
