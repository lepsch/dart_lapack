      SUBROUTINE SGBT01( M, N, KL, KU, A, LDA, AFAC, LDAFAC, IPIV, WORK,
     $                   RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            KL, KU, LDA, LDAFAC, M, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), AFAC( LDAFAC, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, I2, IL, IP, IW, J, JL, JU, JUA, KD, LENJ
      REAL               ANORM, EPS, T
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH
      EXTERNAL           SASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick exit if M = 0 or N = 0.
*
      RESID = ZERO
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*     Determine EPS and the norm of A.
*
      EPS = SLAMCH( 'Epsilon' )
      KD = KU + 1
      ANORM = ZERO
      DO 10 J = 1, N
         I1 = MAX( KD+1-J, 1 )
         I2 = MIN( KD+M-J, KL+KD )
         IF( I2.GE.I1 )
     $      ANORM = MAX( ANORM, SASUM( I2-I1+1, A( I1, J ), 1 ) )
   10 CONTINUE
*
*     Compute one column at a time of L*U - A.
*
      KD = KL + KU + 1
      DO 40 J = 1, N
*
*        Copy the J-th column of U to WORK.
*
         JU = MIN( KL+KU, J-1 )
         JL = MIN( KL, M-J )
         LENJ = MIN( M, J ) - J + JU + 1
         IF( LENJ.GT.0 ) THEN
            CALL SCOPY( LENJ, AFAC( KD-JU, J ), 1, WORK, 1 )
            DO 20 I = LENJ + 1, JU + JL + 1
               WORK( I ) = ZERO
   20       CONTINUE
*
*           Multiply by the unit lower triangular matrix L.  Note that L
*           is stored as a product of transformations and permutations.
*
            DO 30 I = MIN( M-1, J ), J - JU, -1
               IL = MIN( KL, M-I )
               IF( IL.GT.0 ) THEN
                  IW = I - J + JU + 1
                  T = WORK( IW )
                  CALL SAXPY( IL, T, AFAC( KD+1, I ), 1, WORK( IW+1 ),
     $                        1 )
                  IP = IPIV( I )
                  IF( I.NE.IP ) THEN
                     IP = IP - J + JU + 1
                     WORK( IW ) = WORK( IP )
                     WORK( IP ) = T
                  END IF
               END IF
   30       CONTINUE
*
*           Subtract the corresponding column of A.
*
            JUA = MIN( JU, KU )
            IF( JUA+JL+1.GT.0 )
     $         CALL SAXPY( JUA+JL+1, -ONE, A( KU+1-JUA, J ), 1,
     $                     WORK( JU+1-JUA ), 1 )
*
*           Compute the 1-norm of the column.
*
            RESID = MAX( RESID, SASUM( JU+JL+1, WORK, 1 ) )
         END IF
   40 CONTINUE
*
*     Compute norm(L*U - A) / ( N * norm(A) * EPS )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF
*
      RETURN
*
*     End of SGBT01
*
      END