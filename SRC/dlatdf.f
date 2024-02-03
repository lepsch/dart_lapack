      SUBROUTINE DLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV, JPIV )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IJOB, LDZ, N
      DOUBLE PRECISION   RDSCAL, RDSUM
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      DOUBLE PRECISION   RHS( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXDIM
      PARAMETER          ( MAXDIM = 8 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, K
      DOUBLE PRECISION   BM, BP, PMONE, SMINU, SPLUS, TEMP
*     ..
*     .. Local Arrays ..
      INTEGER            IWORK( MAXDIM )
      DOUBLE PRECISION   WORK( 4*MAXDIM ), XM( MAXDIM ), XP( MAXDIM )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGECON, DGESC2, DLASSQ, DLASWP, DSCAL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DASUM, DDOT
      EXTERNAL           DASUM, DDOT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( IJOB.NE.2 ) THEN
*
*        Apply permutations IPIV to RHS
*
         CALL DLASWP( 1, RHS, LDZ, 1, N-1, IPIV, 1 )
*
*        Solve for L-part choosing RHS either to +1 or -1.
*
         PMONE = -ONE
*
         DO 10 J = 1, N - 1
            BP = RHS( J ) + ONE
            BM = RHS( J ) - ONE
            SPLUS = ONE
*
*           Look-ahead for L-part RHS(1:N-1) = + or -1, SPLUS and
*           SMIN computed more efficiently than in BSOLVE [1].
*
            SPLUS = SPLUS + DDOT( N-J, Z( J+1, J ), 1, Z( J+1, J ), 1 )
            SMINU = DDOT( N-J, Z( J+1, J ), 1, RHS( J+1 ), 1 )
            SPLUS = SPLUS*RHS( J )
            IF( SPLUS.GT.SMINU ) THEN
               RHS( J ) = BP
            ELSE IF( SMINU.GT.SPLUS ) THEN
               RHS( J ) = BM
            ELSE
*
*              In this case the updating sums are equal and we can
*              choose RHS(J) +1 or -1. The first time this happens
*              we choose -1, thereafter +1. This is a simple way to
*              get good estimates of matrices like Byers well-known
*              example (see [1]). (Not done in BSOLVE.)
*
               RHS( J ) = RHS( J ) + PMONE
               PMONE = ONE
            END IF
*
*           Compute the remaining r.h.s.
*
            TEMP = -RHS( J )
            CALL DAXPY( N-J, TEMP, Z( J+1, J ), 1, RHS( J+1 ), 1 )
*
   10    CONTINUE
*
*        Solve for U-part, look-ahead for RHS(N) = +-1. This is not done
*        in BSOLVE and will hopefully give us a better estimate because
*        any ill-conditioning of the original matrix is transferred to U
*        and not to L. U(N, N) is an approximation to sigma_min(LU).
*
         CALL DCOPY( N-1, RHS, 1, XP, 1 )
         XP( N ) = RHS( N ) + ONE
         RHS( N ) = RHS( N ) - ONE
         SPLUS = ZERO
         SMINU = ZERO
         DO 30 I = N, 1, -1
            TEMP = ONE / Z( I, I )
            XP( I ) = XP( I )*TEMP
            RHS( I ) = RHS( I )*TEMP
            DO 20 K = I + 1, N
               XP( I ) = XP( I ) - XP( K )*( Z( I, K )*TEMP )
               RHS( I ) = RHS( I ) - RHS( K )*( Z( I, K )*TEMP )
   20       CONTINUE
            SPLUS = SPLUS + ABS( XP( I ) )
            SMINU = SMINU + ABS( RHS( I ) )
   30    CONTINUE
         IF( SPLUS.GT.SMINU ) CALL DCOPY( N, XP, 1, RHS, 1 )
*
*        Apply the permutations JPIV to the computed solution (RHS)
*
         CALL DLASWP( 1, RHS, LDZ, 1, N-1, JPIV, -1 )
*
*        Compute the sum of squares
*
         CALL DLASSQ( N, RHS, 1, RDSCAL, RDSUM )
*
      ELSE
*
*        IJOB = 2, Compute approximate nullvector XM of Z
*
         CALL DGECON( 'I', N, Z, LDZ, ONE, TEMP, WORK, IWORK, INFO )
         CALL DCOPY( N, WORK( N+1 ), 1, XM, 1 )
*
*        Compute RHS
*
         CALL DLASWP( 1, XM, LDZ, 1, N-1, IPIV, -1 )
         TEMP = ONE / SQRT( DDOT( N, XM, 1, XM, 1 ) )
         CALL DSCAL( N, TEMP, XM, 1 )
         CALL DCOPY( N, XM, 1, XP, 1 )
         CALL DAXPY( N, ONE, RHS, 1, XP, 1 )
         CALL DAXPY( N, -ONE, XM, 1, RHS, 1 )
         CALL DGESC2( N, Z, LDZ, RHS, IPIV, JPIV, TEMP )
         CALL DGESC2( N, Z, LDZ, XP, IPIV, JPIV, TEMP )
         IF( DASUM( N, XP, 1 ).GT.DASUM( N, RHS, 1 ) ) CALL DCOPY( N, XP, 1, RHS, 1 )
*
*        Compute the sum of squares
*
         CALL DLASSQ( N, RHS, 1, RDSCAL, RDSUM )
*
      END IF
*
      RETURN
*
*     End of DLATDF
*
      END
