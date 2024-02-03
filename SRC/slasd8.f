      SUBROUTINE SLASD8( ICOMPQ, K, D, Z, VF, VL, DIFL, DIFR, LDDIFR, DSIGMA, WORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            ICOMPQ, INFO, K, LDDIFR
*     ..
*     .. Array Arguments ..
      REAL               D( * ), DIFL( * ), DIFR( LDDIFR, * ), DSIGMA( * ), VF( * ), VL( * ), WORK( * ), Z( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J
      REAL               DIFLJ, DIFRJ, DJ, DSIGJ, DSIGJP, RHO, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLASCL, SLASD4, SLASET, XERBLA
*     ..
*     .. External Functions ..
      REAL               SDOT, SLAMC3, SNRM2
      EXTERNAL           SDOT, SLAMC3, SNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( ( ICOMPQ.LT.0 ) .OR. ( ICOMPQ.GT.1 ) ) THEN
         INFO = -1
      ELSE IF( K.LT.1 ) THEN
         INFO = -2
      ELSE IF( LDDIFR.LT.K ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASD8', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( K.EQ.1 ) THEN
         D( 1 ) = ABS( Z( 1 ) )
         DIFL( 1 ) = D( 1 )
         IF( ICOMPQ.EQ.1 ) THEN
            DIFL( 2 ) = ONE
            DIFR( 1, 2 ) = ONE
         END IF
         RETURN
      END IF
*
*     Book keeping.
*
      IWK1 = 1
      IWK2 = IWK1 + K
      IWK3 = IWK2 + K
      IWK2I = IWK2 - 1
      IWK3I = IWK3 - 1
*
*     Normalize Z.
*
      RHO = SNRM2( K, Z, 1 )
      CALL SLASCL( 'G', 0, 0, RHO, ONE, K, 1, Z, K, INFO )
      RHO = RHO*RHO
*
*     Initialize WORK(IWK3).
*
      CALL SLASET( 'A', K, 1, ONE, ONE, WORK( IWK3 ), K )
*
*     Compute the updated singular values, the arrays DIFL, DIFR,
*     and the updated Z.
*
      DO 40 J = 1, K
         CALL SLASD4( K, J, DSIGMA, Z, WORK( IWK1 ), RHO, D( J ), WORK( IWK2 ), INFO )
*
*        If the root finder fails, report the convergence failure.
*
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
         WORK( IWK3I+J ) = WORK( IWK3I+J )*WORK( J )*WORK( IWK2I+J )
         DIFL( J ) = -WORK( J )
         DIFR( J, 1 ) = -WORK( J+1 )
         DO 20 I = 1, J - 1
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* WORK( IWK2I+I ) / ( DSIGMA( I )- DSIGMA( J ) ) / ( DSIGMA( I )+ DSIGMA( J ) )
   20    CONTINUE
         DO 30 I = J + 1, K
            WORK( IWK3I+I ) = WORK( IWK3I+I )*WORK( I )* WORK( IWK2I+I ) / ( DSIGMA( I )- DSIGMA( J ) ) / ( DSIGMA( I )+ DSIGMA( J ) )
   30    CONTINUE
   40 CONTINUE
*
*     Compute updated Z.
*
      DO 50 I = 1, K
         Z( I ) = SIGN( SQRT( ABS( WORK( IWK3I+I ) ) ), Z( I ) )
   50 CONTINUE
*
*     Update VF and VL.
*
      DO 80 J = 1, K
         DIFLJ = DIFL( J )
         DJ = D( J )
         DSIGJ = -DSIGMA( J )
         IF( J.LT.K ) THEN
            DIFRJ = -DIFR( J, 1 )
            DSIGJP = -DSIGMA( J+1 )
         END IF
         WORK( J ) = -Z( J ) / DIFLJ / ( DSIGMA( J )+DJ )
*
*        Use calls to the subroutine SLAMC3 to enforce the parentheses
*        (x+y)+z. The goal is to prevent optimizing compilers
*        from doing x+(y+z).
*
         DO 60 I = 1, J - 1
            WORK( I ) = Z( I ) / ( SLAMC3( DSIGMA( I ), DSIGJ )-DIFLJ ) / ( DSIGMA( I )+DJ )
   60    CONTINUE
         DO 70 I = J + 1, K
            WORK( I ) = Z( I ) / ( SLAMC3( DSIGMA( I ), DSIGJP )+DIFRJ ) / ( DSIGMA( I )+DJ )
   70    CONTINUE
         TEMP = SNRM2( K, WORK, 1 )
         WORK( IWK2I+J ) = SDOT( K, WORK, 1, VF, 1 ) / TEMP
         WORK( IWK3I+J ) = SDOT( K, WORK, 1, VL, 1 ) / TEMP
         IF( ICOMPQ.EQ.1 ) THEN
            DIFR( J, 2 ) = TEMP
         END IF
   80 CONTINUE
*
      CALL SCOPY( K, WORK( IWK2 ), 1, VF, 1 )
      CALL SCOPY( K, WORK( IWK3 ), 1, VL, 1 )
*
      RETURN
*
*     End of SLASD8
*
      END
