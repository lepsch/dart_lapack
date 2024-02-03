      SUBROUTINE DORT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             ROWCOL;
      int                LDU, LWORK, M, N
      double             RESID;
*     ..
*     .. Array Arguments ..
      double             U( LDU, * ), WORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      String             TRANSU;
      int                I, J, K, LDWORK, MNMIN
      double             EPS, TMP;
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DDOT, DLAMCH, DLANSY;
      EXTERNAL           LSAME, DDOT, DLAMCH, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASET, DSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RESID = ZERO
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
      EPS = DLAMCH( 'Precision' )
      IF( M.LT.N .OR. ( M.EQ.N .AND. LSAME( ROWCOL, 'R' ) ) ) THEN
         TRANSU = 'N'
         K = N
      ELSE
         TRANSU = 'T'
         K = M
      END IF
      MNMIN = MIN( M, N )
*
      IF( ( MNMIN+1 )*MNMIN.LE.LWORK ) THEN
         LDWORK = MNMIN
      ELSE
         LDWORK = 0
      END IF
      IF( LDWORK.GT.0 ) THEN
*
*        Compute I - U*U' or I - U'*U.
*
         CALL DLASET( 'Upper', MNMIN, MNMIN, ZERO, ONE, WORK, LDWORK )
         CALL DSYRK( 'Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK, LDWORK )
*
*        Compute norm( I - U*U' ) / ( K * EPS ) .
*
         RESID = DLANSY( '1', 'Upper', MNMIN, WORK, LDWORK, WORK( LDWORK*MNMIN+1 ) )
         RESID = ( RESID / DBLE( K ) ) / EPS
      ELSE IF( TRANSU.EQ.'T' ) THEN
*
*        Find the maximum element in abs( I - U'*U ) / ( m * EPS )
*
         DO 20 J = 1, N
            DO 10 I = 1, J
               IF( I.NE.J ) THEN
                  TMP = ZERO
               ELSE
                  TMP = ONE
               END IF
               TMP = TMP - DDOT( M, U( 1, I ), 1, U( 1, J ), 1 )
               RESID = MAX( RESID, ABS( TMP ) )
   10       CONTINUE
   20    CONTINUE
         RESID = ( RESID / DBLE( M ) ) / EPS
      ELSE
*
*        Find the maximum element in abs( I - U*U' ) / ( n * EPS )
*
         DO 40 J = 1, M
            DO 30 I = 1, J
               IF( I.NE.J ) THEN
                  TMP = ZERO
               ELSE
                  TMP = ONE
               END IF
               TMP = TMP - DDOT( N, U( J, 1 ), LDU, U( I, 1 ), LDU )
               RESID = MAX( RESID, ABS( TMP ) )
   30       CONTINUE
   40    CONTINUE
         RESID = ( RESID / DBLE( N ) ) / EPS
      END IF
      RETURN
*
*     End of DORT01
*
      END
