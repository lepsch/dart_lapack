      SUBROUTINE SLASQ1( N, D, E, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, N;
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      int                I, IINFO;
      REAL               EPS, SCALE, SAFMIN, SIGMN, SIGMX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SLAS2, SLASCL, SLASQ2, SLASRT, XERBLA
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'SLASQ1', -INFO )
         RETURN
      ELSE IF( N.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      ELSE IF( N.EQ.2 ) THEN
         CALL SLAS2( D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX )
         D( 1 ) = SIGMX
         D( 2 ) = SIGMN
         RETURN
      END IF
*
*     Estimate the largest singular value.
*
      SIGMX = ZERO
      DO 10 I = 1, N - 1
         D( I ) = ABS( D( I ) )
         SIGMX = MAX( SIGMX, ABS( E( I ) ) )
   10 CONTINUE
      D( N ) = ABS( D( N ) )
*
*     Early return if SIGMX is zero (matrix is already diagonal).
*
      IF( SIGMX.EQ.ZERO ) THEN
         CALL SLASRT( 'D', N, D, IINFO )
         RETURN
      END IF
*
      DO 20 I = 1, N
         SIGMX = MAX( SIGMX, D( I ) )
   20 CONTINUE
*
*     Copy D and E into WORK (in the Z format) and scale (squaring the
*     input data makes scaling by a power of the radix pointless).
*
      EPS = SLAMCH( 'Precision' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      SCALE = SQRT( EPS / SAFMIN )
      CALL SCOPY( N, D, 1, WORK( 1 ), 2 )
      CALL SCOPY( N-1, E, 1, WORK( 2 ), 2 )
      CALL SLASCL( 'G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1, IINFO )
*
*     Compute the q's and e's.
*
      DO 30 I = 1, 2*N - 1
         WORK( I ) = WORK( I )**2
   30 CONTINUE
      WORK( 2*N ) = ZERO
*
      CALL SLASQ2( N, WORK, INFO )
*
      IF( INFO.EQ.0 ) THEN
         DO 40 I = 1, N
            D( I ) = SQRT( WORK( I ) )
   40    CONTINUE
         CALL SLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
      ELSE IF( INFO.EQ.2 ) THEN
*
*     Maximum number of iterations exceeded.  Move data from WORK
*     into D and E so the calling subroutine can try to finish
*
         DO I = 1, N
            D( I ) = SQRT( WORK( 2*I-1 ) )
            E( I ) = SQRT( WORK( 2*I ) )
         END DO
         CALL SLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
         CALL SLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, E, N, IINFO )
      END IF
*
      RETURN
*
*     End of SLASQ1
*
      END
