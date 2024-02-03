      SUBROUTINE DLASQ1( N, D, E, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, IINFO;
      double             EPS, SCALE, SAFMIN, SIGMN, SIGMX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAS2, DLASCL, DLASQ2, DLASRT, XERBLA
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
         CALL XERBLA( 'DLASQ1', -INFO )
         RETURN
      } else if ( N.EQ.0 ) {
         RETURN
      } else if ( N.EQ.1 ) {
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      } else if ( N.EQ.2 ) {
         CALL DLAS2( D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX )
         D( 1 ) = SIGMX
         D( 2 ) = SIGMN
         RETURN
      }

      // Estimate the largest singular value.

      SIGMX = ZERO
      DO 10 I = 1, N - 1
         D( I ) = ABS( D( I ) )
         SIGMX = MAX( SIGMX, ABS( E( I ) ) )
   10 CONTINUE
      D( N ) = ABS( D( N ) )

      // Early return if SIGMX is zero (matrix is already diagonal).

      if ( SIGMX.EQ.ZERO ) {
         CALL DLASRT( 'D', N, D, IINFO )
         RETURN
      }

      DO 20 I = 1, N
         SIGMX = MAX( SIGMX, D( I ) )
   20 CONTINUE

      // Copy D and E into WORK (in the Z format) and scale (squaring the
      // input data makes scaling by a power of the radix pointless).

      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SCALE = SQRT( EPS / SAFMIN )
      CALL DCOPY( N, D, 1, WORK( 1 ), 2 )
      CALL DCOPY( N-1, E, 1, WORK( 2 ), 2 )
      CALL DLASCL( 'G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1, IINFO )

      // Compute the q's and e's.

      DO 30 I = 1, 2*N - 1
         WORK( I ) = WORK( I )**2
   30 CONTINUE
      WORK( 2*N ) = ZERO

      CALL DLASQ2( N, WORK, INFO )

      if ( INFO.EQ.0 ) {
         DO 40 I = 1, N
            D( I ) = SQRT( WORK( I ) )
   40    CONTINUE
         CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
      } else if ( INFO.EQ.2 ) {

      // Maximum number of iterations exceeded.  Move data from WORK
      // into D and E so the calling subroutine can try to finish

         DO I = 1, N
            D( I ) = SQRT( WORK( 2*I-1 ) )
            E( I ) = SQRT( WORK( 2*I ) )
         END DO
         CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
         CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, E, N, IINFO )
      }

      RETURN

      // End of DLASQ1

      }
