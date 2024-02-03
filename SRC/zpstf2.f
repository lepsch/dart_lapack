      SUBROUTINE ZPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             TOL;
      int                INFO, LDA, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      double             WORK( 2*N );
      int                PIV( N );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         ZTEMP
      double             AJJ, DSTOP, DTEMP;
      int                I, ITEMP, J, PVT;
      bool               UPPER;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      bool               LSAME, DISNAN;
      // EXTERNAL DLAMCH, LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZGEMV, ZLACGV, ZSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCONJG, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZPSTF2', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Initialize PIV

      DO 100 I = 1, N
         PIV( I ) = I
  100 CONTINUE

      // Compute stopping value

      DO 110 I = 1, N
         WORK( I ) = DBLE( A( I, I ) )
  110 CONTINUE
      PVT = MAXLOC( WORK( 1:N ), 1 )
      AJJ = DBLE( A( PVT, PVT ) )
      if ( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) {
         RANK = 0
         INFO = 1
         GO TO 200
      }

      // Compute stopping value if not supplied

      if ( TOL.LT.ZERO ) {
         DSTOP = N * DLAMCH( 'Epsilon' ) * AJJ
      } else {
         DSTOP = TOL
      }

      // Set first half of WORK to zero, holds dot products

      DO 120 I = 1, N
         WORK( I ) = 0
  120 CONTINUE

      if ( UPPER ) {

         // Compute the Cholesky factorization P**T * A * P = U**H* U

         DO 150 J = 1, N

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            DO 130 I = J, N

               if ( J.GT.1 ) {
                  WORK( I ) = WORK( I ) + DBLE( DCONJG( A( J-1, I ) )* A( J-1, I ) )
               }
               WORK( N+I ) = DBLE( A( I, I ) ) - WORK( I )

  130       CONTINUE

            if ( J.GT.1 ) {
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               if ( AJJ.LE.DSTOP.OR.DISNAN( AJJ ) ) {
                  A( J, J ) = AJJ
                  GO TO 190
               }
            }

            if ( J.NE.PVT ) {

               // Pivot OK, so can now swap pivot rows and columns

               A( PVT, PVT ) = A( J, J )
               CALL ZSWAP( J-1, A( 1, J ), 1, A( 1, PVT ), 1 )
               IF( PVT.LT.N ) CALL ZSWAP( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA )
               DO 140 I = J + 1, PVT - 1
                  ZTEMP = DCONJG( A( J, I ) )
                  A( J, I ) = DCONJG( A( I, PVT ) )
                  A( I, PVT ) = ZTEMP
  140          CONTINUE
               A( J, PVT ) = DCONJG( A( J, PVT ) )

               // Swap dot products and PIV

               DTEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = DTEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            }

            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ

            // Compute elements J+1:N of row J

            if ( J.LT.N ) {
               CALL ZLACGV( J-1, A( 1, J ), 1 )
               CALL ZGEMV( 'Trans', J-1, N-J, -CONE, A( 1, J+1 ), LDA, A( 1, J ), 1, CONE, A( J, J+1 ), LDA )
               CALL ZLACGV( J-1, A( 1, J ), 1 )
               CALL ZDSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            }

  150    CONTINUE

      } else {

         // Compute the Cholesky factorization P**T * A * P = L * L**H

         DO 180 J = 1, N

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            DO 160 I = J, N

               if ( J.GT.1 ) {
                  WORK( I ) = WORK( I ) + DBLE( DCONJG( A( I, J-1 ) )* A( I, J-1 ) )
               }
               WORK( N+I ) = DBLE( A( I, I ) ) - WORK( I )

  160       CONTINUE

            if ( J.GT.1 ) {
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               if ( AJJ.LE.DSTOP.OR.DISNAN( AJJ ) ) {
                  A( J, J ) = AJJ
                  GO TO 190
               }
            }

            if ( J.NE.PVT ) {

               // Pivot OK, so can now swap pivot rows and columns

               A( PVT, PVT ) = A( J, J )
               CALL ZSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
               IF( PVT.LT.N ) CALL ZSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 )
               DO 170 I = J + 1, PVT - 1
                  ZTEMP = DCONJG( A( I, J ) )
                  A( I, J ) = DCONJG( A( PVT, I ) )
                  A( PVT, I ) = ZTEMP
  170          CONTINUE
               A( PVT, J ) = DCONJG( A( PVT, J ) )

               // Swap dot products and PIV

               DTEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = DTEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            }

            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ

            // Compute elements J+1:N of column J

            if ( J.LT.N ) {
               CALL ZLACGV( J-1, A( J, 1 ), LDA )
               CALL ZGEMV( 'No Trans', N-J, J-1, -CONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 )
               CALL ZLACGV( J-1, A( J, 1 ), LDA )
               CALL ZDSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            }

  180    CONTINUE

      }

      // Ran to completion, A has full rank

      RANK = N

      GO TO 200
  190 CONTINUE

      // Rank is number of steps completed.  Set INFO = 1 to signal
      // that the factorization cannot be used to solve a system.

      RANK = J - 1
      INFO = 1

  200 CONTINUE
      RETURN

      // End of ZPSTF2

      }
