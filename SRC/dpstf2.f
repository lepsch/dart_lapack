      SUBROUTINE DPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             TOL;
      int                INFO, LDA, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), WORK( 2*N );
      int                PIV( N );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
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
      // EXTERNAL DGEMV, DSCAL, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT, MAXLOC
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
         CALL XERBLA( 'DPSTF2', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Initialize PIV

      DO 100 I = 1, N
         PIV( I ) = I
  100 CONTINUE

      // Compute stopping value

      PVT = 1
      AJJ = A( PVT, PVT )
      DO I = 2, N
         if ( A( I, I ).GT.AJJ ) {
            PVT = I
            AJJ = A( PVT, PVT )
         }
      END DO
      if ( AJJ.LE.ZERO.OR.DISNAN( AJJ ) ) {
         RANK = 0
         INFO = 1
         GO TO 170
      }

      // Compute stopping value if not supplied

      if ( TOL.LT.ZERO ) {
         DSTOP = N * DLAMCH( 'Epsilon' ) * AJJ
      } else {
         DSTOP = TOL
      }

      // Set first half of WORK to zero, holds dot products

      DO 110 I = 1, N
         WORK( I ) = 0
  110 CONTINUE

      if ( UPPER ) {

         // Compute the Cholesky factorization P**T * A * P = U**T * U

         DO 130 J = 1, N

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            DO 120 I = J, N

               if ( J.GT.1 ) {
                  WORK( I ) = WORK( I ) + A( J-1, I )**2
               }
               WORK( N+I ) = A( I, I ) - WORK( I )

  120       CONTINUE

            if ( J.GT.1 ) {
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               if ( AJJ.LE.DSTOP.OR.DISNAN( AJJ ) ) {
                  A( J, J ) = AJJ
                  GO TO 160
               }
            }

            if ( J.NE.PVT ) {

               // Pivot OK, so can now swap pivot rows and columns

               A( PVT, PVT ) = A( J, J )
               CALL DSWAP( J-1, A( 1, J ), 1, A( 1, PVT ), 1 )
               IF( PVT.LT.N ) CALL DSWAP( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA )
               CALL DSWAP( PVT-J-1, A( J, J+1 ), LDA, A( J+1, PVT ), 1 )

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
               CALL DGEMV( 'Trans', J-1, N-J, -ONE, A( 1, J+1 ), LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            }

  130    CONTINUE

      } else {

         // Compute the Cholesky factorization P**T * A * P = L * L**T

         DO 150 J = 1, N

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            DO 140 I = J, N

               if ( J.GT.1 ) {
                  WORK( I ) = WORK( I ) + A( I, J-1 )**2
               }
               WORK( N+I ) = A( I, I ) - WORK( I )

  140       CONTINUE

            if ( J.GT.1 ) {
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               if ( AJJ.LE.DSTOP.OR.DISNAN( AJJ ) ) {
                  A( J, J ) = AJJ
                  GO TO 160
               }
            }

            if ( J.NE.PVT ) {

               // Pivot OK, so can now swap pivot rows and columns

               A( PVT, PVT ) = A( J, J )
               CALL DSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
               IF( PVT.LT.N ) CALL DSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 )
               CALL DSWAP( PVT-J-1, A( J+1, J ), 1, A( PVT, J+1 ), LDA )

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
               CALL DGEMV( 'No Trans', N-J, J-1, -ONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            }

  150    CONTINUE

      }

      // Ran to completion, A has full rank

      RANK = N

      GO TO 170
  160 CONTINUE

      // Rank is number of steps completed.  Set INFO = 1 to signal
      // that the factorization cannot be used to solve a system.

      RANK = J - 1
      INFO = 1

  170 CONTINUE
      RETURN

      // End of DPSTF2

      }
