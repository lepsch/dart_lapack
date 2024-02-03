      SUBROUTINE DPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )

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
      int                I, ITEMP, J, JB, K, NB, PVT;
      bool               UPPER;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      int                ILAENV;
      bool               LSAME, DISNAN;
      // EXTERNAL DLAMCH, ILAENV, LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DPSTF2, DSCAL, DSWAP, DSYRK, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT, MAXLOC
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

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
         CALL XERBLA( 'DPSTRF', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Get block size

      NB = ILAENV( 1, 'DPOTRF', UPLO, N, -1, -1, -1 )
      if ( NB.LE.1 .OR. NB.GE.N ) {

         // Use unblocked code

         CALL DPSTF2( UPLO, N, A( 1, 1 ), LDA, PIV, RANK, TOL, WORK, INFO )
         GO TO 200

      } else {

      // Initialize PIV

         DO 100 I = 1, N
            PIV( I ) = I
  100    CONTINUE

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
            GO TO 200
         }

      // Compute stopping value if not supplied

         if ( TOL.LT.ZERO ) {
            DSTOP = N * DLAMCH( 'Epsilon' ) * AJJ
         } else {
            DSTOP = TOL
         }


         if ( UPPER ) {

            // Compute the Cholesky factorization P**T * A * P = U**T * U

            DO 140 K = 1, N, NB

               // Account for last block not being NB wide

               JB = MIN( NB, N-K+1 )

               // Set relevant part of first half of WORK to zero,
               // holds dot products

               DO 110 I = K, N
                  WORK( I ) = 0
  110          CONTINUE

               DO 130 J = K, K + JB - 1

               // Find pivot, test for exit, else swap rows and columns
               // Update dot products, compute possible pivots which are
               // stored in the second half of WORK

                  DO 120 I = J, N

                     if ( J.GT.K ) {
                        WORK( I ) = WORK( I ) + A( J-1, I )**2
                     }
                     WORK( N+I ) = A( I, I ) - WORK( I )

  120             CONTINUE

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

                  // Compute elements J+1:N of row J.

                  if ( J.LT.N ) {
                     CALL DGEMV( 'Trans', J-K, N-J, -ONE, A( K, J+1 ), LDA, A( K, J ), 1, ONE, A( J, J+1 ), LDA )
                     CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
                  }

  130          CONTINUE

               // Update trailing matrix, J already incremented

               if ( K+JB.LE.N ) {
                  CALL DSYRK( 'Upper', 'Trans', N-J+1, JB, -ONE, A( K, J ), LDA, ONE, A( J, J ), LDA )
               }

  140       CONTINUE

         } else {

         // Compute the Cholesky factorization P**T * A * P = L * L**T

            DO 180 K = 1, N, NB

               // Account for last block not being NB wide

               JB = MIN( NB, N-K+1 )

               // Set relevant part of first half of WORK to zero,
               // holds dot products

               DO 150 I = K, N
                  WORK( I ) = 0
  150          CONTINUE

               DO 170 J = K, K + JB - 1

               // Find pivot, test for exit, else swap rows and columns
               // Update dot products, compute possible pivots which are
               // stored in the second half of WORK

                  DO 160 I = J, N

                     if ( J.GT.K ) {
                        WORK( I ) = WORK( I ) + A( I, J-1 )**2
                     }
                     WORK( N+I ) = A( I, I ) - WORK( I )

  160             CONTINUE

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

                  // Compute elements J+1:N of column J.

                  if ( J.LT.N ) {
                     CALL DGEMV( 'No Trans', N-J, J-K, -ONE, A( J+1, K ), LDA, A( J, K ), LDA, ONE, A( J+1, J ), 1 )
                     CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
                  }

  170          CONTINUE

               // Update trailing matrix, J already incremented

               if ( K+JB.LE.N ) {
                  CALL DSYRK( 'Lower', 'No Trans', N-J+1, JB, -ONE, A( J, K ), LDA, ONE, A( J, J ), LDA )
               }

  180       CONTINUE

         }
      }

      // Ran to completion, A has full rank

      RANK = N

      GO TO 200
  190 CONTINUE

      // Rank is the number of steps completed.  Set INFO = 1 to signal
     t // hat the factorization cannot be used to solve a system.

      RANK = J - 1
      INFO = 1

  200 CONTINUE
      RETURN

      // End of DPSTRF

      }
