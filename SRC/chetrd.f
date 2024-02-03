      SUBROUTINE CHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, NBMIN, NX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHER2K, CHETD2, CLATRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LWORK.LT.1 && .NOT.LQUERY ) {
         INFO = -9
      }

      if ( INFO == 0 ) {

         // Determine the block size.

         NB = ILAENV( 1, 'CHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO != 0 ) {
         xerbla('CHETRD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      NX = N
      IWS = 1
      if ( NB.GT.1 && NB.LT.N ) {

         // Determine when to cross over from blocked to unblocked code
         // (last block is always handled by unblocked code).

         NX = MAX( NB, ILAENV( 3, 'CHETRD', UPLO, N, -1, -1, -1 ) )
         if ( NX.LT.N ) {

            // Determine if workspace is large enough for blocked code.

            LDWORK = N
            IWS = LDWORK*NB
            if ( LWORK.LT.IWS ) {

               // Not enough workspace to use optimal NB:  determine the
               // minimum value of NB, and reduce NB or force use of
               // unblocked code by setting NX = N.

               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'CHETRD', UPLO, N, -1, -1, -1 )
               if (NB.LT.NBMIN) NX = N;
            }
         } else {
            NX = N
         }
      } else {
         NB = 1
      }

      if ( UPPER ) {

         // Reduce the upper triangle of A.
         // Columns 1:kk are handled by the unblocked method.

         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB

            // Reduce columns i:i+nb-1 to tridiagonal form and form the
            // matrix W which is needed to update the unreduced part of
            // the matrix

            clatrd(UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK, LDWORK );

            // Update the unreduced submatrix A(1:i-1,1:i-1), using an
            // update of the form:  A := A - V*W**H - W*V**H

            cher2k(UPLO, 'No transpose', I-1, NB, -CONE, A( 1, I ), LDA, WORK, LDWORK, ONE, A, LDA );

            // Copy superdiagonal elements back into A, and diagonal
            // elements into D

            for (J = I; J <= I + NB - 1; J++) { // 10
               A( J-1, J ) = E( J-1 )
               D( J ) = REAL( A( J, J ) )
            } // 10
         } // 20

         // Use unblocked code to reduce the last or only block

         chetd2(UPLO, KK, A, LDA, D, E, TAU, IINFO );
      } else {

         // Reduce the lower triangle of A

         DO 40 I = 1, N - NX, NB

            // Reduce columns i:i+nb-1 to tridiagonal form and form the
            // matrix W which is needed to update the unreduced part of
            // the matrix

            clatrd(UPLO, N-I+1, NB, A( I, I ), LDA, E( I ), TAU( I ), WORK, LDWORK );

            // Update the unreduced submatrix A(i+nb:n,i+nb:n), using
            // an update of the form:  A := A - V*W**H - W*V**H

            cher2k(UPLO, 'No transpose', N-I-NB+1, NB, -CONE, A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE, A( I+NB, I+NB ), LDA );

            // Copy subdiagonal elements back into A, and diagonal
            // elements into D

            for (J = I; J <= I + NB - 1; J++) { // 30
               A( J+1, J ) = E( J )
               D( J ) = REAL( A( J, J ) )
            } // 30
         } // 40

         // Use unblocked code to reduce the last or only block

         chetd2(UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ), TAU( I ), IINFO );
      }

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of CHETRD

      }
