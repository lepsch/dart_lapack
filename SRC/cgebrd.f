      SUBROUTINE CGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      const              ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IINFO, J, LDWRKX, LDWRKY, LWKMIN, LWKOPT, MINMN, NB, NBMIN, NX, WS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEBD2, CGEMM, CLABRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      MINMN = MIN( M, N )
      if ( MINMN == 0 ) {
         LWKMIN = 1
         LWKOPT = 1
      } else {
         LWKMIN = MAX( M, N )
         NB = MAX( 1, ILAENV( 1, 'CGEBRD', ' ', M, N, -1, -1 ) )
         LWKOPT = ( M+N )*NB
      }
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      LQUERY = ( LWORK == -1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -10
      }
      if ( INFO.LT.0 ) {
         xerbla('CGEBRD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      WS = MAX( M, N )
      LDWRKX = M
      LDWRKY = N

      if ( NB.GT.1 .AND. NB.LT.MINMN ) {

         // Set the crossover point NX.

         NX = MAX( NB, ILAENV( 3, 'CGEBRD', ' ', M, N, -1, -1 ) )

         // Determine when to switch from blocked to unblocked code.

         if ( NX.LT.MINMN ) {
            WS = LWKOPT
            if ( LWORK.LT.WS ) {

               // Not enough work space for the optimal NB, consider using
               // a smaller block size.

               NBMIN = ILAENV( 2, 'CGEBRD', ' ', M, N, -1, -1 )
               if ( LWORK.GE.( M+N )*NBMIN ) {
                  NB = LWORK / ( M+N )
               } else {
                  NB = 1
                  NX = MINMN
               }
            }
         }
      } else {
         NX = MINMN
      }

      DO 30 I = 1, MINMN - NX, NB

         // Reduce rows and columns i:i+ib-1 to bidiagonal form and return
         // the matrices X and Y which are needed to update the unreduced
         // part of the matrix

         clabrd(M-I+1, N-I+1, NB, A( I, I ), LDA, D( I ), E( I ), TAUQ( I ), TAUP( I ), WORK, LDWRKX, WORK( LDWRKX*NB+1 ), LDWRKY );

         // Update the trailing submatrix A(i+ib:m,i+ib:n), using
         // an update of the form  A := A - V*Y**H - X*U**H

         cgemm('No transpose', 'Conjugate transpose', M-I-NB+1, N-I-NB+1, NB, -ONE, A( I+NB, I ), LDA, WORK( LDWRKX*NB+NB+1 ), LDWRKY, ONE, A( I+NB, I+NB ), LDA );
         cgemm('No transpose', 'No transpose', M-I-NB+1, N-I-NB+1, NB, -ONE, WORK( NB+1 ), LDWRKX, A( I, I+NB ), LDA, ONE, A( I+NB, I+NB ), LDA );

         // Copy diagonal and off-diagonal elements of B back into A

         if ( M.GE.N ) {
            for (J = I; J <= I + NB - 1; J++) { // 10
               A( J, J ) = D( J )
               A( J, J+1 ) = E( J )
            } // 10
         } else {
            for (J = I; J <= I + NB - 1; J++) { // 20
               A( J, J ) = D( J )
               A( J+1, J ) = E( J )
            } // 20
         }
      } // 30

      // Use unblocked code to reduce the remainder of the matrix

      cgebd2(M-I+1, N-I+1, A( I, I ), LDA, D( I ), E( I ), TAUQ( I ), TAUP( I ), WORK, IINFO );
      WORK( 1 ) = SROUNDUP_LWORK( WS )
      RETURN

      // End of CGEBRD

      }
