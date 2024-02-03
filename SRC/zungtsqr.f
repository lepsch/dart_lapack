      SUBROUTINE ZUNGTSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int               INFO, LDA, LDT, LWORK, M, N, MB, NB;
      // ..
      // .. Array Arguments ..
      COMPLEX*16        A( LDA, * ), T( LDT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CONE, CZERO
      const              CONE = ( 1.0D+0, 0.0D+0 ), CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                IINFO, LDC, LWORKOPT, LC, LW, NBLOCAL, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZLAMTSQR, ZLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      LQUERY  = LWORK == -1
      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 .OR. M.LT.N ) {
         INFO = -2
      } else if ( MB.LE.N ) {
         INFO = -3
      } else if ( NB.LT.1 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -6
      } else if ( LDT.LT.MAX( 1, MIN( NB, N ) ) ) {
         INFO = -8
      } else {

         // Test the input LWORK for the dimension of the array WORK.
         // This workspace is used to store array C(LDC, N) and WORK(LWORK)
         // in the call to ZLAMTSQR. See the documentation for ZLAMTSQR.

         if ( LWORK.LT.2 && (.NOT.LQUERY) ) {
            INFO = -10
         } else {

            // Set block size for column blocks

            NBLOCAL = MIN( NB, N )

            // LWORK = -1, then set the size for the array C(LDC,N)
            // in ZLAMTSQR call and set the optimal size of the work array
            // WORK(LWORK) in ZLAMTSQR call.

            LDC = M
            LC = LDC*N
            LW = N * NBLOCAL

            LWORKOPT = LC+LW

            if ( ( LWORK.LT.MAX( 1, LWORKOPT ) ) && (.NOT.LQUERY) ) {
               INFO = -10
            }
         }

      }

      // Handle error in the input parameters and return workspace query.

      if ( INFO != 0 ) {
         xerbla('ZUNGTSQR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = DCMPLX( LWORKOPT )
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N ) == 0 ) {
         WORK( 1 ) = DCMPLX( LWORKOPT )
         RETURN
      }

      // (1) Form explicitly the tall-skinny M-by-N left submatrix Q1_in
      // of M-by-M orthogonal matrix Q_in, which is implicitly stored in
      // the subdiagonal part of input array A and in the input array T.
      // Perform by the following operation using the routine ZLAMTSQR.

          // Q1_in = Q_in * ( I ), where I is a N-by-N identity matrix,
                         // ( 0 )        0 is a (M-N)-by-N zero matrix.

      // (1a) Form M-by-N matrix in the array WORK(1:LDC*N) with ones
      // on the diagonal and zeros elsewhere.

      zlaset('F', M, N, CZERO, CONE, WORK, LDC );

      // (1b)  On input, WORK(1:LDC*N) stores ( I );
                                           // ( 0 )

            // On output, WORK(1:LDC*N) stores Q1_in.

      zlamtsqr('L', 'N', M, N, N, MB, NBLOCAL, A, LDA, T, LDT, WORK, LDC, WORK( LC+1 ), LW, IINFO );

      // (2) Copy the result from the part of the work array (1:M,1:N)
      // with the leading dimension LDC that starts at WORK(1) into
      // the output array A(1:M,1:N) column-by-column.

      for (J = 1; J <= N; J++) {
         zcopy(M, WORK( (J-1)*LDC + 1 ), 1, A( 1, J ), 1 );
      }

      WORK( 1 ) = DCMPLX( LWORKOPT )
      RETURN

      // End of ZUNGTSQR

      }
