      SUBROUTINE DORGTSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int               INFO, LDA, LDT, LWORK, M, N, MB, NB;
      // ..
      // .. Array Arguments ..
      double            A( LDA, * ), T( LDT, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                IINFO, LDC, LWORKOPT, LC, LW, NBLOCAL, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAMTSQR, DLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      LQUERY  = LWORK.EQ.-1
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. M.LT.N ) THEN
         INFO = -2
      ELSE IF( MB.LE.N ) THEN
         INFO = -3
      ELSE IF( NB.LT.1 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDT.LT.MAX( 1, MIN( NB, N ) ) ) THEN
         INFO = -8
      ELSE

         // Test the input LWORK for the dimension of the array WORK.
         // This workspace is used to store array C(LDC, N) and WORK(LWORK)
         // in the call to DLAMTSQR. See the documentation for DLAMTSQR.

         IF( LWORK.LT.2 .AND. (.NOT.LQUERY) ) THEN
            INFO = -10
         ELSE

            // Set block size for column blocks

            NBLOCAL = MIN( NB, N )

            // LWORK = -1, then set the size for the array C(LDC,N)
            // in DLAMTSQR call and set the optimal size of the work array
            // WORK(LWORK) in DLAMTSQR call.

            LDC = M
            LC = LDC*N
            LW = N * NBLOCAL

            LWORKOPT = LC+LW

            IF( ( LWORK.LT.MAX( 1, LWORKOPT ) ).AND.(.NOT.LQUERY) ) THEN
               INFO = -10
            END IF
         END IF

      END IF

      // Handle error in the input parameters and return workspace query.

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGTSQR', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
         WORK( 1 ) = DBLE( LWORKOPT )
         RETURN
      END IF

      // Quick return if possible

      IF( MIN( M, N ).EQ.0 ) THEN
         WORK( 1 ) = DBLE( LWORKOPT )
         RETURN
      END IF

      // (1) Form explicitly the tall-skinny M-by-N left submatrix Q1_in
      // of M-by-M orthogonal matrix Q_in, which is implicitly stored in
     t // he subdiagonal part of input array A and in the input array T.
      // Perform by the following operation using the routine DLAMTSQR.

          // Q1_in = Q_in * ( I ), where I is a N-by-N identity matrix,
                         // ( 0 )        0 is a (M-N)-by-N zero matrix.

      // (1a) Form M-by-N matrix in the array WORK(1:LDC*N) with ones
      // on the diagonal and zeros elsewhere.

      CALL DLASET( 'F', M, N, ZERO, ONE, WORK, LDC )

      // (1b)  On input, WORK(1:LDC*N) stores ( I );
                                           // ( 0 )

            // On output, WORK(1:LDC*N) stores Q1_in.

      CALL DLAMTSQR( 'L', 'N', M, N, N, MB, NBLOCAL, A, LDA, T, LDT, WORK, LDC, WORK( LC+1 ), LW, IINFO )

      // (2) Copy the result from the part of the work array (1:M,1:N)
      // with the leading dimension LDC that starts at WORK(1) into
     t // he output array A(1:M,1:N) column-by-column.

      DO J = 1, N
         CALL DCOPY( M, WORK( (J-1)*LDC + 1 ), 1, A( 1, J ), 1 )
      END DO

      WORK( 1 ) = DBLE( LWORKOPT )
      RETURN

      // End of DORGTSQR

      END
